#pragma once

#include <cmath>
#include <type_traits>
#include <cassert>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <random>

// Template class for vectors.
template<typename TScalarType, size_t TDim>
class Vec
{
public:
	using TScalar = TScalarType;
	static constexpr size_t Dim = TDim;
	Vec() { for (size_t i=0; i<Dim; ++i) mData[i] = 0; }
	template<class TOther>
	Vec(const TOther& other) { 
		for (size_t i = 0; i < std::min(TDim, TOther::Dim); ++i) mData[i] = static_cast<TScalarType>(other[i]); 
		for (size_t i = TOther::Dim; i < TDim; ++i) mData[i] = 0;
	}
	explicit Vec(std::initializer_list<TScalar> list) {
		size_t i = 0;
		for (auto elem : list) { if (i >= TDim) return; mData[i++] = elem; }
		for (; i < TDim; ++i) mData[i] = 0;
	}
	Vec operator+(const Vec& other) const { Vec result; for (size_t i = 0; i < TDim; ++i) result[i] = mData[i] + other[i]; return result; }
	Vec operator-(const Vec& other) const { Vec result; for (size_t i = 0; i < TDim; ++i) result[i] = mData[i] - other[i]; return result; }
	Vec operator*(const Vec& other) const { Vec result; for (size_t i = 0; i < TDim; ++i) result[i] = mData[i] * other[i]; return result; }
	Vec operator/(const Vec& other) const { Vec result; for (size_t i = 0; i < TDim; ++i) result[i] = mData[i] / other[i]; return result; }
	Vec operator*(const TScalar& other) const { Vec result; for (size_t i = 0; i < TDim; ++i) result[i] = mData[i] * other; return result; }
	Vec operator/(const TScalar& other) const { Vec result; for (size_t i = 0; i < TDim; ++i) result[i] = mData[i] / other; return result; }
	bool operator==(const Vec& other) const { for (size_t i = 0; i < TDim; ++i) if (mData[i] != other[i]) return false; return true; }
	Vec operator-() const { Vec result; for (size_t i = 0; i < TDim; ++i) result[i] = -mData[i]; return result; }
	Vec& operator+=(const Vec& other) { for (size_t i = 0; i < TDim; ++i) mData[i] += other[i]; return *this; }
	Vec& operator-=(const Vec& other) { for (size_t i = 0; i < TDim; ++i) mData[i] -= other[i]; return *this; }
	Vec& operator*=(const TScalar& other) { for (size_t i = 0; i < TDim; ++i) mData[i] *= other; return *this; }
	Vec& operator/=(const TScalar& other) { for (size_t i = 0; i < TDim; ++i) mData[i] /= other; return *this; }
	TScalar& operator[](size_t i) { assert(i < TDim); if (i < TDim) return mData[i]; else return mData[0]; }
	const TScalar& operator[](size_t i) const { assert(i < TDim); if (i < TDim) return mData[i]; else return mData[0]; }
	TScalar lengthSquared() const { TScalar result(0); for (size_t i = 0; i < TDim; ++i) result += mData[i] * mData[i]; return result; }
	static Vec ones() { Vec result; for (size_t i = 0; i < TDim; ++i) result[i] = TScalar(1); return result; }
protected:
	TScalar mData[TDim];
};
template<typename TScalarType, size_t TDim>
inline const Vec<TScalarType, TDim> operator*(const TScalarType& val, const Vec<TScalarType, TDim>& v) { return v * static_cast<TScalarType>(val); }
typedef Vec<float, 2> Vec2f;
typedef Vec<float, 3> Vec3f; 
typedef Vec<double, 2> Vec2d;
typedef Vec<double, 3> Vec3d;
typedef Vec<int, 2> Vec2i;
typedef Vec<int, 3> Vec3i;

// Bounding Box in domain coordinates.
template<typename TVector>
class BoundingBox {
public:
	BoundingBox(const TVector& minimum, const TVector& maximum) : Min(minimum), Max(maximum) {}

	// Checks if a point is inside the box.
	bool Contains(const TVector& x) const {
		for (size_t i = 0; i < TVector::Dim; ++i)
			if (x[i] < Min[i] || Max[i] < x[i])
				return false;
		return true;
	}

	// Clamps values outside of the bounding box back on to the bounding box
	TVector ClampToDomain(const TVector& pnt) const {
		TVector result = pnt;
		for (int i = 0; i < TVector::Dim; ++i)
			result[i] = std::min(std::max(Min[i], pnt[i]), Max[i]);
		return result;
	}
	TVector Min, Max;
};
typedef BoundingBox<Vec2d> BoundingBox2d;
typedef BoundingBox<Vec3d> BoundingBox3d;

// Stores a field on a regular grid.
template<typename TValueType, size_t TDim>
class RegularGrid
{
public:
	using TGridCoord = Vec<int, TDim>;
	static constexpr size_t Dim = TDim;
	using TValue = TValueType;
	using TBoundingBox = BoundingBox<Vec<double, TDim>>;
	using TDomainCoord = Vec<double, TDim>;

	// Allocates own data.
	RegularGrid(const TGridCoord& res, const TBoundingBox& domain) : mResolution(res), mDomain(domain), mExternData(NULL) { 
		int numElements = 1;
		for (int d = 0; d < this->Dim; ++d) numElements *= res[d];
		mData.resize(numElements); 
	}
	// Uses external data without taking ownership.
	RegularGrid(const TGridCoord& res, const TBoundingBox& domain, TValue* data) : mResolution(res), mDomain(domain), mExternData(data) {}
	RegularGrid(const RegularGrid&) = delete;
	virtual ~RegularGrid() {}

	// Linearly samples the field.
	TValue Sample(const TDomainCoord& coord) const
	{
		TDomainCoord position = this->mDomain.ClampToDomain(coord);
		TDomainCoord vfTex = (position - mDomain.Min) / (mDomain.Max - mDomain.Min);
		TDomainCoord vfSample = vfTex * (static_cast<TDomainCoord>(mResolution) - TDomainCoord::ones());
		TGridCoord viSampleBase0, viSampleBase1;
		for (int i = 0; i < TDomainCoord::Dim; ++i) {
			viSampleBase0[i] = std::min(std::max(0, (int)vfSample[i]), mResolution[i] - 1);
			viSampleBase1[i] = std::min(std::max(0, viSampleBase0[i] + 1), mResolution[i] - 1);
		}
		TDomainCoord vfSampleInterpol = vfSample - static_cast<TDomainCoord>(viSampleBase0);
		size_t numCorners = (size_t)std::pow(2, TDomainCoord::Dim);
		TValue result({ 0 });
		for (size_t i = 0; i < numCorners; ++i) {
			typename TDomainCoord::TScalar weight(1);
			TGridCoord gridCoord;
			for (size_t d = 0; d < TDomainCoord::Dim; ++d) {
				if (i & (size_t(1) << (TDomainCoord::Dim - size_t(1) - d))) {
					gridCoord[d] = viSampleBase1[d];
					weight *= vfSampleInterpol[d];
				}
				else {
					gridCoord[d] = viSampleBase0[d];
					weight *= 1 - vfSampleInterpol[d];
				}
			}
			result += static_cast<TValueType>(GetVertexDataAt(gridCoord)) * float(weight);
		}
		return result;
	}

	// Gets the linear array index based on a grid coordinate index.
	int GetLinearIndex(const TGridCoord& gridCoord) const {
		int stride = 1;
		int linearIndex = gridCoord[0];
		for (int d = 1; d < this->Dim; ++d) {
			stride *= mResolution[(int64_t)d - 1];
			linearIndex += gridCoord[d] * stride;
		}
		return linearIndex;
	}

	// Gets the spatial location of a grid vertex.
	TDomainCoord GetCoordAt(const TGridCoord& gridCoord) const {
		TDomainCoord s;
		for (int i = 0; i < TDomainCoord::Dim; ++i)
			s[i] = mResolution[i] < 2 ? 0.5 : gridCoord[i] / (mResolution[i] - typename TDomainCoord::TScalar(1.));
		return mDomain.Min + s * (mDomain.Max - mDomain.Min);
	}

	// Gets the grid coordinate based on the linear array index.
	TGridCoord GetGridCoord(const size_t& linearIndex) const {
		TGridCoord result;
		int stride = 1;
		for (int d = 0; d < TGridCoord::Dim - 1; ++d)
			stride *= mResolution[d];
		int t = (int)linearIndex;
		for (int d = TGridCoord::Dim - 1; d >= 0; --d) {
			result[d] = t / stride;
			t = t % stride;
			if (d > 0)
				stride /= mResolution[d - 1ll];
		}
		return result;
	}

	// Gets the grid coordinate of the voxel that contains the domain coordinate.
	TGridCoord GetGridCoord(const TDomainCoord& domainCoord) const {
		TDomainCoord relativeCoord = (domainCoord - mDomain.Min) / (mDomain.Max - mDomain.Min);
		TGridCoord result = relativeCoord * mResolution;
		for (int d = 0; d < TGridCoord::Dim; ++d)
			result[d] = std::min(std::max(0, result[d]), mResolution[d] - 1);
		return result;
	}

	// Gets the vertex data stored at a certain grid coordinate.
	TValue GetVertexDataAt(const TGridCoord& gridCoord) const {
		for (int dim = 0; dim < TDim; ++dim)
			assert(gridCoord[dim] >= 0 && gridCoord[dim] < mResolution[dim]);
		int addr = GetLinearIndex(gridCoord);
		return mExternData ? mExternData[addr] : mData[addr];
	}

	// Sets the vertex data at a certain grid coordinate.
	void SetVertexDataAt(const TGridCoord& gridCoord, const TValue& value) {
		for (int dim = 0; dim < TDim; ++dim)
			assert(gridCoord[dim] >= 0 && gridCoord[dim] < mResolution[dim]);
		int addr = GetLinearIndex(gridCoord);
		if (mExternData) mExternData[addr] = value; else mData[addr] = value;
	}

	const TValue* GetData() const { return mExternData ? mExternData : mData.data(); }
	const TGridCoord& GetResolution() const { return mResolution; }
	const TBoundingBox& GetDomain() const { return mDomain; }
	TDomainCoord GetVoxelSize() const { return (mDomain.Max - mDomain.Min) / (static_cast<TDomainCoord>(mResolution) - TDomainCoord::ones()); }
protected:
	TValue* mExternData;
	std::vector<TValue> mData;
	TGridCoord mResolution;
	TBoundingBox mDomain;
};
typedef RegularGrid<float, 2> SteadyScalarField2d;
typedef RegularGrid<Vec2f, 2> SteadyVectorField2d;
typedef RegularGrid<float, 3> UnsteadyScalarField2d;
typedef RegularGrid<Vec2f, 3> UnsteadyVectorField2d;

// Particle tracer.
class Tracer
{
public:
	// Get the position reached after a certain integration duration. Points are specified in space-time.
	template<typename TVectorField>
	static typename TVectorField::TDomainCoord FlowMap(const TVectorField& vectorField, typename TVectorField::TDomainCoord state, double stepSize, double duration, bool& indomain) {
		return Integrate(vectorField, state, stepSize, duration, EIntegralCurve::Pathline, indomain);
	}

	// Trace a streamline. Points are specified in space-time.
	template<typename TVectorField>
	static void Streamline(const TVectorField& vectorField, typename TVectorField::TDomainCoord state, double stepSize, double duration, int numSteps, std::vector<typename TVectorField::TDomainCoord>& curve) {
		TangentCurve(vectorField, state, stepSize, duration, numSteps, EIntegralCurve::Streamline, curve);
	}

	// Trace a pathline. Points are specified in space-time.
	template<typename TVectorField>
	static void Pathline(const TVectorField& vectorField, typename TVectorField::TDomainCoord state, double stepSize, double duration, int numSteps, std::vector<typename TVectorField::TDomainCoord>& curve) {
		TangentCurve(vectorField, state, stepSize, duration, numSteps, EIntegralCurve::Pathline, curve);
	}

	// Trace a streakline without refinement. Points are specified in space-time.
	template<typename TVectorField>
	static void Streakline(const TVectorField& vectorField, const typename TVectorField::TDomainCoord& state, double stepSize, double duration, int numSteps, std::vector<typename TVectorField::TDomainCoord>& curve) {
		Vec3d seed = state;
		curve.clear();
		curve.reserve(numSteps);
		for (int i = 0; i < numSteps; ++i) {
			Vec3d particle({ seed[0], seed[1], seed[2] + (stepSize > 0 ? 1 : -1) * duration * (1.0 - i / (numSteps - 1.0)) });
			bool indomain = true;
			typename TVectorField::TDomainCoord reached = FlowMap(vectorField, particle, stepSize, duration + (stepSize > 0 ? 1 : -1) * (seed[2] - particle[2]), indomain);
			if (indomain) curve.push_back(reached);
		}
	}

	// Trace a timeline without refinement. Points are specified in space-time.
	template<typename TVectorField>
	static void Timeline(const TVectorField& vectorField, const typename TVectorField::TDomainCoord& stateA, const typename TVectorField::TDomainCoord& stateB, double stepSize, double duration, int numSteps, std::vector<typename TVectorField::TDomainCoord>& curve) {
		curve.clear();
		curve.reserve(numSteps);
		for (int i = 0; i < numSteps; ++i) {
			typename TVectorField::TDomainCoord seed = stateA + (stateB - stateA) * i / (numSteps - 1.0);
			bool indomain = true;
			typename TVectorField::TDomainCoord reached = FlowMap(vectorField, seed, stepSize, duration, indomain);
			if (indomain) curve.push_back(reached);
		}
	}

private:
	Tracer() = delete;

	enum class EIntegralCurve {
		Streamline,
		Pathline
	};

	// Sample the space-time lifted field.
	template<typename TVectorField>
	static typename TVectorField::TDomainCoord Sample(const TVectorField& vectorField, const typename TVectorField::TDomainCoord& state, EIntegralCurve integralCurve, bool& indomain) {
		indomain = vectorField.GetDomain().Contains(state); if (!indomain) return typename TVectorField::TDomainCoord();
		typename TVectorField::TValue velocity = vectorField.Sample(state);
		if (integralCurve == EIntegralCurve::Streamline) {
			double norm = std::sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1]);
			return typename TVectorField::TDomainCoord({ norm != 0 ? velocity[0] / norm : 0, norm != 0 ? velocity[1] / norm : 0 });
		} else return typename TVectorField::TDomainCoord({ velocity[0], velocity[1], 1.0 });
	}

	// Perform fourth-order Runge-Kutta integration step.
	template<typename TVectorField>
	static typename TVectorField::TDomainCoord Step(const TVectorField& vectorField, const typename TVectorField::TDomainCoord& state, double s, EIntegralCurve integralCurve, bool& indomain) {
		auto v1 = Sample(vectorField, state, integralCurve, indomain);					if (!indomain) return typename TVectorField::TDomainCoord();
		auto v2 = Sample(vectorField, state + v1 * s / 2, integralCurve, indomain);		if (!indomain) return typename TVectorField::TDomainCoord();
		auto v3 = Sample(vectorField, state + v2 * s / 2, integralCurve, indomain);		if (!indomain) return typename TVectorField::TDomainCoord();
		auto v4 = Sample(vectorField, state + v3 * s, integralCurve, indomain);			if (!indomain) return typename TVectorField::TDomainCoord();
		return state + (v1 / 6 + v2 / 3 + v3 / 3 + v4 / 6) * s;
	}

	// Integrates a curve segment in the vector field.
	template<typename TVectorField>
	static typename TVectorField::TDomainCoord Integrate(const TVectorField& vectorField, typename TVectorField::TDomainCoord state, double stepSize, double duration, EIntegralCurve integralCurve, bool& indomain)
	{
		double tau = 0;
		while ((stepSize > 0 ? 1 : -1) * tau < duration && indomain) {
			double dt = (stepSize > 0 ? 1 : -1) ? std::min(stepSize, duration - tau) : std::max(stepSize, tau - duration);
			state = Step(vectorField, state, dt, integralCurve, indomain);
			tau += dt;
		}
		return state;
	}

	// Trace tangent curve with equidistant sampling in time.
	template<typename TVectorField>
	static void TangentCurve(const TVectorField& vectorField, typename TVectorField::TDomainCoord state, double stepSize, double duration, int numSteps, EIntegralCurve integralCurve, std::vector<typename TVectorField::TDomainCoord>& curve) {
		curve.clear();
		curve.reserve(numSteps);
		curve.push_back(state);
		bool indomain = true;
		for (int i = 0; (i < numSteps - 1) && indomain; ++i) {
			state = Integrate(vectorField, state, stepSize, duration / (numSteps - 1.0), integralCurve, indomain);
			if (indomain) curve.push_back(state);
		}
	}
};

// Computation of vorticity for a 2D vector field.
class Vorticity2d
{
public:
	static void Compute(const SteadyVectorField2d& flow, SteadyScalarField2d& result) { ComputeSlice(flow, result, 0); }
	static void Compute(const UnsteadyVectorField2d& flow, SteadyScalarField2d& result, int slice) { ComputeSlice(flow, result, slice); }
	static void Compute(const UnsteadyVectorField2d& flow, UnsteadyScalarField2d& result) { for (int slice = 0; slice < result.GetResolution()[2]; ++slice) ComputeSlice(flow, result, slice); }
private:
	friend class InstantaneousVorticityDeviation2d;
	template<size_t DimIn, size_t DimOut>
	static void ComputeSlice(const RegularGrid<Vec2f, DimIn>& flow, RegularGrid<float, DimOut>& result, int slice) {
		for (int dim = 0; dim < std::min(DimIn, DimOut); ++dim) {
			assert(flow.GetResolution()[dim] == result.GetResolution()[dim]);
			assert(flow.GetDomain().Min[dim] == result.GetDomain().Min[dim]);
			assert(flow.GetDomain().Max[dim] == result.GetDomain().Max[dim]);
		}
		Vec2i res = result.GetResolution();
		#pragma omp parallel for
		for (int i = 0; i < res[0] * res[1]; ++i) {
			Vec3i gridCoord = flow.GetGridCoord(i + size_t(slice) * res[0] * res[1]);
			int ix = gridCoord[0], iy = gridCoord[1], ix0 = std::max(0, ix - 1), ix1 = std::min(ix + 1, res[0] - 1), iy0 = std::max(0, iy - 1), iy1 = std::min(iy + 1, res[1] - 1);
			double uy = ((double)flow.GetVertexDataAt(Vec3i({ ix, iy1, slice }))[0] - flow.GetVertexDataAt(Vec3i({ ix, iy0, slice }))[0]) / (((double)iy1 - iy0) * flow.GetVoxelSize()[1]);
			double vx = ((double)flow.GetVertexDataAt(Vec3i({ ix1, iy, slice }))[1] - flow.GetVertexDataAt(Vec3i({ ix0, iy, slice }))[1]) / (((double)ix1 - ix0) * flow.GetVoxelSize()[0]);
			result.SetVertexDataAt(gridCoord, vx - uy);
		}
	}
};

// Computation of the finite-time Lyapunov exponent (FTLE) for a 2D vector field.
// Shadden, S. C., Lekien, F., & Marsden, J. E. (2005). Definition and properties of Lagrangian coherent structures from finite-time Lyapunov exponents in two-dimensional aperiodic flows. Physica D: Nonlinear Phenomena, 212(3-4), 271-304.
class FiniteTimeLyapunovExponent2d
{
public:
	static void Compute(const SteadyVectorField2d& flow, double stepSize, double duration, SteadyScalarField2d& result) { ComputeSlice(flow, stepSize, duration, result, 0); }
	static void Compute(const UnsteadyVectorField2d& flow, double stepSize, double duration, SteadyScalarField2d& result, int slice) { ComputeSlice(flow, stepSize, duration, result, slice); }
	static void Compute(const UnsteadyVectorField2d& flow, double stepSize, double duration, UnsteadyScalarField2d& result) { for (int slice = 0; slice < result.GetResolution()[2]; ++slice) ComputeSlice(flow, stepSize, duration, result, slice); }
private:
	template<size_t DimIn, size_t DimOut>
	static void ComputeSlice(const RegularGrid<Vec2f, DimIn>& flow, double stepSize, double duration, RegularGrid<float, DimOut>& result, int slice) {
		Vec2i res = result.GetResolution();
		SteadyVectorField2d phi(res, BoundingBox2d(Vec2d({ result.GetDomain().Min[0], result.GetDomain().Min[1]}), Vec2d({ result.GetDomain().Max[0], result.GetDomain().Max[1] })));
		RegularGrid<bool, 2> valid(res, phi.GetDomain());
		#pragma omp parallel for schedule(guided)
		for (int i = 0; i < res[0] * res[1]; ++i) {
			Vec2i gridCoord = phi.GetGridCoord(i);
			Vec2d coord2 = phi.GetCoordAt(gridCoord);
			Vec3d coord3 = flow.GetCoordAt(Vec3i({ 0,0,slice }));
			Vec3d coord({ coord2[0], coord2[1], coord3[2] });
			bool indomain = true;
			phi.SetVertexDataAt(gridCoord, Tracer::FlowMap(flow, coord, stepSize, duration, indomain));
			valid.SetVertexDataAt(gridCoord, indomain);
		}
		#pragma omp parallel for schedule(guided)
		for (int i = 0; i < res[0] * res[1]; ++i) {
			Vec2i gridCoord2 = phi.GetGridCoord(i);
			Vec3i gridCoord3({ gridCoord2[0], gridCoord2[1], slice });
			int ix = gridCoord2[0], iy = gridCoord2[1], ix0 = std::max(0, ix - 1), ix1 = std::min(ix + 1, res[0] - 1), iy0 = std::max(0, iy - 1), iy1 = std::min(iy + 1, res[1] - 1);
			if (valid.GetVertexDataAt(Vec2i({ ix1, iy })) && valid.GetVertexDataAt(Vec2i({ ix0, iy })) && valid.GetVertexDataAt(Vec2i({ ix, iy1 })) && valid.GetVertexDataAt(Vec2i({ ix, iy0 })) && valid.GetVertexDataAt(Vec2i({ ix, iy }))) {
				Vec2d phix = (phi.GetVertexDataAt(Vec2i({ ix1, iy })) - phi.GetVertexDataAt(Vec2i({ ix0, iy }))) / (((double)ix1 - ix0) * phi.GetVoxelSize()[0]);
				Vec2d phiy = (phi.GetVertexDataAt(Vec2i({ ix, iy1 })) - phi.GetVertexDataAt(Vec2i({ ix, iy0 }))) / (((double)iy1 - iy0) * phi.GetVoxelSize()[1]);
				double a = phix[0], b = phiy[0], c = phix[1], d = phiy[1];
				double lambdaMax = (std::sqrt(d * d * d * d + (2 * c * c + 2 * b * b - 2 * a * a) * d * d + 8 * a * b * c * d + c * c * c * c + (2 * a * a - 2 * b * b) * c * c + b * b * b * b + 2 * a * a * b * b + a * a * a * a) + d * d + c * c + b * b + a * a) / 2;
				if (lambdaMax > 0)
					result.SetVertexDataAt(gridCoord3, 1.0 / std::abs(duration) * std::log(std::sqrt(lambdaMax)));
				else result.SetVertexDataAt(gridCoord3, 0);
			} else result.SetVertexDataAt(gridCoord3, 0);
		}
	}
};

// Computation of the instantaneous vorticity deviation (IVD) for a 2D vector field.
// Haller, G., Hadjighasem, A., Farazmand, M., & Huhn, F. (2016). Defining coherent vortices objectively from the vorticity. Journal of Fluid Mechanics, 795, 136-173.
class InstantaneousVorticityDeviation2d
{
public:
	static void Compute(const SteadyVectorField2d& flow, int windowSize, SteadyScalarField2d& result) { ComputeSlice(flow, windowSize, result, 0); }
	static void Compute(const UnsteadyVectorField2d& flow, int windowSize, SteadyScalarField2d& result, int slice) { ComputeSlice(flow, windowSize, result, slice); }
	static void Compute(const UnsteadyVectorField2d& flow, int windowSize, UnsteadyScalarField2d& result) { for (int slice = 0; slice < result.GetResolution()[2]; ++slice) ComputeSlice(flow, windowSize, result, slice); }
private:
	template<size_t DimIn, size_t DimOut>
	static void ComputeSlice(const RegularGrid<Vec2f, DimIn>& flow, int windowSize, RegularGrid<float, DimOut>& result, int slice) {
		Vec2i res = result.GetResolution();
		SteadyScalarField2d vorticity(res, BoundingBox2d(Vec2d({ result.GetDomain().Min[0], result.GetDomain().Min[1] }), Vec2d({ result.GetDomain().Max[0], result.GetDomain().Max[1] })));
		Vorticity2d::ComputeSlice(flow, vorticity, slice);
		#pragma omp parallel for
		for (int i = 0; i < res[0] * res[1]; ++i) {
			Vec2i gridCoord = vorticity.GetGridCoord(i);
			int ix = gridCoord[0], iy = gridCoord[1], ix0 = std::max(0, ix - windowSize), ix1 = std::min(ix + windowSize, res[0] - 1), iy0 = std::max(0, iy - windowSize), iy1 = std::min(iy + windowSize, res[1] - 1);
			double avg = 0;
			for (int iyy = iy0; iyy <= iy1; ++iyy)
				for (int ixx = ix0; ixx <= ix1; ++ixx)
					avg += vorticity.GetVertexDataAt(Vec2i({ ixx, iyy }));
			avg /= ((double)ix1 - ix0 + 1) * ((double)iy1 - iy0 + 1);
			result.SetVertexDataAt(Vec3i({ gridCoord[0], gridCoord[1], slice }), vorticity.GetVertexDataAt(gridCoord) - avg);
		}
	}
};

// Integrates a quantity along a trajectory using time parameterization.
class LagrangianAveraging
{
public:
	static void Compute(const SteadyVectorField2d& flow, const SteadyScalarField2d& scalar, double stepSize, double duration, SteadyScalarField2d& result) { ComputeSlice(flow, scalar, stepSize, duration, result, 0); }
	static void Compute(const UnsteadyVectorField2d& flow, const UnsteadyScalarField2d& scalar, double stepSize, double duration, SteadyScalarField2d& result, int slice) { ComputeSlice(flow, scalar, stepSize, duration, result, slice); }
	static void Compute(const UnsteadyVectorField2d& flow, const UnsteadyScalarField2d& scalar, double stepSize, double duration, UnsteadyScalarField2d& result) { for (int slice = 0; slice < result.GetResolution()[2]; ++slice) ComputeSlice(flow, scalar, stepSize, duration, result, slice); }
private:
	template<size_t DimIn, size_t DimOut>
	static void ComputeSlice(const RegularGrid<Vec2f, DimIn>& flow, const RegularGrid<float, DimIn>& scalar, double stepSize, double duration, RegularGrid<float, DimOut>& result, int slice) {
		Vec2i res = result.GetResolution();
		#pragma omp parallel for schedule(guided)
		for (int i = 0; i < res[0] * res[1]; ++i) {
			Vec2i gridCoord = result.GetGridCoord(i);
			Vec2d coord2 = result.GetCoordAt(gridCoord);
			Vec3d coord3 = flow.GetCoordAt(Vec3i({ 0,0,slice }));
			Vec3d coord({ coord2[0], coord2[1], coord3[2] });
			std::vector<typename RegularGrid<Vec2d, DimIn>::TDomainCoord> curve;
			Tracer::Pathline(flow, coord, stepSize, duration, (int)(duration / stepSize), curve);
			double avg = 0;
			for (size_t vid = 0; vid < curve.size(); ++vid)
				avg += scalar.Sample(curve[vid]);
			if (!curve.empty()) avg /= curve.size();
			result.SetVertexDataAt(Vec3i({ gridCoord[0], gridCoord[1], slice }), avg);
		}
	}
};

// Class that fills a 2D scalar field with uniform noise between [0,1].
class Noise
{
public:
	static void Compute(SteadyScalarField2d& result, unsigned int randomSeed = 1337) { 
		std::default_random_engine rng(randomSeed);
		std::uniform_real_distribution<float> rnd(0, 1);
		for (int iy = 0; iy < result.GetResolution()[1]; ++iy)
			for (int ix = 0; ix < result.GetResolution()[0]; ++ix)
				result.SetVertexDataAt(Vec2i({ ix, iy }), rnd(rng));
	}
};

// Computes a line integral convolution (LIC).
// Cabral, B., & Leedom, L. C. (1993, September). Imaging vector fields using line integral convolution. In Proceedings of the 20th annual conference on Computer graphics and interactive techniques (pp. 263-270).
class LineIntegralConvolution
{
public:
	static void Compute(const SteadyVectorField2d& flow, const SteadyScalarField2d& scalar, double stepSize, double duration, SteadyScalarField2d& result) { ComputeSlice(flow, scalar, stepSize, duration, result, 0); }
	static void Compute(const UnsteadyVectorField2d& flow, const SteadyScalarField2d& scalar, double stepSize, double duration, SteadyScalarField2d& result, int slice) { ComputeSlice(flow, scalar, stepSize, duration, result, slice); }
	static void Compute(const UnsteadyVectorField2d& flow, const SteadyScalarField2d& scalar, double stepSize, double duration, UnsteadyScalarField2d& result) { for (int slice = 0; slice < result.GetResolution()[2]; ++slice) ComputeSlice(flow, scalar, stepSize, duration, result, slice); }
private:
	template<size_t DimIn, size_t DimOut>
	static void ComputeSlice(const RegularGrid<Vec2f, DimIn>& flow, const SteadyScalarField2d& noise, double stepSize, double duration, RegularGrid<float, DimOut>& result, int slice) {
		Vec2i res = result.GetResolution();
		#pragma omp parallel for schedule(guided)
		for (int i = 0; i < res[0] * res[1]; ++i) {
			Vec2i gridCoord = result.GetGridCoord(i);
			Vec2d coord2 = result.GetCoordAt(gridCoord);
			Vec3d coord3 = flow.GetCoordAt(Vec3i({ 0,0,slice }));
			Vec3d coord({ coord2[0], coord2[1], coord3[2] });
			const double sigma = 0.3;
			double count = 0, sum = 0;
			std::vector<typename RegularGrid<Vec2d, DimIn>::TDomainCoord> curve;
			Tracer::Streamline(flow, coord, stepSize, duration, (int)(duration / stepSize), curve);
			for (size_t vid = 0; vid < curve.size() - 1; ++vid) {
				double stepLength = std::sqrt((curve[vid] - curve[vid + 1]).lengthSquared());
				double x = vid / (curve.size() - 1.0);
				double kernel = 1.0 / (2.5006 * sigma) * exp(-(x * x) / (2 * sigma * sigma));
				sum += noise.Sample(curve[vid + 1]) * kernel * stepLength;
				count += kernel * stepLength;
			}
			Tracer::Streamline(flow, coord, -stepSize, duration, (int)(duration / stepSize), curve);
			for (size_t vid = 0; vid < curve.size() - 1; ++vid) {
				double stepLength = std::sqrt((curve[vid] - curve[vid + 1]).lengthSquared());
				double x = vid / (curve.size() - 1.0);
				double kernel = 1.0 / (2.5006 * sigma) * exp(-(x * x) / (2 * sigma * sigma));
				sum += noise.Sample(curve[vid + 1]) * kernel * stepLength;
				count += kernel * stepLength;
			}
			// scale the gray value range
			double value = count > 0 ? std::min(std::max(0.0, sum / count), 1.0) : 0.5;
			value = std::min(std::max(0.0, (value - 0.5) * (std::min(std::max(0.0, count / 1.57), 1.0) * 3) + 0.5), 1.0);
			result.SetVertexDataAt(Vec3i({ gridCoord[0], gridCoord[1], slice }), value);
		}
	}
};

// Advects a given scalar texture over time by fetching it backwards.
class TextureAdvection
{
public:
	static void Compute(const SteadyVectorField2d& flow, const SteadyScalarField2d& scalar, double stepSize, double duration, SteadyScalarField2d& result) { ComputeSlice(flow, scalar, stepSize, duration, result, 0); }
	static void Compute(const UnsteadyVectorField2d& flow, const SteadyScalarField2d& scalar, double stepSize, double duration, SteadyScalarField2d& result, int slice) { ComputeSlice(flow, scalar, stepSize, duration, result, slice); }
	static void Compute(const UnsteadyVectorField2d& flow, const SteadyScalarField2d& scalar, double stepSize, double duration, UnsteadyScalarField2d& result) { for (int slice = 0; slice < result.GetResolution()[2]; ++slice) ComputeSlice(flow, scalar, stepSize, duration, result, slice); }
private:
	template<size_t DimIn, size_t DimOut>
	static void ComputeSlice(const RegularGrid<Vec2f, DimIn>& flow, const SteadyScalarField2d& tex, double stepSize, double duration, RegularGrid<float, DimOut>& result, int slice) {
		Vec2i res = result.GetResolution();
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < res[0] * res[1]; ++i) {
			Vec3i gridCoord({ result.GetGridCoord(i)[0], result.GetGridCoord(i)[1], slice });
			Vec2d coord2 = result.GetCoordAt(gridCoord);
			Vec3d coord3 = flow.GetCoordAt(Vec3i({ 0,0,slice }));
			Vec3d coord({ coord2[0], coord2[1], coord3[2] });
			bool indomain = true;
			Vec3d phi = Tracer::FlowMap(flow, coord, -stepSize, duration, indomain);
			if (indomain)
				result.SetVertexDataAt(Vec3i({ gridCoord[0], gridCoord[1], slice }), tex.Sample(phi));
			else result.SetVertexDataAt(Vec3i({ gridCoord[0], gridCoord[1], slice }), 0);
		}
	}
};