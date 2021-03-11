#include "mex.h"
#include "vislcs.hpp"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	Vec3i resolution;   // resolution of the space-time grid
    Vec3d minDomain;    // min corner of the physical space-time domain
    Vec3d maxDomain;    // max corner of the physical space-time domain
    Vec2f* data;        // velocity components at the grid points
    double stepSize;    // integration step size
    double duration;    // integration duration
    Vec3i out_res;      // resolution of the output grid
    int slice;          // time slice to compute for
    float* result;      // pointer to the output array

    if (nrhs < 7) {
        mexErrMsgIdAndTxt("nrhs","At least seven inputs required.");
        return;
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("nlhs","One output required.");
        return;
    }

	/* Input arguments */
	resolution = *(Vec3i*)mxGetPr(prhs[0]);
    minDomain = *(Vec3d*)mxGetPr(prhs[1]);
    maxDomain = *(Vec3d*)mxGetPr(prhs[2]);
    data = (Vec2f*)mxGetPr(prhs[3]);
    stepSize = mxGetScalar(prhs[4]);
    duration = mxGetScalar(prhs[5]);
    out_res = *(Vec3i*)mxGetPr(prhs[6]);
    slice = nrhs>=8 ? (int)mxGetScalar(prhs[7]) : 0;
    
    /* Output arguments */
    plhs[0] = mxCreateNumericMatrix(out_res[0], out_res[1] * (nrhs>=8 ? 1 : out_res[2]), mxSINGLE_CLASS, mxREAL);
    result = (float*) mxGetData(plhs[0]);

    /* Perform the computation */
    UnsteadyVectorField2d inField(resolution, BoundingBox3d(minDomain, maxDomain), data);
    if (nrhs>=8) {
        SteadyScalarField2d outField(Vec2i({ out_res[0], out_res[1] }), BoundingBox2d(Vec2d({ minDomain[0], minDomain[1] }), Vec2d({ maxDomain[0], maxDomain[1] })), result);
        FiniteTimeLyapunovExponent2d::Compute(inField, stepSize, duration, outField, slice);
    }
    else {
        UnsteadyScalarField2d outField(out_res, BoundingBox3d(minDomain, maxDomain), result);
        FiniteTimeLyapunovExponent2d::Compute(inField, stepSize, duration, outField);	
    }
}
