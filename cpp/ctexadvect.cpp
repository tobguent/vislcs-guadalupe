#include "mex.h"
#include "vislcs.hpp"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	Vec3i resolution;   // resolution of the space-time grid
    Vec3d minDomain;    // min corner of the physical space-time domain
    Vec3d maxDomain;    // max corner of the physical space-time domain
    Vec2f* data;        // velocity components at the grid points
    float* tex_data;    // scalar texture at the grid points
    Vec2i tex_res;      // noise resolution of the space grid
    double stepSize;    // integration step size
    double duration;    // integration duration
    Vec3i out_res;      // resolution of the output grid
    int slice;          // time slice to compute for
    float* result;      // pointer to the output array

    if (nrhs < 9) {
        mexErrMsgIdAndTxt("nrhs","At least nine inputs required.");
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
    tex_data = (float*)mxGetPr(prhs[4]);
    tex_res = *(Vec2i*)mxGetPr(prhs[5]);
    stepSize = mxGetScalar(prhs[6]);
    duration = mxGetScalar(prhs[7]);
    out_res = *(Vec3i*)mxGetPr(prhs[8]);
    slice = nrhs>=10 ? (int)mxGetScalar(prhs[9]) : 0;
    
    /* Output arguments */
    plhs[0] = mxCreateNumericMatrix(out_res[0], out_res[1] * (nrhs>=10 ? 1 : out_res[2]), mxSINGLE_CLASS, mxREAL);
    result = (float*) mxGetData(plhs[0]);

    /* Perform the computation */
    UnsteadyVectorField2d inField(resolution, BoundingBox3d(minDomain, maxDomain), data);
    SteadyScalarField2d tex(tex_res, BoundingBox2d(Vec2d({ minDomain[0],minDomain[1] }), Vec2d({ maxDomain[0],maxDomain[1] })), tex_data);
    if (nrhs>=10) {
        SteadyScalarField2d outField(Vec2i({ out_res[0], out_res[1] }), BoundingBox2d(Vec2d({ minDomain[0], minDomain[1] }), Vec2d({ maxDomain[0], maxDomain[1] })), result);
        TextureAdvection::Compute(inField, tex, stepSize, duration, outField, slice);
    }
    else {
        UnsteadyScalarField2d outField(out_res, BoundingBox3d(minDomain, maxDomain), result);
        TextureAdvection::Compute(inField, tex, stepSize, duration, outField);	
    }
}
