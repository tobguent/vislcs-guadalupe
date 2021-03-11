#include "mex.h"
#include "vislcs.hpp"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	Vec3i resolution;   // resolution of the space-time grid
    Vec3d minDomain;    // min corner of the physical space-time domain
    Vec3d maxDomain;    // max corner of the physical space-time domain
    Vec2f* data;        // velocity components at the grid points
    int windowSize;     // size of the neighborhood region
    int slice;          // time slice to compute for
    float* result;      // pointer to the output array

    if (nrhs < 5) {
        mexErrMsgIdAndTxt("nrhs","At least five inputs required.");
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
    windowSize = (int)mxGetScalar(prhs[4]);
    slice = nrhs>=6 ? (int)mxGetScalar(prhs[5]) : 0;
    
    /* Output arguments */
    plhs[0] = mxCreateNumericMatrix(resolution[0], resolution[1] * (nrhs>=6 ? 1 : resolution[2]), mxSINGLE_CLASS, mxREAL);
    result = (float*) mxGetData(plhs[0]);

    /* Perform the computation */
    UnsteadyVectorField2d inField(resolution, BoundingBox3d(minDomain, maxDomain), data);
    if (nrhs>=6) {
        SteadyScalarField2d outField(Vec2i({ resolution[0], resolution[1] }), BoundingBox2d(Vec2d({ minDomain[0], minDomain[1] }), Vec2d({ maxDomain[0], maxDomain[1] })), result);
        InstantaneousVorticityDeviation2d::Compute(inField, windowSize, outField, slice);
    }
    else {
        UnsteadyScalarField2d outField(resolution, BoundingBox3d(minDomain, maxDomain), result);
        InstantaneousVorticityDeviation2d::Compute(inField, windowSize, outField);	
    }
}
