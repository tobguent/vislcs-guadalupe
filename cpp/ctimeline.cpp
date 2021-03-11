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
    Vec3d seedA;        // beginning of the seed line
    Vec3d seedB;        // end of the seed line
    int numSteps;       // number of time steps
    double* result;     // pointer to the output array

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
    stepSize = mxGetScalar(prhs[4]);
    duration = mxGetScalar(prhs[5]);
    seedA = *(Vec3d*)mxGetPr(prhs[6]);
    seedB = *(Vec3d*)mxGetPr(prhs[7]);
    numSteps = (int)mxGetScalar(prhs[8]);
    
    /* Perform the computation */
    UnsteadyVectorField2d inField(resolution, BoundingBox3d(minDomain, maxDomain), data);
    std::vector<Vec3d> curve;
    Tracer::Timeline(inField, seedA, seedB, stepSize, duration, numSteps, curve);
    
    /* Output arguments */
    plhs[0] = mxCreateDoubleMatrix(curve.size(), 3, mxREAL);
    result = (double*) mxGetData(plhs[0]);
    memcpy(result, &curve.data()[0], sizeof(Vec3d) * curve.size());
}
