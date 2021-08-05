#include "mex.h"
#include "vislcs.hpp"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	Vec3i resolution;   // resolution of the space-time grid
    Vec3d minDomain;    // min corner of the physical space-time domain
    Vec3d maxDomain;    // max corner of the physical space-time domain
    Vec2f* data;        // velocity components at the grid points
    double stepSize;    // integration step size
    double t0;          // time slice to compute streamlines for
    double duration;    // integration duration
    double dTest;       // lines will not get closer too each other than
    double dSep;        // new lines are seeded with an offset of
    double* result0;    // pointer to the output vertices
    int* result1;       // pointer to the output start offsets
    int* result2;       // pointer to the output number of vertices

    if (nrhs < 9) {
        mexErrMsgIdAndTxt("nrhs","At least nine inputs required.");
        return;
    }
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("nlhs","Three output required.");
        return;
    }

	/* Input arguments */
	resolution = *(Vec3i*)mxGetPr(prhs[0]);
    minDomain = *(Vec3d*)mxGetPr(prhs[1]);
    maxDomain = *(Vec3d*)mxGetPr(prhs[2]);
    data = (Vec2f*)mxGetPr(prhs[3]);
    stepSize = mxGetScalar(prhs[4]);
    t0 = mxGetScalar(prhs[5]);
    duration = mxGetScalar(prhs[6]);
    dTest = mxGetScalar(prhs[7]);
    dSep = mxGetScalar(prhs[8]);
    
    
    /* Perform the computation */
    UnsteadyVectorField2d inField(resolution, BoundingBox3d(minDomain, maxDomain), data);
    std::vector<Vec2d> vertices;
	std::vector<int> offset;
	std::vector<int> length;
	EvenlySpacedStreamlines::Compute(inField, stepSize, t0, duration, dTest, dSep, vertices, offset, length);
    
    /* Output arguments */
    plhs[0] = mxCreateDoubleMatrix(vertices.size(), 2, mxREAL);
    result0 = (double*) mxGetData(plhs[0]);
    memcpy(result0, &vertices.data()[0], sizeof(Vec2d) * vertices.size());
    
    plhs[1] = mxCreateNumericMatrix(offset.size(), 1, mxINT32_CLASS, mxREAL);
    result1 = (int*) mxGetData(plhs[1]);
    memcpy(result1, &offset.data()[0], sizeof(int) * offset.size());
    
    plhs[2] = mxCreateNumericMatrix(length.size(), 1, mxINT32_CLASS, mxREAL);
    result2 = (int*) mxGetData(plhs[2]);
    memcpy(result2, &length.data()[0], sizeof(int) * length.size());
}
