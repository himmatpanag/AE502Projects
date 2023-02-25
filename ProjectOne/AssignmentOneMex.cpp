#include "mex.h"
#include "UniversalPropagator.h"
#include <vector>
//To run this function, type the following into terminal:
// ' /Applications/MATLAB_R2021b.app/bin/mex AssignmentOneMex.cpp UniversalPropagator.cpp'

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check number of inputs and outputs
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("my_mex_file:nrhs", "Two inputs expected.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("my_mex_file:nlhs", "One output required.");
    }
    
    int targetNum = mxGetScalar(prhs[0]);
    bool rendezvous = mxIsLogicalScalarTrue(prhs[1]);

    std::vector<std::vector<double>> result = AssignmentOne(targetNum, rendezvous);

    // Write the result to the output mxArray
    mwSize m = result.size(); //Number of rows
    mwSize n = result[0].size(); // Number of columns
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);

    // copy the result to the output mxArray
    double* output = mxGetPr(plhs[0]);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            output[j * m + i] = result[i][j];
        }
    }

}