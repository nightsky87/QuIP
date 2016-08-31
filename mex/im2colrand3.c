#include "mex.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize h, w;
    mwSize bsizeX, bsizeY, psize, numPatches;
    mwIndex n, x, y, i;

    double *bsize, *Y, *Cb, *Cr, *PY, *PCb, *PCr;
    
    // Extract image dimensions
    h = mxGetM(prhs[0]);
    w = mxGetN(prhs[0]) / 3;
    
    // Get the pointer to the data
    Y = mxGetPr(prhs[0]);
    Cb = Y + h * w;
    Cr = Cb + h * w;
    bsize = mxGetPr(prhs[1]);
    numPatches = mxGetScalar(prhs[2]);
    
    // Copy the block size parameters
    bsizeY = bsize[0];
    bsizeX = bsize[1];
    psize = bsizeX * bsizeY;
    
    // Create the output matrix
    plhs[0] = mxCreateDoubleMatrix(3 * psize, numPatches, mxREAL);
    PY = mxGetPr(plhs[0]);
    PCb = PY + psize;
    PCr = PCb + psize;
    
    // Seed the random number generator
    srand(0);
    
    for (n = 0; n < numPatches; n++)
    {
        x = rand() % (w - bsizeX + 1);
        y = rand() % (h - bsizeY + 1);
        
        for (i = 0; i < bsizeX; i++)
        {
            memcpy(PY, &Y[h*(x+i)+y], bsizeY * sizeof(double));
            memcpy(PCb, &Cb[h*(x+i)+y], bsizeY * sizeof(double));
            memcpy(PCr, &Cr[h*(x+i)+y], bsizeY * sizeof(double));
            
            PY += bsizeY;
            PCb += bsizeY;
            PCr += bsizeY;
        }
        
        PY += 2 * psize;
        PCb += 2 * psize;
        PCr += 2 * psize;
    }
}