#include "mex.h"
#include <string.h>

typedef unsigned char  uint8;
typedef unsigned short uint16;
typedef unsigned long  uint32;

#define SCALE   (128)
#define MAXRATE (1.5 * numDims * SCALE)

#define D(n,q)  (distTab[(q << bs) + n])
#define R(n,q)  (rateTab[(q << bs) + n])
#define LD(n,s) (leastDist[(s << bs) + n])
#define QC(n,s) (qChoice[(s << bs) + n])

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize numDims, numCols;
    mwIndex i, j, q;
    mwIndex bs;
	
	uint16 s, qVal;
    uint16 *P, *qChoice, *rateTab;
	uint32 hist[32768], histQ[32768];
	
	double p;
	double *leastDist, *distTab;
	
    // Extract signal dimensions
    numDims = mxGetM(prhs[0]);
    numCols = mxGetN(prhs[0]);
	
	// Ensure that the patches are either 4x4 or 8x8
	if (numDims == 16)
	{
		bs = 4;
	}
	else if (numDims == 64)
	{
		bs = 6;
	}
	else
	{
		mexErrMsgTxt("Invalid block size used. Only 4x4 and 8x8 blocks are allowed.");
	}
	
	// Ensure that the DCT coefficients are of uint16 type
	if (mxGetClassID(prhs[0]) != mxUINT16_CLASS)
	{
		mexErrMsgTxt("DCT coefficients have to be in uint16 format (discarding sign information).");
	}
	
    // Get the pointer to the patch data
    P = (uint16 *)mxGetData(prhs[0]);
	
	// Construct the output tables
	plhs[0] = mxCreateDoubleMatrix(numDims, MAXRATE, mxREAL); // leastDist
	plhs[1] = mxCreateNumericMatrix(numDims, MAXRATE, mxUINT16_CLASS, mxREAL); // qChoice
	plhs[2] = mxCreateNumericMatrix(numDims, 4096, mxUINT16_CLASS, mxREAL); // rateTab
	plhs[3] = mxCreateDoubleMatrix(numDims, 4096, mxREAL); // distTab
	
	// Copy pointers
	leastDist = mxGetPr(plhs[0]);
	qChoice = (uint16 *)mxGetData(plhs[1]);
	rateTab = (uint16 *)mxGetData(plhs[2]);
	distTab = mxGetPr(plhs[3]);
	
	// Process each dimension
	for (i = 0; i < numDims; i++)
	{
		// Initialize the symbol histograms
		for (j = 0; j < 32768; j++)
		{
			hist[j] = 0;
		}
		
		// Calculate the unquantized symbol histogram
		for (j = 0; j < numCols; j++)
		{
			hist[P[(j << bs) + i]]++;
		}
		
		// Characterize the rate and distortion response at different quantization levels
		for (q = 1; q <= 4096; q++)
		{
			// Initialize the quantized symbol histograms
			for (j = 0; j < 32768; j++)
			{
				histQ[j] = 0;
			}
			
			// Process each quantized value
			for (j = 0; j < 32768; j++)
			{
				// Calculate the quantized value
				qVal = floor((j + (q >> 1)) / q);
				
				// Update the quantized histogram
				histQ[qVal] += hist[j];
				
				// Update the error table
				D(i, q-1) += ((q * qVal - j) * (q * qVal - j) * hist[j]);
			}
			
			// Calculate the rate
			for (j = 0; j < 32768; j++)
			{
				if (histQ[j] > 0)
				{
					// Calculate the probability
					p = (double)histQ[j] / (double)numCols;
					R(i, q-1) += (uint16)(-SCALE * p * log2(p));
				}
			}
		}
	}
	
	// Initialize the least distortion table
	for (i = 0; i < MAXRATE * numDims; i++)
	{
		leastDist[i] = 1e99;
 	}
	
	// Generate the first row of the least distortion table
	for (q = 0; q < 4096; q++)
	{
		if (D(0, q) < LD(0, R(0, q)))
		{
			LD(0, R(0, q)) = D(0, q);
			QC(0, R(0, q)) = q;
		}
	}
	
	// Generate the least distortion table
	for (i = 1; i < numDims; i++)
	{
		for (q = 0; q < 4096; q++)
		{
			for (j = 0; j < MAXRATE; j++)
			{
				s = j + R(i, q);
				
				if (s > MAXRATE)
					continue;
				
				if ((D(i, q) + LD(i-1, j)) < LD(i, s))
				{
					LD(i, s) = D(i, q) + LD(i-1, j);
					QC(i, s) = q;
				}
			}
		}
	}
	
	// Adjust for 1-indexing in MATLAB
	for (i = 0; i < MAXRATE * numDims; i++)
	{
		qChoice[i]++;
 	}
}