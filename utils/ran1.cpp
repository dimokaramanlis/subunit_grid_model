/*
 * Copyright (c) 2008, Christian Mendl
 * All rights reserved.
 *
 */

// Pseudo random number generator 'ran1' used in the stimulation program.

#include <mex.h>
#include <memory.h>


const long NTAB = 32;

// seed also contains 'static' variables used in 'ran1'
struct Seed
{
	long idum;
	long iy;
	long iv[NTAB];

	Seed()
	: idum(-1), iy(0) {
	}
};

// function declarations
double ran1(Seed& seed);	// random generator
mxArray *Seed2Mat(const Seed& seed);
void Mat2Seed(const mxArray *mat, Seed& seed);

/* x = ran1(seed)
 * x = ran1(seed,n)
 * [x,seed] = ran1(seed)
 * [x,seed] = ran1(seed,n)
 * seed = ran1(seed,n,0)
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs == 0 || nrhs > 3) mexErrMsgTxt("One, two or three input arguments required.");
	if (nlhs > 2) mexErrMsgTxt("Wrong number of output arguments.");
	if (nlhs > 1 && nrhs == 3) mexErrMsgTxt("One output argument expected given three input arguments.");

	Seed seed;
	if (mxIsDouble(prhs[0]))
		seed.idum = (long)mxGetScalar(prhs[0]);
	else if (mxIsInt32(prhs[0]))
		memcpy(&seed.idum, mxGetData(prhs[0]), sizeof(long));
	else if (mxIsStruct(prhs[0]))
		Mat2Seed(prhs[0], seed);
	else
		mexErrMsgTxt("Invalid class of input argument 'seed'.");

	// number of generated random numbers
	mwSize n;
	if (nrhs < 2) n = 1;
	else {
		n = (mwSize)mxGetScalar(prhs[1]);
		if (n < 1) mexErrMsgTxt("'n' must be a positive integer.");
	}

	if (nrhs <= 2)
	{
		// generate random numbers
		double *x = new double[n];
		for (mwSize i = 0; i < n; i++)
			x[i] = ran1(seed);

		// copy output: random numbers
		plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
		memcpy(mxGetPr(plhs[0]), x, n*sizeof(double));
		delete [] x;
		// seed
		if (nlhs == 2) {
			plhs[1] = Seed2Mat(seed);
		}
	}
	else
	{
		// omit random numbers, just advance the seed value
		for (mwSize i = 0; i < n; i++)
			ran1(seed);

		// copy output: seed
		plhs[0] = Seed2Mat(seed);
	}
}


// convert Matlab 'mxArray' to the 'seed' structure and vice versa

inline void Mat2Seed(const mxArray *mat, Seed& seed)
{
	// copy 'idum'
	mxArray *idum = mxGetField(mat, 0, "idum");
	if (!idum) mexErrMsgTxt("Structure field 'idum' not found.");
	if (!mxIsInt32(idum)) mexErrMsgTxt("Structure field 'idum' should be of class 'int32'.");
	memcpy(&seed.idum, mxGetData(idum), sizeof(seed.idum));
	// copy 'iy'
	mxArray *iy = mxGetField(mat, 0, "iy");
	if (!iy) mexErrMsgTxt("Structure field 'iy' not found.");
	if (!mxIsInt32(iy)) mexErrMsgTxt("Structure field 'iy' should be of class 'int32'.");
	memcpy(&seed.iy, mxGetData(iy), sizeof(seed.iy));
	// copy 'iv'
	mxArray *iv = mxGetField(mat, 0, "iv");
	if (!iv) mexErrMsgTxt("Structure field 'iv' not found.");
	if (!mxIsInt32(iv)) mexErrMsgTxt("Structure field 'iv' should be of class 'int32'.");
	if (!(mxGetM(iv) == NTAB && mxGetN(iv) == 1) && !(mxGetM(iv) == 1 && mxGetN(iv) == NTAB))
		mexErrMsgTxt("Wrong dimension of structure field 'iv'.");
	memcpy(&seed.iv, mxGetData(iv), sizeof(seed.iv));
}

inline mxArray *Seed2Mat(const Seed& seed)
{
	// create seed structure fields
	mxArray *idum = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
	memcpy(mxGetData(idum), &seed.idum, sizeof(seed.idum));
	mxArray *iy = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
	memcpy(mxGetData(iy), &seed.iy, sizeof(seed.iy));
	mxArray *iv = mxCreateNumericMatrix(1, NTAB, mxINT32_CLASS, mxREAL);
	memcpy(mxGetData(iv), &seed.iv, sizeof(seed.iv));

	// create return value structure
	const char *fieldnames[] = { "idum", "iy", "iv" };
	mxArray *ret = mxCreateStructMatrix(1, 1, 3, fieldnames);
	mxSetFieldByNumber(ret, 0, 0, idum);
	mxSetFieldByNumber(ret, 0, 1, iy);
	mxSetFieldByNumber(ret, 0, 2, iv);

	return ret;
}


const long IA = 16807;
const long IM = 2147483647;
const double AM = 1.0/IM;
const long IQ = 127773;
const long IR = 2836;
const long NDIV = 1+(IM-1)/NTAB;
const double EPS = 1.2e-7;
const double RNMX = 1.0-EPS;

double ran1(Seed& seed)
{
	int j;
	long k;
	// static long iy=0;
	// static long iv[NTAB];
	double temp;

	if (seed.idum <= 0 || !seed.iy) {
		if (-seed.idum < 1) seed.idum=1;
		else seed.idum = -seed.idum;
		for (j=NTAB+7;j>=0;j--) {
			k=seed.idum/IQ;
			seed.idum=IA*(seed.idum-k*IQ)-IR*k;
			if (seed.idum < 0) seed.idum += IM;
			if (j < NTAB) seed.iv[j] = seed.idum;
		}
		seed.iy=seed.iv[0];
	}
	k=(seed.idum)/IQ;
	seed.idum=IA*(seed.idum-k*IQ)-IR*k;
	if (seed.idum < 0) seed.idum += IM;
	j=seed.iy/NDIV;
	seed.iy=seed.iv[j];
	seed.iv[j] = seed.idum;
	if ((temp=AM*seed.iy) > RNMX) return RNMX;
	else return temp;
}
