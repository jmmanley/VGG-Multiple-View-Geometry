#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

//void ludcmp(double **a, int n, int *indx, float *d)
void ludcmp(double **a, int n, double *d, double *vv)
{
	int i,imax,j,k;
	double big,dum,sum,temp;

	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) mexErrMsgTxt("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		//indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=1.0e-20;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}

}


/* D = detn(X) vectorized determinant; D(i1,...,in) = det(X(:,:,i1,...,in)). */
void mexFunction( int nargout,
                  mxArray **argout,
                  int nargin,
                  const mxArray **argin
                )
{
  int ndims, n, i, j, k, N;
  const int *dim;
  double *X, *x, *D, **a, *vv, *d;

  if ( nargin != 1 ) mexErrMsgTxt("One input parameter required.");

  if ( mxGetClassID(argin[0]) != mxDOUBLE_CLASS ) mexErrMsgTxt("X must be double.");
  ndims = mxGetNumberOfDimensions(argin[0]);
  if ( ndims < 2 ) mexErrMsgTxt("X must have 2+ dimensions.");
  dim = mxGetDimensions(argin[0]);
  if ( dim[0]!=dim[1] ) mexErrMsgTxt("X must have first two dimensions equal.");
  X = (double*)mxGetPr(argin[0]);

  if ( ndims>2 )
    argout[0] = mxCreateNumericArray(ndims-2,dim+2,mxDOUBLE_CLASS,mxREAL);
  else
    argout[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  if ( argout[0]==0 ) mexErrMsgTxt("Out of memory.");
  D = (double*)mxGetPr(argout[0]);

  // Allocate aux. array a
  a = (double**)malloc((dim[0]+1)*sizeof(double*));
  for ( i = 0; i <= dim[0]; i++ ) a[i] = (double*)malloc((dim[1]+1)*sizeof(double));

  vv=(double*)malloc((dim[0]+1)*sizeof(double));

  N = 1; for ( i = 2; i < ndims; i++ ) N *= dim[i]; // N = number of determinants
  for ( n = 0, x = X, d = D; n < N; n++, x+=dim[0]*dim[1], d++ ) {
    if ( dim[0]==2 )
      *d = x[0]*x[3]-x[2]*x[1];
    else if ( dim[0]==3 )
      *d = x[0]*x[4]*x[8]+x[2]*x[3]*x[7]+x[1]*x[5]*x[6]-x[2]*x[4]*x[6]-x[1]*x[3]*x[8]-x[0]*x[5]*x[7];
    else {
      k = 0; for ( i = 1; i <= dim[0]; i++) for ( j = 1; j <= dim[1]; j++ ) a[j][i] = x[k++];
      ludcmp(a,dim[0],d,vv);
      for ( i=1;i<=dim[0];i++) *d *= a[i][i];
    }
  }  

  free(vv);
  for ( i = 0; i < dim[0]; i++ ) free(a[i]);
  free(a);
}
