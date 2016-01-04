#include "eigen.h"

#define   SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
#define eps 0.0001

void reduce(float **a, int n, float *d, float *e)
{
	register int    l, k, j, i;
	float          scale, hh, h, g, f;

	for (i = n - 1; i > 0; i--) {
		l = i - 1;
		h = scale = 0.0;
		if (l > 0) {
			for (k = 0; k <= l; k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)               // skip transformation
				e[i] = a[i][l];
			else {
				for (k = 0; k <= l; k++) {
					a[i][k] /= scale;          // used scaled a's for transformation
					h += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i] = scale*g;
				h -= f * g;
				a[i][l] = f - g;
				f = 0.0;

				for (j = 0; j <= l; j++) {
					a[j][i] = a[i][j] / h;       // can be omitted if eigenvector not wanted
					g = 0.0;
					for (k = 0; k <= j; k++) {
						g += a[j][k] * a[i][k];
					}
					for (k = j + 1; k <= l; k++)
						g += a[k][j] * a[i][k];
					e[j] = g / h;
					f += e[j] * a[i][j];
				}
				hh = f / (h + h);
				for (j = 0; j <= l; j++) {
					f = a[i][j];
					e[j] = g = e[j] - hh*f;
					for (k = 0; k <= j; k++)
						a[j][k] -= (f*e[k] + g*a[i][k]);
				}
			}  // end k-loop
		}  // end if-loop for l > 1
		else {
			e[i] = a[i][l];
		}
		d[i] = h;
	}  // end i-loop
	d[0] = 0.0;
	e[0] = 0.0;

	/* Contents of this loop can be omitted if eigenvectors not
	* 	 ** wanted except for statement d[i]=a[i][i];
	* 	          */

	for (i = 0; i < n; i++) {
		l = i - 1;
		if (d[i]) {
			for (j = 0; j <= l; j++) {
				g = 0.0;
				for (k = 0; k <= l; k++) {
					g += a[i][k] * a[k][j];
				}
				for (k = 0; k <= l; k++) {
					a[k][j] -= g * a[k][i];
				}
			}
		}
		d[i] = a[i][i];
		a[i][i] = 1.0;
		for (j = 0; j <= l; j++) {
			a[j][i] = a[i][j] = 0.0;
		}
	}
} 

float pythag(float a, float b)
{
	float absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return
		absa*sqrt(1.0 + (absb / absa)*(absb / absa));
	else return
		(absb == 0.0 ? 0.0 : absb*sqrt(1.0 + (absa / absb)*(absa / absb)));
}

void findEigen(float *d, float *e, int n, float **z)
{

	register int   m, l, iter, i, k;
	float         s, r, p, g, f, dd, c, b;


	for (i = 1; i < n; i++) e[i - 1] = e[i];
	e[n] = 0.0;
	for (l = 0; l < n; l++) {
		iter = 0;
		do {
			for (m = l; m < n - 1; m++) {
				dd = fabs(d[m]) + fabs(d[m + 1]);
				if ((double)(fabs(e[m]) + dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) {
					std::cerr << "\n\nToo many iterations in findEigen.\n";
					exit(1);
				}
				g = (d[l + 1] - d[l]) / (2.0 * e[l]);
				r = pythag(g, 1.0);
				g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
				s = c = 1.0;
				p = 0.0;
				for (i = m - 1; i >= l; i--) {
					f = s * e[i];
					b = c*e[i];
					e[i + 1] = (r = pythag(f, g));
					if (r == 0.0) {
						d[i + 1] -= p;
						e[m] = 0.0;
						break;
					}
					s = f / r;
					c = g / r;
					g = d[i + 1] - p;
					r = (d[i] - g) * s + 2.0 * c * b;
					d[i + 1] = g + (p = s * r);
					g = c * r - b;
					for (k = 0; k < n; k++) {
						f = z[k][i + 1];
						z[k][i + 1] = s * z[k][i] + c * f;
						z[k][i] = c * z[k][i] - s * f;
					} /* end k-loop */
				} /* end i-loop */
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			} /* end if-loop for m != 1 */
		} while (m != l);
	}
	/* end l-loop */
} 

void sortEigen(float *d, float **v, int n) {
	int  k, j, i;
	float p;
	for (i = 0; i < n - 1; i++) {
		p = d[k = i];
		for (j = i + 1; j < n; j++)  if (d[j] >= p) p = d[k = j];
		if (k != i) {
			d[k] = d[i]; d[i] = p;
			for (j = 0; j < n; j++) {
				p = v[j][i]; v[j][i] = v[j][k]; v[j][k] = p;
			}
		}
	}
}

void SVD(Matrix<float> &M, Matrix<float> &U, Matrix<float> &V)
{
	Matrix <float> MtM = M.transpose() * M;
	Matrix <float> MMt = M * M.transpose();

	//Starting transformations for M*Mt and Mt*M
	float* d = new float[MMt.getRows()];
	float* e = new float[MMt.getRows()];
	float* dw = new float[MtM.getRows()];
	float* ew = new float[MtM.getRows()];

	//After transformations M*Mt is no longer valid
	//instead after QL iterations it contains eigenvectors of M*Mt
	reduce(MMt, MMt.getRows(), d, e);
	findEigen(d, e, MMt.getRows(), MMt);
	sortEigen(d, MMt, MMt.getRows());

	U = MMt;

	//Do the same for Mt*M matrix
	reduce(MtM, MtM.getRows(), dw, ew);
	findEigen(dw, ew, MtM.getRows(), MtM);
	sortEigen(dw, MtM, MtM.getRows());

	V = MtM;

	//correcting signes of vectors in V
	for (int i = 0; i < V.getRows(); i++)
		for (int j = 1; j < V.getColumns(); j++)
			V[i][j] = -V[i][j];

	//Matrix D contains singular values in descending order
	//on it's diagonal and zeros on other places
	Matrix<float> D(U.getRows(), V.getRows());
	for (int i = 0; i < D.getRows(); i++)
		for (int j = 0; j < D.getColumns(); j++)
		{
			if (i == j && d[i] > eps) {
				D[i][j] = sqrt(d[i]);
			}
			else D[i][j] = 0;
		}
	M = D;
}