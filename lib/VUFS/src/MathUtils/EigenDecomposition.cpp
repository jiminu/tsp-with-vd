#include "nrutil.h"

#include "EigenDecomposition.h"


void elmhes(rg_REAL **a, rg_INT n)
//Reduction to Hessenberg form by the elimination method. The real, nonsymmetric matrix
//a[1..n][1..n] is replaced by an upper Hessenberg matrix with identical eigenvalues. Rec-
//ommended, but not required, is that this routine be preceded by balanc. On output, the
//Hessenberg matrix is in elements a[i][j] with i . j+1. Elements with i > j+1 are to be
//thought of as zero, but are returned with random values.
{
	rg_INT m, j, i;
	rg_REAL y, x;
	for(m = 1;m < n - 1;m++) //m is called r + 1 in the text.
	{
		x = 0.0;
		i = m;
		for(j = m;j <= n - 1;j++) // Find the pivot.
		{
			if(fabs(a[ j ][m - 1]) > fabs(x))
			{
				x = a[ j ][m - 1];
				i = j;
			}
		}

		if(i != m) // Interchange rows and columns.
		{
			for(j = m - 1;j <= n - 1;j++)
				SWAP(a[ i ][ j ], a[ m ][ j ]);
			for(j = 0;j <= n - 1;j++)
				SWAP(a[ j ][ i ], a[ j ][ m ]);
		}
		
		if( x ) // Carry out the elimination.
		{
			for(i = m + 1;i <= n - 1;i++)
			{
				if((y = a[ i ][m - 1]) != 0.0)
				{
					y /= x;
					a[ i ][m - 1] = y;
					for(j = m;j <= n - 1;j++)
						a[ i ][ j ] -= y * a[ m ][ j ];
					for(j = 0;j <= n - 1;j++)
						a[ j ][ m ] += y * a[ j ][ i ];
				}
			}
		}
	}
}

/*
void hqr(rg_REAL **a, rg_INT n, rg_REAL wr[], rg_REAL wi[])
//Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n]. On input a can be
//exactly as output from elmhes x 11.5; on output it is destroyed. The real and imaginary parts
//of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.
{
	rg_INT nn, m, l, k, j, its, i, mmin;
	rg_REAL z, y, x, w, v, u, t, s, r, q, p, anorm;
	anorm = 0.0;	          // Compute matrix rg_MathFunc::norm for possible use in 
	for (i = 1;i <= n;i++)    // locating single small subdiagonal element.
		for(j = IMAX(i - 1, 1);j <= n;j++)
			anorm += fabs(a[ i ][ j ]);
	nn = n;
	t = 0.0;        // Gets changed only by an exceptional shift.
    while (nn >= 1) // Begin search for next eigenvalue.
	{
		its = 0;
		do
		{
			for(l = nn;l >= 2;l--) //Begin iteration: look for single small subdiagonal element.
			{  
				s = fabs(a[l - 1][l - 1]) + fabs(a[ l ][ l ]);
				if (s == 0.0) 
					s = anorm;
				if ((rg_FLOAT)(fabs(a[ l ][l - 1]) + s) == s) 
					break;
			}
			x = a[ nn ][ nn ];

            if(l == nn) // One rg_MathFunc::root found.
			{
				wr[nn] = x + t;
				wi[nn--] = 0.0;
			}
			else
			{
				y = a[nn - 1][nn - 1];
				w = a[nn][nn - 1] * a[nn - 1][nn];
                if (l == (nn - 1)) // Two roots found...
				{
					p = 0.5 * (y - x);
					q = p * p + w;
					z = sqrt(fabs(q));
					x += t;
					if (q >= 0.0) //  ...a real pair.
					{
						z = p + SIGN(z, p);
						wr[nn - 1] = wr[nn] = x + z;
						if( z ) wr[nn] = x - w/z;
						wi[nn - 1] = wi[nn] = 0.0;
					}
					else { // ...a complex pair.
						wr[nn - 1] = wr[nn] = x + p;
						wi[nn - 1]= - (wi[nn] = z);
					}
					nn -= 2;
				}
				else // No roots found. Continue iteration.
				{
					if(its == 30) 
						//nrerror("Too many iterations in hqr"); // ???????????
					if(its == 10 || its == 20)  // Form exceptional shift.
					{
						t += x;
						for(i = 1;i <= nn;i++) a[ i ][ i ] -= x;
						s = fabs(a[nn][nn - 1]) + fabs(a[nn - 1][nn - 2]);
						y = x = 0.75 * s;
						w = -0.4375 * s * s;
					}
					++its;
					for(m = (nn - 2);m >= l;m--) // Form shift and then look for 2 consecutive small subdiagonal elements.
					{
						z = a[ m ][ m ];
						r = x - z;
						s = y - z;
						p = (r * s - w)/a[m + 1][ m ] + a[ m ][m + 1]; // Equation (11.6.23).
						q = a[m + 1][m + 1] - z - r - s;
						r = a[m + 2][m + 1];
						s = fabs(p) + fabs(q) + fabs(r); // Scale to prevent overflow or underflow. 
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						u = fabs(a[ m ][m - 1]) * (fabs(q) + fabs(r));
						v = fabs(p) * (fabs(a[m - 1][m - 1]) + fabs(z) + fabs(a[m + 1][m + 1]));
						if ((rg_FLOAT)(u + v) == v) break; //Equation (11.6.26).
					}
					for(i = m + 2;i <= nn;i++)
					{
						a[ i ][i - 2] = 0.0;
						if (i != (m + 2)) a[ i ][i - 3] = 0.0;
					}
					for(k = m;k <= nn - 1;k++) // Double QR step on rows l to nn and columns m to nn.
					{
						if(k != m)
						{
							p = a[ k ][k - 1]; // Begin setup of Householder vector.
							q = a[k + 1][k - 1];
							r = 0.0;
							if(k != (nn - 1))
								r = a[k + 2][k - 1];
							if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0)
							{
								p /= x; // Scale to prevent over ow orunder ow. 
								q /= x;
								r /= x;
							}
						}
						if((s = SIGN(sqrt(p * p + q * q + r * r), p)) != 0.0)
						{
							if(k == m)
							{
								if(l != m)
									a[ k ][k - 1] = - a[k][k - 1];
							}
							else
								a[ k ][k - 1] = -s*x;
							p += s; // Equations (11.6.24).
							x = p/s;
							y = q/s;
							z = r/s;
							q /= p;
							r /= p;
							for (j = k;j <= nn;j++) // Row modication.
							{
								p = a[ k ][ j ] + q * a[k + 1][j];
								if(k != (nn - 1))
								{
									p += r * a[k + 2][ j ];
									a[k + 2][ j ] -= p * z;
								}
								a[k + 1][ j ] -= p * y;
								a[ k ][ j ] -= p * x;
							}
							mmin = nn< k + 3 ? nn : k + 3;
							for(i = l;i <= mmin;i++) // Column modication.
							{
								p = x * a[ i ][ k ] + y * a [ i ][k + 1];
								if(k != (nn - 1))
								{
									p += z * a[ i ][k + 2];
									a[ i ][k + 2] -= p * r;
								}
								a[ i ][k + 1] -= p * q;
								a[ i ][ k ] -= p;
							}
						}
					}
				}
			}
		}while(l < nn-1);
	}
}
*/

void hqr(rg_REAL **a, rg_INT n, rg_REAL wr[], rg_REAL wi[])
//Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n]. On input a can be
//exactly as output from elmhes x 11.5; on output it is destroyed. The real and imaginary parts
//of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.
{
	rg_INT nn, m, l, k, j, its, i, mmin;
	rg_REAL z, y, x, w, v, u, t, s, r, q, p, anorm;
	anorm = 0.0;	          // Compute matrix rg_MathFunc::norm for possible use in 
	for (i = 0;i <= n - 1;i++)    // locating single small subdiagonal element.
		for(j = IMAX(i - 1, 1);j <= n - 1;j++)
			anorm += fabs(a[ i ][ j ]);
	nn = n - 1;
	t = 0.0;        // Gets changed only by an exceptional shift.
    while (nn >= 0) // Begin search for next eigenvalue.
	{
		its = 0;
		do
		{
			for(l = nn;l >= 1;l--) //Begin iteration: look for single small subdiagonal element.
			{  
				s = fabs(a[l - 1][l - 1]) + fabs(a[ l ][ l ]);
				if (s == 0.0) 
					s = anorm;
				if ((rg_FLOAT)(fabs(a[ l ][l - 1]) + s) == s) 
					break;
			}
			x = a[ nn ][ nn ];

            if(l == nn) // One rg_MathFunc::root found.
			{
				wr[nn] = x + t;
				wi[nn--] = 0.0;
			}
			else
			{
				y = a[nn - 1][nn - 1];
				w = a[nn][nn - 1] * a[nn - 1][nn];
                if (l == (nn - 1)) // Two roots found...
				{
					p = 0.5 * (y - x);
					q = p * p + w;
					z = sqrt(fabs(q));
					x += t;
					if (q >= 0.0) //  ...a real pair.
					{
						z = p + SIGN(z, p);
						wr[nn - 1] = wr[nn] = x + z;
						if( z ) wr[nn] = x - w/z;
						wi[nn - 1] = wi[nn] = 0.0;
					}
					else { // ...a complex pair.
						wr[nn - 1] = wr[nn] = x + p;
						wi[nn - 1]= - (wi[nn] = z);
					}
					nn -= 2;
				}
				else // No roots found. Continue iteration.
				{
					if(its == 30) 
						//nrerror("Too many iterations in hqr"); // ???????????
					if(its == 10 || its == 20)  // Form exceptional shift.
					{
						t += x;
						for(i = 0;i <= nn;i++) a[ i ][ i ] -= x;
						s = fabs(a[nn][nn - 1]) + fabs(a[nn - 1][nn - 2]);
						y = x = 0.75 * s;
						w = -0.4375 * s * s;
					}
					++its;
					for(m = (nn - 2);m >= l;m--) // Form shift and then look for 2 consecutive small subdiagonal elements.
					{
						z = a[ m ][ m ];
						r = x - z;
						s = y - z;
						p = (r * s - w)/a[m + 1][ m ] + a[ m ][m + 1]; // Equation (11.6.23).
						q = a[m + 1][m + 1] - z - r - s;
						r = a[m + 2][m + 1];
						s = fabs(p) + fabs(q) + fabs(r); // Scale to prevent overflow or underflow. 
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						u = fabs(a[ m ][m - 1]) * (fabs(q) + fabs(r));
						v = fabs(p) * (fabs(a[m - 1][m - 1]) + fabs(z) + fabs(a[m + 1][m + 1]));
						if ((rg_FLOAT)(u + v) == v) break; //Equation (11.6.26).
					}
					for(i = m + 2;i <= nn;i++)
					{
						a[ i ][i - 2] = 0.0;
						if (i != (m + 2)) a[ i ][i - 3] = 0.0;
					}
					for(k = m;k <= nn - 1;k++) // Double QR step on rows l to nn and columns m to nn.
					{
						if(k != m)
						{
							p = a[ k ][k - 1]; // Begin setup of Householder vector.
							q = a[k + 1][k - 1];
							r = 0.0;
							if(k != (nn - 1))
								r = a[k + 2][k - 1];
							if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0)
							{
								p /= x; // Scale to prevent over ow orunder ow. 
								q /= x;
								r /= x;
							}
						}
						if((s = SIGN(sqrt(p * p + q * q + r * r), p)) != 0.0)
						{
							if(k == m)
							{
								if(l != m)
									a[ k ][k - 1] = - a[k][k - 1];
							}
							else
								a[ k ][k - 1] = - s * x;
							p += s; // Equations (11.6.24).
							x = p/s;
							y = q/s;
							z = r/s;
							q /= p;
							r /= p;
							for (j = k;j <= nn;j++) // Row modication.
							{
								p = a[ k ][ j ] + q * a[k + 1][j];
								if(k != (nn - 1))
								{
									p += r * a[k + 2][ j ];
									a[k + 2][ j ] -= p * z;
								}
								a[k + 1][ j ] -= p * y;
								a[ k ][ j ] -= p * x;
							}
							mmin = nn < k + 3 ? nn : k + 3;
							for(i = l;i <= mmin;i++) // Column modication.
							{
								p = x * a[ i ][ k ] + y * a [ i ][k + 1];
								if(k != (nn - 1))
								{
									p += z * a[ i ][k + 2];
									a[ i ][k + 2] -= p * r;
								}
								a[ i ][k + 1] -= p * q;
								a[ i ][ k ] -= p;
							}
						}
					}
				}
			}
		}while(l < nn - 1);
	}
}


