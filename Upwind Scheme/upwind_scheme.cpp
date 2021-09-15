# include<iostream>
# include<cmath>
# include<fstream>
using namespace std;
// to find maximum in an array
double maximum(double A[], int size)
{
	int i;
	double t;
	t=A[1];
	for(i=2;i<size-1;i++)
    {
		if(A[i]>t)
		t=A[i];
	}
	return(t);
}
// First term of array
double first(double* U)
{
	double a;
	a = U[0];
	return a;
}
// Second term of array
double second(double* U)
{
	double a;
	a = U[1];
	return a;
}
// Third term of array
double third(double* U)
{
	double a;
	a = U[2];
	return a;
}
// maximum of two numbers
double max(double a, double b)
{
	double m;
	if(a>b)
	{
		m = a;
	}
	else
	{
		m = b;
	}
	return m;
}
// minimum of two numbers
double min(double a,double b)
{
	double m;
	if(a<b)
	{
		m = a;
	}
	else
	{
		m = b;
	}
	return m;
}
// set of 3 numbers
double* vector(double a, double b, double c)
{
	double* A = new double[3];
	A[0] = a;
	A[1] = b;
	A[2] = c;
	return A;
}
// Function of lambda(+)
double* max_lambda(double* L)
{
	double* B = new double[3];
	for(int i=0; i<3; i++)
	{
		B[i]= max(L[i],0);
	}
	return B; 
}
// Function of lambda(-)
double* min_lambda(double* L)
{
	double* B = new double[3];
	for(int i=0; i<3; i++)
	{
		B[i]= min(L[i],0);
	}
	return B;
}
// Function of matrix K
double** K_matrix(double a, double u, double H)
{
	double** K = new double*[3];
	for(int i=0;i<3;i++)
	{
		K[i] = new double[3];
	}
	K[0][0] = 1;
	K[1][0] = u-a;
	K[2][0] = H-u*a;
	
	K[0][1] = 1;
	K[1][1] = u;
	K[2][1] = 0.5*u*u;
	
	K[0][2] = 1;
	K[1][2] = u+a;
	K[2][2] = H+u*a;
	return K;
}
// Function of matrix K inverse
double** K_inv(double a, double u, double H)
{
	double** K = new double*[3];
	for(int i=0;i<3;i++)
	{
		K[i] = new double[3];	
	}	
	double det_K = 2*a*H - u*u*a;
	K[0][0] = (1/det_K)*(u*H + 0.5*u*u*a - 0.5*u*u*u);
	K[1][0] = (1/det_K)*2*(a*H - a*u*u);
	K[2][0] = (1/det_K)*(-u*H + 0.5*u*u*a + 0.5*u*u*u);
			
	K[0][1] = (1/det_K)*(0.5*u*u - H - u*a);
	K[1][1] = (1/det_K)*2*u*a;
	K[2][1] = (1/det_K)*(-0.5*u*u + H - u*a);
			
	K[0][2] = (1/det_K)*a;
	K[1][2] = -(1/det_K)*2*a;
	K[2][2] = (1/det_K)*a;
	return K;
}
// Function of matrix A
double** A_matrix(double* L, double** K, double** K_inv)
{
	double** R = new double*[3];
	for(int i=0; i<3; i++)
	{
		R[i] = new double[3];
	}
	double** A = new double*[3];
	for(int i=0; i<3; i++)
	{
		A[i] = new double[3];
	}
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			R[i][j] = L[i]*K_inv[i][j];
		}
	}
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			A[i][j] = 0;
			for(int k=0; k<3; k++)
			{
				A[i][j] += K[i][k]*R[k][j];	
			}			
		}
	}
	return A;
}
// matrix vector multiplication
double* mat_vec(double* A, double** B)
{
	double* C = new double[3];
	for(int i=0; i<3; i++)
	{
		C[i] = 0;
	}
	for(int i=0;i<3;i++)
	{
		for(int j=0; j<3; j++)
		{		
			{
				C[i] += B[i][j]*A[j];	
			}			
		}
	}
	return C;
}
// Function for calculation of flux
double* Flux(double* U1, double* U2, double** A)
{
	double* R = new double[3];
	double* S = new double[3];
	for(int i=0; i<3; i++)
	{
		R[i] = U1[i] - U2[i];
	}
	S = mat_vec(R,A); 
	return S;
}

int main()
{
	int i, j, iter = 0;
	int N = 402;								// No. of cells + 2(ghost cells)
	double rho_L = 1.0, rho_R = 0.125;			// initial densities
	double u_L = 0.0, u_R = 0.0;				// initial velocities
	double P_L = 1.0, P_R = 0.1;				// initial pressures
	double L = 1.0;								// length of tube
	double t=0.0, t_end = 0.1;					// time domain
	double gamma = 1.4, C=0.8;
	double dx, dt, delt[N], L_max[N];
	dx = L/(N-2);
	double x[N], x_half = dx/2;
	double rho[N], u[N], P[N], e[N], H[N], IE[N], a[N];
	double U1[N], U2[N], U3[N];
	double Un1[N], Un2[N], Un3[N];
	double F11[N], F12[N], F21[N], F22[N], F31[N], F32[N];
	double L1[N], L2[N], L3[N];
	double* U[N]; double* U_new[N];
	double* F1[N]; double* F2[N];					
	double* lambda[N]; double* Lmax[N]; double* Lmin[N];
	double** K_mat[N]; double** K_inverse[N]; double** A_plus[i]; double** A_minus[N];
	for(i=1; i<N-1; i++)
	{
		x[i] = -x_half + i*dx;
		if(x[i]<0.5)
		{
			rho[i] = rho_L;
			u[i] = u_L;
			P[i] = P_L;
		}
		else
		{
			rho[i] = rho_R;
			u[i] = u_R;
			P[i] = P_R;
		}
		e[i] = 0.5*pow(u[i],2) + P[i]/((gamma-1)*rho[i]);
		H[i] = (rho[i]*e[i] + P[i])/rho[i];	
		a[i] = sqrt(gamma*(P[i]/rho[i]));
		IE[i] = P[i]/((gamma-1)*rho[i]);
	}	
	rho[0] = rho[1];
	u[0] = -u[1];
	P[0] = P[1];
	e[0] = e[1];
	rho[N-1] = rho[N-2];
	u[N-1] = -u[N-2];
	P[N-1] = P[N-2];
	e[N-1] = e[N-2];
	
	while(t<t_end)
	{
		iter++;
		for(i=0; i<N; i++)
		{
		U1[i] = rho[i];
		U2[i] = rho[i]*u[i];
		U3[i] = rho[i]*e[i];
		}
		for(i=0; i<N; i++)
		{
			U[i] = vector(U1[i], U2[i], U3[i]);
		}
		for(i=0; i<N; i++)
		{
			U_new[i] = U[i];
		}
		// calculation of fluxes
		for(i=1; i<N-1; i++)
		{
			lambda[i] = vector(u[i]-a[i],u[i],u[i]+a[i]);
			L_max[i] = fabs(u[i])+a[i];
			Lmax[i] = max_lambda(lambda[i]);
			Lmin[i] = min_lambda(lambda[i]);
			K_mat[i] = K_matrix(a[i], u[i], H[i]);
			K_inverse[i] = K_inv(a[i], u[i], H[i]);
			A_plus[i] = A_matrix(Lmax[i], K_mat[i] ,K_inverse[i]);
			A_minus[i] = A_matrix(Lmin[i], K_mat[i], K_inverse[i]);
			F1[i] = Flux(U_new[i], U_new[i-1], A_plus[i]);
			F2[i] = Flux(U_new[i+1], U_new[i], A_minus[i]);	
 		}
 		dt = C*dx/maximum(L_max,N);
 		// Upwind Scheme
 		for(i=1; i<N-1; i++)
 		{
			Un1[i] = first(U_new[i]); Un2[i] = second(U_new[i]); Un3[i] = third(U_new[i]);
			F11[i] = first(F1[i]); F12[i] = first(F2[i]);
			F21[i] = second(F1[i]); F22[i] = second(F2[i]);
			F31[i] = third(F1[i]); F32[i] = third(F2[i]);
			Un1[i] = Un1[i] - (dt/dx)*(F11[i]+F12[i]);
			Un2[i] = Un2[i] - (dt/dx)*(F21[i]+F22[i]);
			Un3[i] = Un3[i] - (dt/dx)*(F31[i]+F32[i]);
			U1[i] = Un1[i];
			U2[i] = Un2[i];
			U3[i] = Un3[i];
			rho[i] = U1[i];
			u[i] = U2[i]/U1[i];
			e[i] = U3[i]/U1[i];
			P[i] = (e[i] - 0.5*u[i]*u[i])*rho[i]*(gamma-1);
			H[i] = (rho[i]*e[i] + P[i])/rho[i];	
			a[i] = sqrt(gamma*(P[i]/rho[i]));
			IE[i] = P[i]/((gamma-1)*rho[i]);
		}
		rho[0] = rho[1];
		u[0] = -u[1];
		P[0] = P[1];
		e[0] = e[1];
		rho[N-1] = rho[N-2];
		u[N-1] = -u[N-2];
		P[N-1] = P[N-2];
		e[N-1] = e[N-2];
		t += dt; 	
	}
	ofstream myfile("1DShocktube(upwind_scheme).txt"); 
	for(i=1; i<N-1; i++)
	{
		cout<<x[i]<<" "<<rho[i]<<" "<<u[i]<<" "<<P[i]<<" "<<IE[i]<<endl;
		myfile<<x[i]<<" "<<rho[i]<<" "<<u[i]<<" "<<P[i]<<" "<<IE[i]<<endl;
	}
	return 0;
}
