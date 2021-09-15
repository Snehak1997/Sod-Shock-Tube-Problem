# include<iostream>
# include<cmath>
using namespace std;
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
int main()
{
	int i, j, iter = 0;
	int N = 402;								// No. of cells + 2(ghost cells)
	double rho_L = 1.0, rho_R = 0.125;			// initial densities
	double u_L = 0.0, u_R = 0.0;				// initial velocities
	double P_L = 1.0, P_R = 0.1;				// initial pressures
	double L = 1.0, t=0.0, t_end = 0.1;								// length of tube
	double gamma = 1.4, C=0.8;
	double dx, dt, delt[N], L_max[N];
	dx = L/(N-2);
	double x[N], x_half = dx/2;
	double rho[N], u[N], P[N], e[N], H[N], IE[N], a[N];
	double U1[N], U2[N], U3[N];
	double L1[N], L2[N], L3[N];					// eigenvalues
	double K1[3], K2[3], K3[3], det_K;					// eigenvectors
	double K1_inv[3], K2_inv[3], K3_inv[3];
	double L1_plus, L1_minus;
	double L2_plus, L2_minus;
	double L3_plus, L3_minus;
	double A1_plus[3], A1_minus[3];
	double A2_plus[3], A2_minus[3];
	double A3_plus[3], A3_minus[3];
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
	//	cout<<i<<" "<<x[i]<<" "<<rho[i]<<" "<<u[i]<<" "<<P[i]<<endl;
	}	
	rho[0] = rho[1];
	u[0] = -u[1];
	P[0] = P[1];
	e[0] = 0.5*pow(u[0],2) + P[0]/((gamma-1)*rho[0]);
	rho[401] = rho[400];
	u[401] = -u[400];
	P[401] = P[400];
	e[401] = 0.5*pow(u[401],2) + P[401]/((gamma-1)*rho[401]);
	while(t<t_end)
	{
		iter++;
		for(i=0; i<N; i++)
		{
			U1[i] = rho[i];
			U2[i] = rho[i]*u[i];
			U3[i] = rho[i]*e[i];
		//	cout<<i<<" "<<U1[i]<<" "<<U2[i]<<" "<<U3[i]<<endl;
		}
		
		for(i=1; i<N-1; i++)
		{
			L1[i] = u[i]-a[i];
			L2[i] = u[i];
			L3[i] = u[i]+a[i];
			L_max[i] = fabs(u[i])+a[i];
		//	cout<<L1[i]<<" "<<L2[i]<<" "<<L3[i]<<endl;
		}
		dt = C*(dx/maximum(L_max, N));
	//	cout<<"max lambda "<<maximum(L_max, N);
		for(i=1; i<N-1; i++)
		{
			L1_plus = max(L1[i],0);
			L1_minus = min(L1[i],0);
			L2_plus = max(L2[i],0);
			L2_minus = min(L2[i],0);
			L3_plus = max(L3[i],0);
			L3_minus = min(L3[i],0);
			
			K1[0] = 1;
			K1[1] = u[i]-a[i];
			K1[2] = H[i]-u[i]*a[i];
			
			K2[0] = 1;
			K2[1] = u[i];
			K2[2] = 0.5*u[i]*u[i];
			
			K3[0] = 1;
			K3[1] = u[i]+a[i];
			K3[2] = H[i]+u[i]*a[i];
			
			det_K = 2*a[i]*H[i] - u[i]*u[i]*a[i];
			
			K1_inv[0] = (1/det_K)*(u[i]*H[i] + 0.5*u[i]*u[i]*a[i] - 0.5*pow(u[i],3));
			K1_inv[1] = (1/det_K)*(2*(a[i]*H[i] - a[i]*u[i]*u[i]));
			K1_inv[2] = (1/det_K)*(-u[i]*H[i] + 0.5*u[i]*u[i]*a[i] + 0.5*pow(u[i],3));
			
			K2_inv[0] = (1/det_K)*(0.5*u[i]*u[i] - H[i] - u[i]*a[i]);
			K2_inv[1] = (1/det_K)*(2*u[i]*a[i]);
			K2_inv[2] = (1/det_K)*(-0.5*u[i]*u[i] + H[i] - u[i]*a[i]);
			
			K3_inv[0] = (1/det_K)*a[i];
			K3_inv[1] = -(1/det_K)*2*a[i];
			K3_inv[2] = (1/det_K)*a[i];
			
			A1_plus[0] = K1[0]*L1_plus*K1_inv[0] + K2[0]*L2_plus*K1_inv[1] + K3[0]*L3_plus*K1_inv[2];	
			A1_plus[1] = K1[0]*L1_plus*K2_inv[0] + K2[0]*L2_plus*K2_inv[1] + K3[0]*L3_plus*K2_inv[2];	
			A1_plus[2] = K1[0]*L1_plus*K3_inv[0] + K2[0]*L2_plus*K3_inv[1] + K3[0]*L3_plus*K3_inv[2];
			
			A2_plus[0] = K1[1]*L1_plus*K1_inv[0] + K2[1]*L2_plus*K1_inv[1] + K3[1]*L3_plus*K1_inv[2];	
			A2_plus[1] = K1[1]*L1_plus*K2_inv[0] + K2[1]*L2_plus*K2_inv[1] + K3[1]*L3_plus*K2_inv[2];	
			A2_plus[2] = K1[1]*L1_plus*K3_inv[0] + K2[1]*L2_plus*K3_inv[1] + K3[1]*L3_plus*K3_inv[2];
		
			A3_plus[0] = K1[2]*L1_plus*K1_inv[0] + K2[2]*L2_plus*K1_inv[1] + K3[2]*L3_plus*K1_inv[2];	
			A3_plus[1] = K1[2]*L1_plus*K2_inv[0] + K2[2]*L2_plus*K2_inv[1] + K3[2]*L3_plus*K2_inv[2];	
			A3_plus[2] = K1[2]*L1_plus*K3_inv[0] + K2[2]*L2_plus*K3_inv[1] + K3[2]*L3_plus*K3_inv[2];
		
			A1_minus[0] = K1[0]*L1_minus*K1_inv[0] + K2[0]*L2_minus*K1_inv[1] + K3[0]*L3_minus*K1_inv[2];	
			A1_minus[1] = K1[0]*L1_minus*K2_inv[0] + K2[0]*L2_minus*K2_inv[1] + K3[0]*L3_minus*K2_inv[2];	
			A1_minus[2] = K1[0]*L1_minus*K3_inv[0] + K2[0]*L2_minus*K3_inv[1] + K3[0]*L3_minus*K3_inv[2];
			
			A2_minus[0] = K1[1]*L1_minus*K1_inv[0] + K2[1]*L2_minus*K1_inv[1] + K3[1]*L3_minus*K1_inv[2];	
			A2_minus[1] = K1[1]*L1_minus*K2_inv[0] + K2[1]*L2_minus*K2_inv[1] + K3[1]*L3_minus*K2_inv[2];	
			A2_minus[2] = K1[1]*L1_minus*K3_inv[0] + K2[1]*L2_minus*K3_inv[1] + K3[1]*L3_minus*K3_inv[2];
		
			A3_minus[0] = K1[2]*L1_minus*K1_inv[0] + K2[2]*L2_minus*K1_inv[1] + K3[2]*L3_minus*K1_inv[2];	
			A3_minus[1] = K1[2]*L1_minus*K2_inv[0] + K2[2]*L2_minus*K2_inv[1] + K3[2]*L3_minus*K2_inv[2];	
			A3_minus[2] = K1[2]*L1_minus*K3_inv[0] + K2[2]*L2_minus*K3_inv[1] + K3[2]*L3_minus*K3_inv[2];
		
			U1[i] = U1[i]-(dt/dx)*(A1_plus[0]*(U1[i]-U1[i-1]) + A1_plus[1]*(U2[i]-U2[i-1]) + A1_plus[2]*(U3[i]-U3[i-1]));
			U1[i] = U1[i]-(dt/dx)*(A1_minus[0]*(U1[i]-U1[i-1]) + A1_minus[1]*(U2[i]-U2[i-1]) + A1_minus[2]*(U3[i]-U3[i-1]));
			
			U2[i] = U2[i]-(dt/dx)*(A2_plus[0]*(U1[i]-U1[i-1]) + A2_plus[1]*(U2[i]-U2[i-1]) + A2_plus[2]*(U3[i]-U3[i-1]));
			U2[i] = U2[i]-(dt/dx)*(A2_minus[0]*(U1[i]-U1[i-1]) + A2_minus[1]*(U2[i]-U2[i-1]) + A2_minus[2]*(U3[i]-U3[i-1]));
		
			U3[i] = U3[i]-(dt/dx)*(A3_plus[0]*(U1[i]-U1[i-1]) + A3_plus[1]*(U2[i]-U2[i-1]) + A3_plus[2]*(U3[i]-U3[i-1]));
			U3[i] = U3[i]-(dt/dx)*(A3_minus[0]*(U1[i]-U1[i-1]) + A3_minus[1]*(U2[i]-U2[i-1]) + A3_minus[2]*(U3[i]-U3[i-1]));
		}
		for(i=1; i<N-1; i++)
		{
			rho[i] = U1[i];
			u[i] = U2[i]/rho[i];
			e[i] = U3[i]/rho[i];
			P[i] = rho[i]*(gamma-1)*(e[i] - 0.5*u[i]*u[i]);
			H[i] = (rho[i]*e[i] + P[i])/rho[i];	
			a[i] = sqrt(gamma*(P[i]/rho[i]));	
		}
		rho[0] = rho[1];
		u[0] = -u[1];
		P[0] = P[1];
		e[0] = 0.5*pow(u[0],2) + P[0]/((gamma-1)*rho[0]);
		rho[401] = rho[400];
		u[401] = -u[400];
		P[401] = P[400];
		e[401] = 0.5*pow(u[401],2) + P[401]/((gamma-1)*rho[401]);
		t = t + dt;
	}
	for(i=1; i<N-1; i++)
	{
		cout<<i<<" "<<rho[i]<<" "<<u[i]<<" "<<P[i]<<endl;
	}
	return 0;	
}


