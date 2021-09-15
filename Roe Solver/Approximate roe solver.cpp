# include<iostream>
# include<fstream>
# include<cmath>
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
	double dx, dt, L_max[N];
	dx = L/(N-2);
	double x[N], x_half = dx/2;

	double rho[N], u[N], P[N], e[N], H[N], IE[N], a[N];

	double U1[N], U2[N], U3[N];
	double Q1[N], Q2[N], Q3[N];
	double F1[N], F2[N], F3[N];
	double Fl_1[N], Fl_2[N], Fl_3[N];
	double Fr_1[N], Fr_2[N], Fr_3[N];
	double Un1[N], Un2[N], Un3[N];

	double dU1[N], dU2[N], dU3[N];

	double a_[N], u_[N], H_[N];
	double L1_[N], L2_[N], L3_[N];
	double K1[3][N], K2[3][N], K3[3][N];
    double del_1[N],del_2[N], del_3[N];

	// initial conditions
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

	//boundary condition
	rho[0] = rho[1];
	u[0] = -u[1];
	P[0] = P[1];
	e[0] = e[1];
	H[0] = H[1];
	rho[N-1] = rho[N-2];
	u[N-1] = -u[N-2];
	P[N-1] = P[N-2];
	e[N-1] = e[N-2];
	H[N-1] = H[N-2];

	while(t<t_end)
    {
        // updating U
        for(i=0; i<N; i++)
        {
            U1[i] = rho[i];
            U2[i] = U1[i]*u[i];
            U3[i] = U1[i]*e[i];

            Q1[i] = sqrt(rho[i]);
            Q2[i] = Q1[i]*u[i];
            Q3[i] = Q1[i]*H[i];

            F1[i] = U2[i];
            F2[i] = (0.5*(3-gamma)*U2[i]*U2[i])/U1[i] + (gamma-1)*U3[i];
            F3[i] = (gamma*U2[i]*U3[i])/U1[i] - 0.5*(gamma-1)*((U2[i]*U2[i]*U2[i])/(U1[i]*U1[i]));

            Un1[i] = U1[i];
            Un2[i] = U2[i];
            Un3[i] = U3[i];
        }

        for(i=1; i<N-1; i++)
        {
            L_max[i] = fabs(u[i]) + a[i];
        }

        dt = C*dx/maximum(L_max, N);

        // flux at right face
        for(i=1; i<N-1; i++)
        {
            u_[i] = (Q2[i]+Q2[i+1])/(Q1[i]+Q1[i+1]);
            H_[i] = (Q3[i]+Q3[i+1])/(Q1[i]+Q1[i+1]);
            a_[i] = (gamma-1)*pow((H_[i]-0.5*(u_[i]*u_[i])),0.5);

            L1_[i] = u_[i] - a_[i];
            L2_[i] = u_[i];
            L3_[i] = u_[i] + a_[i];

            K1[0][i] = 1; K1[1][i] = u_[i] - a_[i]; K1[2][i] = H_[i] - u_[i]*a_[i];
            K2[0][i] = 1; K2[1][i] = u_[i]; K2[2][i] = 0.5*u_[i]*u_[i];
            K3[0][i] = 1; K3[1][i] = u_[i] + a_[i]; K3[2][i] = H_[i] + u_[i]*a_[i];

            dU1[i] = Un1[i+1] - Un1[i];
            dU2[i] = Un2[i+1] - Un2[i];
            dU3[i] = Un3[i+1] - Un3[i];

            del_2[i] = ((gamma-1)/(a_[i]*a_[i]))*(dU1[i]*(H_[i] - u_[i]*u_[i]) + u_[i]*dU2[i] - dU3[i]);
            del_1[i] = (0.5/a_[i])*(dU1[i]*(u_[i]+a_[i]) - dU2[i] - a_[i]*del_2[i]);
            del_3[i] = dU1[i] - (del_1[i] + del_2[i]);

            Fr_1[i] = 0.5*(F1[i] + F1[i+1]) - 0.5*(fabs(L1_[i])*del_1[i]*K1[0][i] + fabs(L2_[i])*del_2[i]*K2[0][i] + fabs(L3_[i])*del_3[i]*K3[0][i]);
            Fr_2[i] = 0.5*(F2[i] + F2[i+1]) - 0.5*(fabs(L1_[i])*del_1[i]*K1[1][i] + fabs(L2_[i])*del_2[i]*K2[1][i] + fabs(L3_[i])*del_3[i]*K3[1][i]);
            Fr_3[i] = 0.5*(F3[i] + F3[i+1]) - 0.5*(fabs(L1_[i])*del_1[i]*K1[2][i] + fabs(L2_[i])*del_2[i]*K2[2][i] + fabs(L3_[i])*del_3[i]*K3[2][i]);
        }
        // for left face
        for(i=1; i<N-1; i++)
        {
            u_[i] = (Q2[i]+Q2[i-1])/(Q1[i]+Q1[i-1]);
            H_[i] = (Q3[i]+Q3[i-1])/(Q1[i]+Q1[i-1]);
            a_[i] = (gamma-1)*pow((H_[i]-0.5*(u_[i]*u_[i])),0.5);

            L1_[i] = u_[i] - a_[i];
            L2_[i] = u_[i];
            L3_[i] = u_[i] + a_[i];

            K1[0][i] = 1; K1[1][i] = u_[i] - a_[i]; K1[2][i] = H_[i] - u_[i]*a_[i];
            K2[0][i] = 1; K2[1][i] = u_[i]; K2[2][i] = 0.5*u_[i]*u_[i];
            K3[0][i] = 1; K3[1][i] = u_[i] + a_[i]; K3[2][i] = H_[i] + u_[i]*a_[i];

            dU1[i] = Un1[i] - Un1[i-1];
            dU2[i] = Un2[i] - Un2[i-1];
            dU3[i] = Un3[i] - Un3[i-1];

            del_2[i] = ((gamma-1)/(a_[i]*a_[i]))*(dU1[i]*(H_[i] - u_[i]*u_[i]) + u_[i]*dU2[i] - dU3[i]);
            del_1[i] = (0.5/a_[i])*(dU1[i]*(u_[i]+a_[i]) - dU2[i] - a_[i]*del_2[i]);
            del_3[i] = dU1[i] - (del_1[i] + del_2[i]);

            Fl_1[i] = 0.5*(F1[i] + F1[i-1]) - 0.5*(fabs(L1_[i])*del_1[i]*K1[0][i] + fabs(L2_[i])*del_2[i]*K2[0][i] + fabs(L3_[i])*del_3[i]*K3[0][i]);
            Fl_2[i] = 0.5*(F2[i] + F2[i-1]) - 0.5*(fabs(L1_[i])*del_1[i]*K1[1][i] + fabs(L2_[i])*del_2[i]*K2[1][i] + fabs(L3_[i])*del_3[i]*K3[1][i]);
            Fl_3[i] = 0.5*(F3[i] + F3[i-1]) - 0.5*(fabs(L1_[i])*del_1[i]*K1[2][i] + fabs(L2_[i])*del_2[i]*K2[2][i] + fabs(L3_[i])*del_3[i]*K3[2][i]);
        }
        // calculating U
        for(i=1; i<N-1; i++)
        {
            Un1[i] = Un1[i] + (dt/dx)*(Fl_1[i] - Fr_1[i]);
            Un2[i] = Un2[i] + (dt/dx)*(Fl_2[i] - Fr_2[i]);
            Un3[i] = Un3[i] + (dt/dx)*(Fl_3[i] - Fr_3[i]);

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
        // updating boundary conditions
        rho[0] = rho[1];
		u[0] = -u[1];
		P[0] = P[1];
		e[0] = e[1];
		H[0] = H[1];
		rho[N-1] = rho[N-2];
		u[N-1] = -u[N-2];
		P[N-1] = P[N-2];
		e[N-1] = e[N-2];
        H[N-1] = H[N-2];
        t += dt;
    }
    ofstream myfile("1DShockTube(Approximate roe solver).txt");
    for(i=1; i<N-1; i++)
    {
        cout<<i<<" "<<rho[i]<<" "<<u[i]<<" "<<P[i]<<" "<<IE[i]<<endl;
        myfile<<x[i]<<"\t"<<rho[i]<<"\t"<<u[i]<<"\t"<<P[i]<<"\t"<<IE[i]<<endl;
    }
	return 0;
}
