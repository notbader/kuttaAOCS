#include <stdio.h>
#include <math.h>

#define TOTAL_TIME		5
#define TIMESTEP		0.01
#define MAXITERATIONS	10000

int n = 0, n1 = 0;
double I1 = 1.9, I2 = 2.2, I3 = 2.0, h, h1;
double x[MAXITERATIONS], y[MAXITERATIONS], z[MAXITERATIONS];
double Psi[MAXITERATIONS], Theta[MAXITERATIONS];
double Phi[MAXITERATIONS];
double Phi0, Theta0, Psi0;
double L, M, N;
double x0, y01, z0, t0, tn = 0.0, tfinal;
double two_pi;

double f(double t, double w1, double w2, double w3)
{
	double temp;
	
	temp = (L + (I2 - I3) * w2 * w3)/I1;
	return temp;
	
}

double g(double t, double w1, double w2, double w3)
{
	double temp;
	
	temp = (M + (I3 - I1) * w1 * w3)/I2;
	return temp;

}

double v(double t, double w1, double w2, double w3)
{
	double temp;
	
	temp = (N + (I1 - I2) * w1 * w2)/I3;
	return temp;
}

void runge(double t0, double x0, double y0,double z0, double tn, double h)
{
	int times, i;
	double kn1, kn2, kn3, kn4, ln1, ln2, ln3, ln4, mn1, mn2, mn3, mn4;
	
	x[n] = x0;
	y[n] = y0;
	z[n] = z0;
	
	times = (tn - t0)/h;
	
	for (i = 0; i < times; i++){
		
		kn1 = f(t0+n*h,x[n],y[n],z[n]);
		ln1 = g(t0+n*h,x[n],y[n],z[n]);
		mn1 = v(t0+n*h,x[n],y[n],z[n]); 
		
		kn2 = f(t0+n*h+1.0/2.0*h, x[n]+1.0/2.0*h*kn1,y[n]+1.0/2.0*h*ln1, z[n]+1.0/2.0*h*mn1);
		ln2 = g(t0+n*h+1.0/2.0*h, x[n]+1.0/2.0*h*kn1,y[n]+1.0/2.0*h*ln1, z[n]+1.0/2.0*h*mn1);
		mn2 = v(t0+n*h+1.0/2.0*h, x[n]+1.0/2.0*h*kn1,y[n]+1.0/2.0*h*ln1, z[n]+1.0/2.0*h*mn1);
		kn3 = f(t0+n*h+1.0/2.0*h, x[n]+1.0/2.0*h*kn2,y[n]+1.0/2.0*h*ln2, z[n]+1.0/2.0*h*mn2);
		ln3 = g(t0+n*h+1.0/2.0*h, x[n]+1.0/2.0*h*kn2,y[n]+1.0/2.0*h*ln2, z[n]+1.0/2.0*h*mn2);
		mn3 = v(t0+n*h+1.0/2.0*h, x[n]+1.0/2.0*h*kn2,y[n]+1.0/2.0*h*ln2, z[n]+1.0/2.0*h*mn2);
		kn4 = f(t0　+　n　*　h　+　h, x[n]　+　h　*　kn3, y[n]　+　h　*　ln3, z[n] + h * mn3);
		ln4 = g(t0　+　n　*　h　+　h, x[n]　+　h　*　kn3, y[n]　+　h　*　ln3, z[n] + h * mn3);
		mn4 = v(t0　+　n　*　h　+　h, x[n]　+　h　*　kn3, y[n]　+　h　*　ln3, z[n] + h　*　mn3);
		
		x[n+1] = x[n] + 1.0/6.0 * h * (kn1 + 2.0 * kn2 + 2.0 * kn3 + kn4);
		y[n+1] = y[n] + 1.0/6.0 * h * (ln1 + 2.0 * ln2 + 2.0 * ln3 + ln4);
		z[n+1] = z[n] + 1.0/6.0 * h * (mn1 + 2.0 * mn2 + 2.0 * mn3 + mn4);
		
		++n;
		
	}
	
}
	double f1(double t, double theta, double phi, double psi)
	{
		double temp;
		int index;
		
		index = (t - t0)/h;
		temp = y[index] * cos(phi) - z[index] * sin(phi);
		return temp;
	}
	
	double g1(double t, double theta, double phi, double psi)
	{
		double temp;
		int index;
		
		index = (t - t0)/h;
		temp = x [index] + y[index]*sin(phi)*tan(theta) + z[index]*cos(phi)*tan(theta);
		return temp;
		
	}
	
	double v1(double t, double theta, double phi, double psi)
	{
		double temp;
		int index;
		
		index = (t - t0)/h;
		temp = (y[index]*sin(phi) + z[index]*cos(phi))/cos(theta);
		return temp;
	
	}
	
	void runge1(double t0, double theta0, double phi0, double psi0, double tn, double h)
	{
		int times, i, round_temp;
		double kn1, kn2, kn3, kn4, ln1, ln2, ln3, ln4, mn1, mn2, mn3, mn4;
		Theta[n1] = theta0;
		Phi[n1] = phi0;
		Psi[n1] = psi0;
		times = (tn-t0)/h;
		
		for (i = 0 ; i < times ; i++) {
			kn1 = f1(t0+n1*h, Theta[n1], Phi[n1], Psi[n1]);
			ln1 = g1(t0+n1*h, Theta[n1], Phi[n1], Psi[n1]);
			mn1 = v1(t0+n1*h, Theta[n1], Phi[n1], Psi[n1]);
			
			kn2 = f1(t0+n1*h+1.0/2.0*h , Theta[n1]+1.0/2.0*h*kn1 , Phi[n1]+1.0/2.0*h*ln1 , Psi[n1]+1.0/2.0*h*mn1);
			ln2 = g1(t0+n1*h+1.0/2.0*h , Theta[n1]+1.0/2.0*h*kn1 , Phi[n1]+1.0/2.0*h*ln1 , Psi[n1]+1.0/2.0*h*mn1);
			mn2 = v1(t0+n1*h+1.0/2.0*h , Theta[n1]+1.0/2.0*h*kn1 , Phi[n1]+1.0/2.0*h*ln1 , Psi[n1]+1.0/2.0*h*mn1);
			kn3 = f1(t0+n1*h+1.0/2.0*h , Theta[n1]+1.0/2.0*h*kn2 , Phi[n1]+1.0/2.0*h*ln2 , Psi[n1]+1.0/2.0*h*mn2);
			ln3 = g1(t0+n1*h+1.0/2.0*h , Theta[n1]+1.0/2.0*h*kn2 , Phi[n1]+1.0/2.0*h*ln2 , Psi[n1]+1.0/2.0*h*mn2);
			mn3 = v1(t0+n1*h+1.0/2.0*h , Theta[n1]+1.0/2.0*h*kn2 , Phi[n1]+1.0/2.0*h*ln2 , Psi[n1]+1.0/2.0*h*mn2);
			kn4 = f1(t0+n1*h+h , Theta[n1]+h*kn3 , Phi[n1]+h*ln3 , Psi[n1]+h*mn3);
			ln4 = g1(t0+n1*h+h , Theta[n1]+h*kn3 , Phi[n1]+h*ln3 , Psi[n1]+h*mn3);
			mn4 = v1(t0+n1*h+h , Theta[n1]+h*kn3 , Phi[n1]+h*ln3 , Psi[n1]+h*mn3);
			
			Theta[n1+1]=Theta[n1] + 1.0/6.0*h*(kn1 + 2.0*kn2 + 2.0*kn3+kn4);
			Phi[n1+1]=Phi[n1]+1.0/6.0*h*(ln1+2.0*ln2+2.0*ln3+ln4);
			Psi[n1+1]=Psi[n1]+1.0/6.0*h*(mn1+2.0*mn2+2.0*mn3+mn4);
			
			++n1;
			
			//if (Psi[n1] > two_pi) {
				//round_temp = Psi[n1]/two_pi;
				//Psi[n1] -= two_pi*round_temp;
				//}
			
		}
	}
	
	int main()
	{
		
		h = 0.01;
		h1 = 2.0 * h;
		two_pi = 2.0 * acos(-1.0);
		tn = TOTAL_TIME;
		
		x0 = 0.1;
		y01 = 1.8;
		z0 = 0.1;
		Theta0 = 0;
		Phi0 = 0;
		Psi0 = 0;
		
		L = 0;
		M = 0;
		N = 0;
		
		n = 0;
		runge(t0, x0, y01, z0, tn, h);
		n1 = 0;
		runge1(t0, Theta0, Phi0, Psi0, tn, h1);
		
		printf("\n At time %lf:\n\n", tn);
		printf("omega1 = %lf		omega2 = %lf		omega3 = %lf\n", x[n], y[n], z[n]);
		printf("theta = %lf		phi = %lf				psi = %lf\n", Theta[n1], Phi[n1], Psi[n]);
	}
