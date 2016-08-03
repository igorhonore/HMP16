#include <cmath>
#include <iostream>
#include <vector>
#include "Diffusion1d.h"
#include "Random.h"
#include <stdlib.h>

using namespace::std;

double Diffusion1d::Measure( double x, Fonction &phi)
{
	double result=0;
	double x_n;
	x_n=x;
	double gamma;
	int i=0;
	gamma =m_gamma0* pow(double(i+1.0),-m_theta);
	result+=gamma*phi.get_value(x_n);
	for(i=1;i<m_n;i++)
	{
		x_n=get_value(x_n,i-1);
		gamma =m_gamma0* pow(double(i+1.),-m_theta);
		result+=gamma*phi.get_value(x_n);
	}
	result/=m_Gamma;
	return result;
}

double Diffusion1d::Run_measure( double x, Fonction &phi)
{
	double result=0;
	double x_n;
	x_n=x;
	double gamma;
	int i=0;
	gamma =m_gamma0* pow(double(i+1.0),-m_theta);
	result+=gamma*Aphi(x_n, phi,0);
	for(i=1;i<m_n;i++)
	{
		x_n=get_value(x_n,i-1);
		gamma =m_gamma0* pow(double(i+1.),-m_theta);
		result+=gamma*Aphi(x_n, phi,i);
	}
	result/=m_Gamma;
	return result;
}

double Diffusion1d::Run_measure_Gn( double x, Fonction &phi)
{
	double result=0;
	double x_n;
	x_n=x;
	double gamma;
	double test;
	int i=0;
	gamma = m_gamma0*pow(double(i+1.), -m_theta);
	result=gamma*Aphi(x_n, phi,0);
	result+=G_n(x_n,i , phi);
	for(i=1;i<m_n;i++)
	{
		x_n=get_value(x_n,i-1);
		gamma = m_gamma0*pow(double(i+1.0), -m_theta);
		result+=gamma*Aphi(x_n, phi,i);
		result+=G_n(x_n,i, phi);
	}
	result/=m_Gamma;
	return result;
}


double Diffusion1d::Run_measure_corrected( double x, Fonction &phi)
{
	double result=0;
	double x_n;
	x_n=x;
	double gamma;
	double test=0;
	int i=0;
	gamma = m_gamma0*pow(double(i+1.0), -m_theta);
	result+=gamma*Aphi(x_n, phi,i);
	result+=G_n(x_n,i , phi);
	result+=D_sigma(x_n,i, phi)+D_b(x_n,i, phi);
	for(i=1;i<m_n;i++)
	{
		x_n=get_value(x_n,i-1);
		gamma = m_gamma0*pow(double(i+1.0), -m_theta);
		result+=gamma*Aphi(x_n, phi,i);
		result+=G_n(x_n,i, phi);
		result+=D_sigma(x_n,i, phi)+D_b(x_n,i, phi);
	}
	result/=m_Gamma;
	return result;
}

double Diffusion1d::Aphi( double x, Fonction &phi, int k)
{
	double result=0;
	double temp_sigma;
	temp_sigma=m_sigma.get_value(x);
	result=m_b.get_value(x) * phi.grad(x)+ 0.5*temp_sigma * temp_sigma*phi.laplacian(x);
	return result;
}

double Diffusion1d::get_value(double x, int k)
{
	double gamma;
	gamma = m_gamma0*pow(double(k+1.0), -m_theta);
	double result;
	result=x+ gamma*m_b.get_value(x)+sqrt(gamma)*m_sigma.get_value(x)*m_U.get_value();
	return result;
}

double Diffusion1d::G_n(double x, int k, Fonction &phi)
{
	int j,l, M1, M2;
	double result, gamma,temp;
	M1=m_T1.get_size();
	M2=m_T2.get_size();
	result=0.;
	gamma=m_gamma0*pow(k+1.,-1.5*m_theta);
	temp=0;
	for(j=1; j<M1+1; j++)
	{
		for(l=1;l<M2+1;l++)
		{
				temp+=m_T1.Uniforme(j)*(1.0-m_T1.Uniforme(j))*Theta_k(x, k, m_T1.Uniforme(j),m_T2.Uniforme(l), phi);
		}
	}
	result=gamma*temp;
	result/=double(M1*M2);
	return result;
}

double Diffusion1d::Theta_k(double x, int k, double T1, double T2, Fonction &phi)
{
	int i,j,l,r,m, M_U;
	double result=0;
	double result_final=0;
	double gamma=0;
	gamma=m_gamma0*pow(k+1., -m_theta);
	double temp;
	double temp1;
	double sigmaU;
	double X_k;
	double temp2;
	double P_U;
	double U_M;

	M_U=m_UM.get_size();
	
	double temp_b;
	double temp_sigma;
	double temp_D3;
	temp_b=m_b.get_value(x);
	temp_sigma=m_sigma.get_value(x);
	
	for(m=0;m<M_U+1;m++)
	{
		result=0;
		X_k=x+gamma*temp_b + T1*T2*sqrt(gamma)* temp_sigma*m_UM.get_quantif(m,1);
		sigmaU= temp_sigma*m_UM.get_quantif(m,1);
		temp= phi.D3(X_k)*sigmaU*sigmaU *sigmaU;
		result_final+=temp*m_UM.get_probability(m);
	return result_final;
	}
}

double Diffusion1d::D_sigma(double x, int k, Fonction &phi)
{
	double result, gamma;
	double Xk;
	result=0;
	gamma=m_gamma0*pow(k+1.,-m_theta);
	double temp_sigma;
	temp_sigma=m_sigma.get_value(x);
	Xk=x+gamma*m_b.get_value(x);
	double temp_laplk;
	temp_laplk=phi.laplacian(Xk);
	double temp_lapl;
	temp_lapl=phi.laplacian(x);
	result=(temp_laplk-temp_lapl)*temp_sigma*temp_sigma;
	result*=0.5*gamma;
	return result;
}

double Diffusion1d::D_b(double x, int k, Fonction &phi)
{
	int i,j,l, M1;
	double result, gamma;
	M1=m_T1.get_size();
	double Xk;
	result=0;
	gamma=m_gamma0*pow(k+1.,-1.5*m_theta);
	double vect_temp;
	double temp_b;
	temp_b=m_b.get_value(x);
	double temp_lapl;
	for(l=1; l<M1+1; l++)
	{
		Xk=x+m_T1.Uniforme(l)*gamma*temp_b;
		temp_lapl=phi.laplacian(Xk);
		vect_temp=temp_lapl*temp_b*temp_b;
		result+=vect_temp*(1.-m_T1.Uniforme(l));
	}
	result*=gamma*gamma;
	result/=double(M1);
	return result;
}

