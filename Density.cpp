#include <cmath>
#include <iostream>
#include <vector>
#include "Random.h"
#include "Density.h"

using namespace::std;

void Density::set_value(double x, double a_inf, double a_sup, int MC)
{
	int i=0;
	double temp_dens1, temp_dens2;
	temp_dens1=a_inf;
	for(i=0;i<m_n+1;i++)
	{
		temp_dens1=a_inf+(a_sup-a_inf)*i/double(m_n);
		temp_dens2=a_inf+(a_sup-a_inf)*(i+1.0)/double(m_n);
		if((x>=temp_dens1)&&(x<=temp_dens2))
			m_vector[i]+=1.;	
	}
}

double Density::operator [] (const int & i) const
{
	return m_vector[i];
}

double Density::get_m()
{
	return m_theo;
}

double Density::get_sigma()
{
	return sigma_theo;
}

double Density::get_an()
{
	return a_n;
}


void Density::set_moments(int MC, int MC_bis, double Gamma, double beta, Random & U_final, Fonction & phi, Fonction & b, Fonction & sigma, Quantification & UM)
{
	int i, j;
	double X0;
	double psi;
	double temp_dens, U_int;
	sigma_theo=0;
	m_theo=0;
	a_n=0;
	for(j=0;j<MC_bis+1;j++)
		{
			a_n+=pow(fabs(UM.get_quantif(j,1)),3.0+beta)*UM.get_probability()[j];
		}
	for(i=0; i< MC; i++)
	{

		U_int=sigma.get_value(0.)*U_final.get_value()/(sqrt(-2.0*b.get_value(1.0)));
		sigma_theo+=sigma.get_value(U_int)*sigma.get_value(U_int)*phi.grad(U_int)*phi.grad(U_int);
		m_theo-=0.5*b.get_value(U_int)*b.get_value(U_int)*phi.laplacian(U_int);
		psi=0;
		temp_dens=0;
		for(j=0;j<MC_bis+1;j++)
		{
			temp_dens=0.5*phi.D3(U_int)*b.get_value(U_int)*sigma.get_value(U_int)*sigma.get_value(U_int)*UM.get_quantif(j,1)*UM.get_quantif(j,1)
			+1.0/(double(24))*phi.D4(U_int)*pow(sigma.get_value(U_int)*UM.get_quantif(j,1),4);
			temp_dens*=UM.get_probability()[j];
			psi+=temp_dens;
		}
		m_theo-=psi;
	}
	a_n/double(MC_bis+1.);
	m_theo/=double(MC);
	sigma_theo/=double(MC);
}

void Density::set_gn(double x, double a_inf, double a_sup)
{
	int j;
	double a;
	for(j=0; j<m_n+1; j++)
	{
		a=a_inf+(a_sup-a_inf)*double(j)/double(m_n);
		gn[j]+=(fabs(std::sqrt(m_Gamma)*x)>=a);
	}	
}

double Density::get_gn(int i)
{
	return gn[i];
}
