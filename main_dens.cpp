#include <stdio.h>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <time.h> 
#include <sstream>
#include <algorithm> 
#include "Density.h"
#include "Diffusion1d.h"

#define _USE_MATH_DEFINES
# define M_PI           3.14159265358979323846  /* pi */

using namespace std;

int main()
{
	ifstream fresult("Result_density/result.data");
	ifstream fresult_Gn("Result_density/result_Gn.data");
	ifstream fresult_corrected("Result_density/result_corrected.data");
	ifstream fm_theo("Result_density/m_theo.data");
	ifstream fsigma_theo("Result_density/sigma_theo.data");
	ofstream fNormal_result("Result_density/Normal_result.data", ios::out);
	
	//Parameters
	double beta=1.;// phi is beta-Lipschitz
	double theta=1.0/double(2.0+beta);// gamma_n= n^{-theta}, becareful to take the same theta as for the realizations
	int d=1;// diffusion dimension
	int r=1;// innovation diffusion
	int n=pow(10,4);// number of diffusion loops
	int MC=pow(10,4); // number of Monte-Carlo loops
	int M=10;// number of quantization in computation of temp_dens_Gn and temp_dens_corrected, becareful the complexity is in M^2
	int M_U=1;// number of points of the quantized random
	int n_dens=100;// number of points for the discretization of the density

	// initialization
	double variance=0;
	double a_inf, a_max;
	double a_inf_G, a_max_G;
	double a_inf_C, a_max_C;
	double temp;
	double m_theo, sigma_theo, gam, a_n=0;
	double gamma_0;
	int i=0;
	double X0;
	double Gamma=0;
	a_n=0;
	gam=0;
	m_theo=0;
	sigma_theo=0;
	double moy=0;
	double moy_G=0;
	double moy_C=0;
	double var=0;
	double var_G=0;
	double var_C=0;

	int j;
	for (j = 1; j < n + 1; j++)
		Gamma = Gamma + pow(double(j), -theta);
	// same mean and variance computing in main_baised.cpp	
	fm_theo>>m_theo;
	fsigma_theo>>sigma_theo;

	//optimization in gamma_0 
	gamma_0=pow(sigma_theo/(m_theo*m_theo*(1.+beta)), 1./(2.+beta));
	m_theo*=pow(gamma_0,(1+beta)*0.5);
	sigma_theo/=gamma_0;
	
	//Computing of the interval values of the empirical measure
		fresult>>a_inf;
		a_max=a_inf*sqrt(Gamma);
		moy=a_inf;
		a_inf*=sqrt(Gamma);
		while(!fresult.eof())
		{
			fresult>>temp;
			if(sqrt(Gamma)*temp<a_inf)
				a_inf=sqrt(Gamma)*temp;
			if(sqrt(Gamma)*temp>a_max)
				a_max=sqrt(Gamma)*temp;
			moy+=temp;
			var+=temp*temp;
		}
		var/=MC;
		var*=Gamma;//gamma_0;
		moy*=sqrt(Gamma)/(double(MC));//*sqrt(gamma_0));
		var-=moy*moy;
		cout<<"moyenne empirique ="<<moy<<endl;
		cout<<"variance empirique ="<<var<<endl;

	//Computing of the interval values of the empirical measure minus G_n

		fresult_Gn>>a_inf_G;
		a_max_C=sqrt(Gamma)*a_inf_G;
		moy_G=a_inf_G;
		a_inf_G*=sqrt(Gamma);
		
		while(!fresult_Gn.eof())
		{
			fresult_Gn>>temp;
			if(sqrt(Gamma)*temp<a_inf_G)
				a_inf_G=sqrt(Gamma)*temp;
			if(sqrt(Gamma)*temp>a_max_G)
				a_max_G=sqrt(Gamma)*temp;
			moy_G+=temp;
			var_G+=temp*temp;
		}
		var_G/=MC;
		var_G*=Gamma;//(gamma_0);
		moy_G*=sqrt(Gamma)/(double(MC));//*sqrt(gamma_0));
		var_G-=moy_G*moy_G;
		cout<<"moyenne empirique corrigee Gn ="<<moy_C<<endl;
		cout<<"variance empirique corrigee Gn="<<var_C<<endl;
		
	//Computing of the interval values of the empirical measure unbaised

		fresult_corrected>>a_inf_C;
		a_max_C=sqrt(Gamma)*a_inf_C;
		moy_C=a_inf_C;
		a_inf_C*=sqrt(Gamma);
		
		while(!fresult_corrected.eof())
		{
			fresult_corrected>>temp;
			if(sqrt(Gamma)*temp<a_inf_C)
				a_inf_C=sqrt(Gamma)*temp;
			if(sqrt(Gamma)*temp>a_max_C)
				a_max_C=sqrt(Gamma)*temp;
			moy_C+=temp;
			var_C+=temp*temp;
		}
		var_C/=MC;
		var_C*=Gamma;//(gamma_0);
		moy_C*=sqrt(Gamma)/(double(MC));//*sqrt(gamma_0));
		var_C-=moy_C*moy_C;
		cout<<"moyenne empirique corrigee Gn ="<<moy_C<<endl;
		cout<<"variance empirique corrigee Gn="<<var_C<<endl;


		fresult.close();
		fresult_Gn.close();
		fresult_corrected.close();
		
		fresult.open("Result_density/result.data");
		fresult_Gn.open("Result_density/result_Gn.data");
		fresult_corrected.open("Result_density/result_corrected.data");
			
		Density dens ( n_dens, Gamma);
		Density dens_Gn ( n_dens, Gamma);
		Density dens_corrected ( n_dens, Gamma);
		while(!fresult.eof())
		{
			fresult>>temp;
			dens.set_value(sqrt(Gamma)*temp,a_inf,a_max,MC);
		}
		while(!fresult_Gn.eof())
		{
			fresult_Gn>>temp;
			dens_Gn.set_value(sqrt(Gamma)*temp,a_inf,a_max,MC);
		}
		while(!fresult_corrected.eof())
		{
			fresult_corrected>>temp;
			dens_corrected.set_value(sqrt(Gamma)*temp,a_inf,a_max,MC);
		}
		fresult.close();
		fresult_Gn.close();
		fresult_corrected.close();
				
		ofstream fdens("Result_density/dens.data");
		ofstream fdens_Gn("Result_density/dens_Gn.data");
		ofstream fdens_corrected("Result_density/dens_corrected.data");
		
		for(i=0;i<n_dens;i++)
		{
			temp=a_inf+(a_max-a_inf)*i/double(n_dens);
			fdens<<a_inf+(a_max-a_inf)*i/double(n_dens)<<"\t"<<dens.operator[](i)*n_dens/(double(MC*(a_max-a_inf)))<<endl;
			fNormal_result<<exp(-(temp-m_theo)*(temp-m_theo)/(2.0*sigma_theo))/(std::sqrt(2.*M_PI*sigma_theo))<<"\t"<<a_inf+(a_max-a_inf)*i/double(n_dens)<<endl;
			
			temp=a_inf_G+(a_max_G-a_inf_G)*i/double(n_dens);
			fdens_corrected<<a_inf_G+(a_max_G-a_inf_G)*i/double(n_dens)<<"\t"<<dens_corrected[i]*n_dens/(double(MC*(a_max_G-a_inf_G)))<<endl;				
			
			temp=a_inf_C+(a_max_C-a_inf_C)*i/double(n_dens);
			fdens_corrected<<a_inf_C+(a_max_C-a_inf_C)*i/double(n_dens)<<"\t"<<dens_corrected[i]*n_dens/(double(MC*(a_max_C-a_inf_C)))<<endl;	
		}
		cout<<"moyenne theorique ="<<m_theo<<endl;
		cout<<"variance theorique ="<<sigma_theo<<endl;

		fdens.close();
		fdens_Gn.close();
		fdens_corrected.close();
		fm_theo.close();
		fsigma_theo.close();
		
		
	return 0;
}

