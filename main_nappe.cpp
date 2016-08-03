#include <stdio.h>
#include <stdlib.h> 
#include "Sinus.h"
#include <iostream>
#include <fstream>
#include "Gaussian.h"
#include "BLineaire.h"
#include "Constant.h"
#include "DeuxX.h"
#include <math.h>
#include <string>
#include <time.h> 
#include "XBeta.h"
#include <sstream>
#include <algorithm> 
#include "Quantification.h"
#include "Density.h"
#include "Diffusion1d.h"
#include "Cosinus.h"
#include "Xcosinus.h"

#define _USE_MATH_DEFINES
# define M_PI           3.14159265358979323846  /* pi */

using namespace std;

// function to convert an integer to a string
std::string to_string(int i)
{
    std::stringstream ss;
    ss << i;
    return ss.str();
}

int main()
{
	ofstream fparabol("Result_nappe/parabol.data", ios::out);
	ofstream fresult;
	ofstream fm_theo("Result_nappe/m_theo.data", ios::out);
	ofstream fsigma_theo("Result_nappe/sigma_theo.data", ios::out);
	ofstream fa_n("Result_nappe/a_n.data", ios::out);
	
	
	//Parameters
	double beta=1.;// phi is beta-Lipschitz
	int nb_theta=5;// number of napps according to different theta
	double theta_min=1.0/double(2.+beta); // the minimum theta of the napps
	double theta_max=1.;// the maximum theta of the napps
	int d=1;// diffusion dimension
	int r=1;// innovation diffusion
	int n=5*pow(10,4);// number of diffusion loops
	int MC=pow(10,4); // number of Monte-Carlo loops
	int M=30;// number of quantization in computation of temp_dens_Gn and temp_dens_corrected, becareful the complexity is in M^2
	int M_U=1;// number of points of the quantized random
	int n_dens=100;// number of points for the discretization of the density
	double a_inf=0.001;// the minimum of the interval of interest for the density
	double a_sup=2.;// the minimum of the interval of interest for the density
	
	// initialization
	double variance=0;
	double proba, moyenne;
	double temp;
	double m_theo, sigma_theo=0;
	double gamma_0=1;
	int i=0;
	double X0;
	double Gamma=0;
	int j;
	double a;
	double a_n=0;
	
	string file_name;
	int k=0;
	std::stringstream oss;
	double theta;
	double temp_dens;
	double gam=0;

	//diffusion functions
	BLineaire b;
		//Constant sigma;
	Cosinus sigma;
	// test function
		//XBeta phi(beta);
	Sinus phi;
	// random variables
	
	Gaussian U;
	Gaussian U_final;	
	Quantification UM(M_U,r);
	Quantification T1(M);
	Quantification T2(M);

	for(k=0;k<nb_theta+1; k++)
	{
		theta=theta_min*(double(nb_theta)-k)/(double(nb_theta))+theta_max*k/(double(nb_theta));// the current theta of the napp
		
			for (j = 1; j < n + 1; j++)
		Gamma = Gamma + pow(double(j), -theta);
		
		file_name="Result_nappe/k"+to_string(k);
		fresult.open(file_name.c_str());
	
		Density dens(n_dens, Gamma);
		
		// estimation of \tilde \gamma
		gam=0;
		for(i=1;i<n+1;i++)
			gam+=pow(i,-(3.0+beta)*0.5*theta);
		gam/=sqrt(Gamma);
		
		// computing of theoretical moments 
		dens.set_moments(MC, M_U ,Gamma, beta, U_final, phi,b ,sigma, UM);
		m_theo=dens.get_m();
		m_theo*=gam;
		sigma_theo=dens.get_sigma();
		fsigma_theo<<sigma_theo;
		fm_theo<<m_theo;
		fm_theo.close();
		fsigma_theo.close();
		
		Diffusion1d X(theta, n, b, sigma, U, M_U, UM,  T1, T2, gamma_0);
		proba=0;
		variance=0;
	
		for(i=0; i< MC; i++)
		{
			X0=U.get_value();
			temp_dens=X.Run_measure(X0, phi);
			proba+=temp_dens;
			variance+=temp_dens*temp_dens;
			dens.set_gn(temp_dens,a_inf,a_sup);
		}
		proba/=MC;
		variance/=MC;
		variance-=moyenne*moyenne;
		cout<<"moyenne="<<moyenne<<endl;
		cout<<"variance="<<variance<<endl;
		
		// computing of a_n
		a_n=dens.get_an();
		a_n*=gam*1/(4.*3.*2.);// Be careful to multiply by \sup \sigma * \sup \nabla \phi, here it is 1
		fa_n<<a_n<<"\t"<<0;
		
		// g_{n,\theta} computing
		for(i=0; i< n_dens; i++)
		{
			a=a_inf+(a_sup-a_inf)*double(i)/double(n_dens);
			fresult<<a<<"\t"<<log(dens.get_gn(i)/double(MC))<<endl;
		}

			fresult.close();
}
		// Sn computing
		for(i=0; i< n_dens; i++)
		{
			a=a_inf+(a_sup-a_inf)*double(i)/double(n_dens);
			fparabol<<a<<"\t"<<-(a-a_n)*(a-a_n)/(2.*(1.)*1.0)<<endl;
		}
		fa_n.close();
		return 0;
}
