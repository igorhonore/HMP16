#include <stdio.h>
#include <stdlib.h> 
#include "Sinus.h"
#include <iostream>
#include <fstream>
#include "BLineaire.h"
#include "Constant.h"
#include "DeuxX.h"
#include <math.h>
#include <string>
#include <time.h> 
#include "XBeta.h"
#include <sstream>
#include "Cosinus.h"
#include "Xcosinus.h"
#include "Cosinus2.h"

#define _USE_MATH_DEFINES
# define M_PI    3.14159265358979323846  /* pi */

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
	ifstream fmeasure("Result_nappe/measure.data");
	
	//Parameters
	double beta=1.;// phi is beta-Lipschitz
	int nb_theta=5;// number of napps according to different theta
	double theta_min=1.0/double(2.+beta); // the minimum theta of the napps
	double theta_max=1.;// the maximum theta of the napps
	int n=pow(10,4);// number of diffusion loops
	int n_dens=100;// number of points for the discretization of the density
	double a_inf=0.001;// the minimum of the interval of interest for the density
	double a_sup=2.;// the minimum of the interval of interest for the density
	int n_rho=5000;// number of iterations for rho optimization
	double grad_phi=1.;//+0.01; // supremum of \nabla \varphi
	double grad_phi_sig=1.;//+0.01; // supremum of \nabla \vartheta, the value is often unknow so we choose a constant with the good order. Anyway, the result does not seem to be very sensitive to it.
	double sigma_inf=1.;// supremum of \sigma
	double rho_max=1.75;// upper bound of the rho maximization interval
	double rho_min=1.+pow(10.,-15.);// lower bound of the rho maximization interval
	
	// initialization
	double measure;
	int i,j,k=0;
	double Gamma=0;
	double a;
	string file_name;
	std::stringstream oss;
	double theta;
	double rho;
	double Bn;
	double An;
	double q, bar_q,p;
	double phi_n;
	double C;
	double P_lambda_inf;
	double P_lambda;
	fmeasure>>measure;
	
	for(k=0;k<nb_theta+1; k++)
	{
		theta=theta_min*(double(nb_theta)-k)/(double(nb_theta))+theta_max*k/(double(nb_theta));// the current theta of the napp
		for (j = 1; j < n + 1; j++)
			Gamma = Gamma + pow(double(j), -theta);

	file_name="Result_nappe/P_lambda_min_k"+to_string(k);
	fresult.open(file_name.c_str());
	q=1.;// you can change the value to be more accurate, a term going to zero with n seems to be apropriate
	bar_q=q;// same remark here
	p=bar_q/(bar_q-1.);
		for(i=0; i<=n_dens; i++)
		{
			a=a_inf+(a_sup-a_inf)*double(i)/double(n_dens);
			fresult<<a;
			P_lambda_inf=0.;
			for(j=0;j<=n_rho;j++)
			{
				rho=rho_min+j/double(n_rho)*(rho_max-rho_min);
				Bn=pow(grad_phi,4)*(sigma_inf*sigma_inf*grad_phi_sig*grad_phi_sig/2.)/4.;// possibly adding a term in p*\tild gamma^2
				An=grad_phi*grad_phi*measure/2.;
				C=a*a/(Bn*Bn*Gamma)+(rho-1.)*pow(2.*An/(3.*Bn),3);
				phi_n=pow(sqrt(C)+a/(sqrt(Gamma)*Bn),1./3.)-pow(sqrt(C)-a/(sqrt(Gamma)*Bn),1./3.);
				P_lambda=-phi_n*pow(rho-1.,1./3.)*sqrt(Gamma)*(3.*a+An*sqrt(Gamma)*phi_n*pow(rho-1.,1./3.))/(8.*rho);// lambda opimization of the polynomial
				if(P_lambda<P_lambda_inf)
					P_lambda_inf=P_lambda;// rho optimization
			}
			fresult<<"\t"<<P_lambda_inf;
			fresult<<endl;
		}
	}

	return 0;
}
