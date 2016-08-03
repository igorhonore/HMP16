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

int main()
{
	ofstream fmeasure("Result_nappe/measure.data", ios::out);
	
	//Parameters
	double beta=0.5;// phi is beta-Lipschitz
	double theta=1.0/double(2.0+beta);// gamma_n= n^{-theta}
	int d=1;// diffusion dimension
	int r=1;// innovation diffusion
	int n=5*pow(10,4);// number of diffusion loops
	int MC=pow(10,4); // number of Monte-Carlo loops
	int M=30;// number of quantization in computation of temp_dens_Gn and temp_dens_corrected, becareful the complexity is in M^2
	int M_U=1;// number of points of the quantized random
	
	// initialization
	double temp,temp_dens;
	double gamma_0=1;
	int i=0;
	double X0;
	double variance=0;
	
	//diffusion functions
	BLineaire b;
		//Constant sigma;
	Cosinus sigma;
	// test function
		//XBeta phi(beta);
	Sinus phi;
	// random variables
	
	Gaussian U_0;
	Gaussian U;
	Quantification UM(M_U,r);
	Quantification T1(M);
	Quantification T2(M);
	
	Diffusion1d X(theta, n, b, sigma, U, M_U, UM,  T1, T2, gamma_0);
	temp=0;
	for(i=0; i< MC; i++)
	{
		X0=U_0.get_value();
		temp=X.Measure(X0, phi); // nu(phi)
		temp_dens+=temp;
		variance+=temp*temp;
	}
	temp_dens/=MC;
	variance/=MC;
	variance-=temp_dens*temp_dens;
	cout<<"mean="<<temp_dens<<endl;
	cout<<"variance="<<variance<<endl;
	fmeasure<<temp_dens<<endl;
	fmeasure.close();
	
	return 0;
}
