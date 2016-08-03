#ifndef Diffusion1d1D_H
#define Diffusion1d1D_H

#include <cmath>
#include <iostream>
#include <vector>
#include <iterator>
#include "Fonction.h"
#include "Random.h"
#include "Quantification.h"

class Diffusion1d
{
	protected:

		Fonction& m_sigma; //stochastic coefficient
		Fonction& m_b; //drift
		Random &m_U; //innovation
		double m_Gamma;
		double m_theta;
		int m_n; //number of diffusion 
		int m_M; // quantification size
		double  m_gamma0; // first time step, it may be optimised in the main
		Quantification & m_T1;// quantification for E_n^M
		Quantification & m_T2;// quantification for E_n^{1,M}
		Quantification & m_UM;// quantification to compute expectation value according to the innovation measure
	public:

		Diffusion1d (double theta, int n, Fonction &b, Fonction &sigma, Random &U, int M, Quantification & UM, Quantification & T1, Quantification & T2, double gamma_0): 
			m_theta(theta), m_n(n), m_b(b), m_sigma(sigma),m_U(U),  m_M(M), m_UM(UM), m_T1(T1), m_T2(T2), m_gamma0(gamma_0)
		{
			int i;
			m_Gamma=0;
			for (i = 1; i < n + 1; i++)
				m_Gamma += m_gamma0*pow(double(i), -m_theta);
		};
		double get_value(double x, int k);
		double Aphi( double x, Fonction &phi, int k); // return A phi ( X_k)
		double Run_measure( double x, Fonction &phi); // return nu_n(A phi)
		double Run_measure_Gn( double x, Fonction &phi); // return nu_n(A phi) +\bar G_n/Gamma_n
		double G_n(double x, int k, Fonction &phi); // return \bar G_n
		double Theta_k(double x, int k, double T1, double T2, Fonction & phi); // return \Theta_k
		double D_sigma(double x, int k, Fonction &phi); // return D_sigma 
		double D_b(double x, int k, Fonction &phi); // return D_b
		double Run_measure_corrected( double x, Fonction &phi);  // return nu_n(A phi) +(\bar G_n + D_{2,n})/Gamma_n
		double Measure( double x, Fonction &phi);// return nu_n(phi)
};

#endif
