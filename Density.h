#ifndef DENSITY_H
#define DENSITY_H

#include <cmath>
#include <iostream>
#include <vector>
#include <iterator>
#include "Fonction.h"
#include "Random.h"
#include "Quantification.h"
class Density
{
	protected:
		int m_n; //number of diffusion 
		double * m_vector;// density vector
		double m_theo;// theoretical mean
		double sigma_theo;// theoretical variance, carrée du champ
		double a_n;// maximization of the bias
		std::vector<double> gn;
		double m_Gamma;// Gamma_n
	public:
		Density(double Gamma):m_Gamma(Gamma)
		{
			m_theo=0;
			sigma_theo=0;
		};

		Density(int n, double Gamma): m_n(n), m_Gamma(Gamma)
		{
			m_vector=(double*)calloc(m_n,sizeof(double));
			gn.resize(m_n+1);
			m_theo=0;
			sigma_theo=0;
			int i;
			for(i=0;i<m_n+1;i++)
				gn[i]=0;
		};
		void set_value(double x, double a_inf, double a_sup, int MC);// compute the density bewteen a_inf and a_sup
		double operator [](const int & i) const;
		double get_m();// return the mean
		double get_sigma(); // return the variance
		double get_an();// return a_n
		double get_gn(int i); // return gn(a_inf+(a_sup-a_inf)*double(i)/double(m_n))
		void set_moments(int MC, int MC_bis, double Gamma, double beta, Random & U_final, Fonction & phi, Fonction & b, Fonction & sigma, Quantification & UM); // compute a_n, the theoretical mean and the theoretical variance
		void set_gn(double x, double a_inf, double a_sup );// feed the gn vector with one realisation
		
};		

#endif
