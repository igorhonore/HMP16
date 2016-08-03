#include "Quantification.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <stdlib.h>


using namespace::std;
	
int Quantification::get_size()
{
	return m_M;
}
double  Quantification::Uniforme(int j)
{	
	double temp;
		temp=(2.*j-1.)/double(2.*m_M);
	return temp;
}
double * Quantification::Uniforme()
{	
	int i;
	double* temp;
	temp= (double*)malloc((m_M+1)*sizeof(double));
	for(i=0;i<m_M+1;i++)
	{
		temp[i]=(i)/double(m_M);
	}
	return temp;
}
double Quantification::get_probability(int i)
{	
	double temp;
		temp=1./double(m_M+1.0);
	return temp;
}
double * Quantification::get_probability()
{	
	int i;
	double* temp;
	temp= (double*)malloc((m_M+1)*sizeof(double));
	for(i=0;i<m_M+1;i++)
	{
		temp[i]=1/double(m_M+1.0);
	}
	return temp;
		free(temp);
}
double ** Quantification::get_quantif()
{
	int i,j;
	double** temp;
	temp= (double**)malloc((m_M+1)*sizeof(double*));
	for(i=0;i<m_M+1;i++)
	{
		temp[i]=(double*)calloc(m_r,sizeof(double));
	}
	for(i=0;i<m_M+1;i++)
	{
		for(j=0;j<m_r;j++)
			temp[i][j]=(-m_M+2*i)*sqrt(3.0/double(m_M*(m_M+2.0)));
	}
	return temp;
	free(temp);
}

double Quantification::get_quantif(int i, int j)
{
	double temp;
	temp=(-m_M+2.0*i)*sqrt(3.0/double(m_M*(m_M+2.0)));
	return temp;
}
