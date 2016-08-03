#include "XBeta.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "Fonction.h"
#include <stdlib.h>

using namespace::std;

	double XBeta::get_value(double x)
{
	return pow(fabs(x),3.0+m_beta)/(1+pow(fabs(x),2.0+m_beta));
}

	double XBeta::grad(double x)
	{
		double result;
		if(x==0)
			result=0;
		else
			result= (3.0+m_beta)*x*pow(fabs(x),1.0+m_beta)/(1+pow(fabs(x),2.0+m_beta))-(2.0+m_beta)*x*pow(fabs(x),3.0+2.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2);
		return result;
	}

	double XBeta::laplacian(double x)
	{
		double result;
		if(x==0)
			result=0;
		else
		result= (3.0+m_beta)*(1.0+m_beta)*pow(fabs(x),1.0+m_beta)/(1+pow(fabs(x),2.0+m_beta))
		+(3.0+m_beta)*pow(fabs(x),1.0+m_beta)/(1+pow(fabs(x),2.0+m_beta))
		-(3.0+m_beta)*(2.0+m_beta)*pow(fabs(x),3.0+2.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2)
		-(2.0+m_beta)*pow(fabs(x),3.0+2.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2)
		-(2.0+m_beta)*(3.0+2.0*m_beta)*x*pow(fabs(x),1.0+2.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2)
		+2.0*pow(2.0+m_beta,2.0)*x*pow(fabs(x),3.0+3.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2);
		return result;
	}
		double XBeta::D3(double x)
	{
		double result;
		if(x==0)
			result=0;
		else
		result= (3.0+m_beta)*pow(1.0+m_beta,2)*x*pow(fabs(x),-1.0+m_beta)/(1+pow(fabs(x),2.0+m_beta))
		-(3.0+m_beta)*(fabs(x),2.0+m_beta)*(1.0+m_beta)*pow(fabs(x),1.0+2.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2)
		+(3.0+m_beta)*(1.0+m_beta)*x*pow(fabs(x),-1.0+m_beta)/(1+pow(fabs(x),2.0+m_beta))
		+(3.0+m_beta)*(2.0+m_beta)*x*pow(fabs(x),1.0+2.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2)
		-(3.0+m_beta)*(2.0+m_beta)*(3.0+2.0*m_beta)*x*pow(fabs(x),2.0+2.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2)
		-(3.0+m_beta)*pow(2.0+m_beta,2)*x*pow(fabs(x),3.0+3.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2)
		-(2.0+m_beta)*(3.0+2.0*m_beta)*x*pow(fabs(x),1.0+2.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2)
		+2.0*pow(2.0+m_beta,2)*x*pow(fabs(x),3.0+3.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),3)
		-(2.0+m_beta)*(3.0+2.0*m_beta)*pow(fabs(x),1.0+2.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2)
		-(2.0+m_beta)*(3.0+2.0*m_beta)*(1.0+2.0*m_beta)*pow(fabs(x),1.0+2.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2)
		+2.0*pow(2.0+m_beta, 2)*(3.0+2.0*m_beta)*pow(fabs(x),3.0+3.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),3)
		+2.0*pow(2.0+m_beta,2.0)*pow(fabs(x),3.0+3.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2)
		+2.0*pow(2.0+m_beta,2.0)*(3.0+3.0*m_beta)*pow(fabs(x),3.0+3.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),2)
		-6.0*pow(2.0+m_beta,3.0)*pow(fabs(x),5.0+4.0*m_beta)/pow(1+pow(fabs(x),2.0+m_beta),3);
		return result;
	}
	
double*** XBeta::D3(double* x)
	{
		int i,j,k,d;
		double temp;
		d=sizeof(x)/sizeof(double);
		for (i = 0; i < d; i++)
			temp=fabs(x[i]);
		double *** result;
		result= (double***)malloc(d*sizeof(double**));
  for (i = 0; i < d; i++)
  {
	  result[i] = (double**) malloc(d * sizeof(double*));
	  for(j=0;j<d;j++)
		result[i][j] = (double*) calloc(d , sizeof(double));
} 
for(i=0;i<d;i++)
{
	for(j=0; j<d; j++)
	{
		for(k=0;k<d;k++)
		{
		if(x[i]*x[j]*x[k]==0)
			result[i][j][k]=1;
		else
		result[i][j][k]=(pow((m_beta+3.0),3.0)*pow(temp,m_beta+3.0))/(x[i]*x[j]*x[k]*(pow(temp,m_beta+2)+1))
		-(3*pow(m_beta,2)*pow(temp,m_beta+3.0))/(x[i]*x[j]*x[k]*(pow(temp,m_beta+2)+1))
		+(2*(m_beta+3)*pow(temp,m_beta+3.0))/(x[i]*x[j]*x[k]*(pow(temp,m_beta+2)+1))
		-((m_beta+2)*pow((2*m_beta+5),2)*pow(temp,(2*m_beta+5)))/(x[i]*x[j]*x[k]*pow(pow(temp,m_beta+2)+1,2))
		-((m_beta+2)*pow(m_beta,2)*pow(temp,(2*m_beta+5)))/(x[i]*x[j]*x[k]*pow(pow(temp,m_beta+2)+1,2))
		-((m_beta+2)*(m_beta+3)*(2*m_beta+5)*pow(temp,(2*m_beta+5)))/(x[i]*x[j]*x[k]*pow(pow(temp,m_beta+2)+1,2))
		+(3*(m_beta+2)*(2*m_beta+5)*pow(temp,(2*m_beta+5)))/(x[i]*x[j]*x[k]*pow(pow(temp,m_beta+2)+1,2))
		+(3*(m_beta+2)*(m_beta+3)*pow(temp,(2*m_beta+5)))/(x[i]*x[j]*x[k]*pow(pow(temp,m_beta+2)+1,2))
		-(2*(m_beta+2)*pow(temp,(2*m_beta+5)))/(x[i]*x[j]*x[k]*pow(pow(temp,m_beta+2)+1,2))
		- 6*pow(m_beta,2)*pow(temp,3*m_beta+7)/(x[i]*x[j]*x[k]*pow(pow(temp,m_beta+2)+1,3))
		+2*pow(m_beta,2)*(3*m_beta+7)*pow(temp,3*m_beta+7)/(x[i]*x[j]*x[k]*pow(pow(temp,m_beta+2)+1,3))
		+2*pow(m_beta,2)*(2*m_beta+5)*pow(temp,3*m_beta+7)/(x[i]*x[j]*x[k]*pow(pow(temp,m_beta+2)+1,3))
		+2*pow(m_beta,2)*(m_beta+3)*pow(temp,3*m_beta+7)/(x[i]*x[j]*x[k]*pow(pow(temp,m_beta+2)+1,3))
		-(6*pow(m_beta+2,3)*pow(temp,4*m_beta+9))/(x[i]*x[j]*x[k]*pow(pow(temp,m_beta+2)+1,4));
		}
	}
}
		return result;
		for (i = 0; i < d; i++)
	{
		  for(j=0;j<d;j++)
			free(result[i][j]);
		free(result[i]);
	} 
	free(result);
	}

		double XBeta::D4(double x)
	{
		double result;
		if(x==0)
			result=0;
		else
		result = (pow((m_beta+3),4)*pow(fabs(x),m_beta+3.0))/(pow(x,4)*pow(fabs(x),m_beta+2.0)+1)
		-(6*pow((m_beta+3.0),3.0)*pow(fabs(x),m_beta+3.0))/(pow(x,4)*pow(fabs(x),m_beta+2.0)+1)
		+(11*pow(m_beta,2)*pow(fabs(x),m_beta+3.0))/(pow(x,4)*pow(fabs(x),m_beta+2.0)+1)
		-(6*(m_beta+3)*pow(fabs(x),m_beta+3.0))/(pow(x,4)*pow(fabs(x),m_beta+2.0)+1)
		-((m_beta+2)*pow(2*m_beta+5,3)*pow(fabs(x),2*m_beta+5))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,2))
		-((m_beta+2)*(m_beta+3)*pow(2*m_beta+5.0,2)*pow(fabs(x),2*m_beta+5))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,2))
		+(6*(m_beta+2)*pow(2*m_beta+5.0,2)*pow(fabs(x),2*m_beta+5))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,2))
		-((m_beta+2)*pow((m_beta+3.0),3.0)*pow(fabs(x),2*m_beta+5))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,2))
		+(6*(m_beta+2)*pow(m_beta,2)*pow(fabs(x),2*m_beta+5))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,2))
		-((m_beta+2)*pow(m_beta,2)*(2*m_beta+5)*pow(fabs(x),2*m_beta+5))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,2))
		+(6*(m_beta+2)*(m_beta+3)*(2*m_beta+5)*pow(fabs(x),2*m_beta+5))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,2))
		-(11*(m_beta+2)*(2*m_beta+5)*pow(fabs(x),2*m_beta+5))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,2))
		-(11*(m_beta+2)*(m_beta+3)*pow(fabs(x),2*m_beta+5))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,2))
		+(6*(m_beta+2)*pow(fabs(x),2*m_beta+5))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,2))
		+(2*pow(m_beta,2)*pow(3*m_beta+7.0,2)*pow(fabs(x),3*m_beta+7))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,3))
		+(2*pow(m_beta,2)*pow(2*m_beta+5.0,2)*pow(fabs(x),3*m_beta+7))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,3))
		+(2*pow(m_beta,2)*pow(m_beta,2)*pow(fabs(x),3*m_beta+7))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,3))
		+(22*pow(m_beta,2)*pow(fabs(x),3*m_beta+7))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,3))
		-(12*pow(m_beta,2)*(3*m_beta+7)*pow(fabs(x),3*m_beta+7))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,3))
		+(2*pow(m_beta,2)*(2*m_beta+5)*(3*m_beta+7)*pow(fabs(x),3*m_beta+7))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,3))
		+(2*pow(m_beta,2)*(m_beta+3)*(3*m_beta+7)*pow(fabs(x),3*m_beta+7))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,3))
		-(12*pow(m_beta,2)*(2*m_beta+5)*pow(fabs(x),3*m_beta+7))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,3))
		+(2*pow(m_beta,2)*(m_beta+3)*(2*m_beta+5)*pow(fabs(x),3*m_beta+7))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,3))
		-(12*pow(m_beta,2)*(m_beta+3)*pow(fabs(x),3*m_beta+7))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,3))
		+(36*pow(m_beta+2.0,3)*pow(fabs(x),4*m_beta+9))/(pow(x,4)*(pow(pow(fabs(x),m_beta+2.0)+1,4)))
		-(6*pow(m_beta+2,3)*(4*m_beta+9)*pow(fabs(x),4*m_beta+9))/(pow(x,4)*(pow(pow(fabs(x),m_beta+2.0)+1,4)))
		-(6*pow(m_beta+2,3)*(3*m_beta+7)*pow(fabs(x),4*m_beta+9))/(pow(x,4)*(pow(pow(fabs(x),m_beta+2.0)+1,4)))
		-(6*pow(m_beta+2,3)*(2*m_beta+5)*pow(fabs(x),4*m_beta+9))/(pow(x,4)*(pow(pow(fabs(x),m_beta+2.0)+1,4)))
		-(6*pow(m_beta+2,3)*(m_beta+3)*pow(fabs(x),4*m_beta+9))/(pow(x,4)*(pow(pow(fabs(x),m_beta+2.0)+1,4)))
		+(24*pow(m_beta,2)*pow(fabs(x),5*m_beta+11))/(pow(x,4)*pow(pow(fabs(x),m_beta+2.0)+1,5));
		return result;
	}
	
	 double XBeta::get_value(std::vector<double> x)
	 {
		 int i,n;
		 n=sizeof(x);
		 double temp=0;
		for(i=1;i<n+1;i++)
			temp+=fabs(x[i]);
		temp=pow(temp,3.0+m_beta)/(1+pow(temp,2.0+m_beta));
		 return temp;
	 }
		 std::vector<double> XBeta::grad(std::vector<double> x)
	{
		 int i,n;
		 n=sizeof(x);
		 std::vector<double> result(n);
		 double temp=0;
		for(i=0;i<n;i++)
			temp+=fabs(x[i]);
		for(i=0;i<n;i++)
		{
			result[i]= (3.0+m_beta)* pow(fabs(temp),1.0+m_beta)/(1+pow(temp,2.0+m_beta)) - (2.0+m_beta)* pow(fabs(temp),3.0+2.0*m_beta)/pow(1+pow(temp,2.0+m_beta),2);
			result[i]*=x[i];
		}
		 return result;
	 }
std::vector<std::vector<double> > XBeta::laplacian(std::vector<double> x)
	{
	
		 int i,j,n;
		 n=sizeof(x);
		 	 std::vector< std::vector<double> > result;
		 	 
		 	 result.resize(n);
		 for(i=0;i<n;i++)
		 	result[i].resize(n);
		 double temp=0;
		for(i=0;i<n;i++)
			temp+=fabs(x[i]);
		for(i=0;i<n;i++)
		{
			result[i][i]= (3.0+m_beta)* pow(fabs(temp),m_beta+1.0)/(1+pow(temp,2.0+m_beta)) - (2.0+m_beta)* pow(fabs(temp),3.0+2.0*m_beta)/pow(1+pow(temp,2.0+m_beta),2);
			for(j=0;j<n;j++)
			{
				result[i][j]+= (3.0+m_beta)* (1.0+m_beta)*pow(fabs(temp),m_beta-1.0)
				- (3.0+m_beta)*(2.0+m_beta)* pow(fabs(temp),1.0+2.0*m_beta)
				- (2.0+m_beta)*(3.0+2.0*m_beta)* pow(fabs(temp),1.0+2.0*m_beta)
				+ 2.0*(2.0+m_beta)*(2.0+m_beta)*x[i]*x[j]* pow(fabs(temp),3.0+3.0*m_beta)/pow(1+pow(temp,2.0+m_beta),2);
				result[i][j]*=x[i]*x[j]/pow(1+pow(temp,2.0+m_beta),2);
			}
		}
		 return result;
	 }
		 std::vector<std::vector<double> > XBeta::get_value( std::vector<double>  x, int d, int r)
	 {
		 std::vector<std::vector<double> > result;
		 return result;
	 }
	 
std::vector<double > XBeta::get_value( std::vector<double>  x, int d)
	 {
		 std::vector<double> result;
		 return result;
	 }
	 
double * XBeta::get_value( double *  x, int d)
	 	{
			int i;
			double * result;
			return result;
		}
double ** XBeta::get_value( double * x, int r,int d)
{
			int i;
			double ** result;
			return result;
}

double *XBeta::grad(double * x)
{
	double * result;
	int i,n;
	n=sizeof(x)/sizeof(double);
	double temp=0;
	result=(double*)malloc(n*sizeof(double));
	for(i=0;i<n;i++)
		temp+=fabs(x[i]);
	for(i=0;i<n;i++)
	{
		if(x[i]==0)
			result[i]=0;
		else
			result[i]= (3.0+m_beta)*pow(fabs(temp),3.0+m_beta)/(x[i]*(1+pow(temp,2.0+m_beta)))- (2.0+m_beta)* pow(fabs(temp),5.0+2.0*m_beta)/(x[i]*pow(1+pow(temp,2.0+m_beta),2));
	}
	return result;
}
double ** XBeta::laplacian(double * x)
{
	int i,j,n;
	n=sizeof(x)/sizeof(double);
	double temp=0;
	for(i=0;i<n;i++)
		temp+=fabs(x[i]);
	double ** result;
	result=(double **) malloc(n*sizeof(double*));
	for (i = 0; i < n; i++)
	{
		result[i]=(double*)calloc(n,sizeof(double));
	} 
	for (i = 0; i < n; i++)
	{
		for(j=0;j<n;j++)
		{
			if(x[i]*x[j]==0)
				result[i][j]=0;
			else
				result[i][j]=(pow(m_beta,2)*pow(temp,m_beta+3.0))/(x[i]*x[j]*pow(temp,m_beta+2.0)+1)
				-((m_beta+3)*pow(temp,m_beta+3.0))/(x[i]*x[j]*pow(temp,m_beta+2.0)+1)
				-((m_beta+2)*(2*m_beta+5)*pow(temp,(2*m_beta+5)))/(x[i]*x[j]*pow(pow(temp,m_beta+2.0)+1,2.0))
				-((m_beta+2)*(m_beta+3)*pow(temp,(2*m_beta+5)))/(x[i]*x[j]*pow(pow(temp,m_beta+2.0)+1,2.0))
				+((m_beta+2)*pow(temp,(2*m_beta+5)))/(x[i]*x[j]*pow(pow(temp,m_beta+2.0)+1,2.0))
				+(2*pow(m_beta,2)*pow(temp,3*m_beta+7)/(x[i]*x[j]*pow(pow(temp,m_beta+2.0)+1,3)));
		}
	} 
	return result;
	for (i = 0; i < n; i++)
	{
		free(result[i]);
	} 
	free(result);
}
		
		double XBeta::get_value( double *  x)
		{
			return pow(fabs(x[0]),3.0+m_beta)/(1+pow(fabs(x[0]),2.0+m_beta));
			
		}
