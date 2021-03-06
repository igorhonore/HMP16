#include "Cosinus2.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <stdlib.h>

using namespace::std;

	double ** Cosinus2::get_value( double*  x, int d, int r)
{	
	int i,j, m;
	m=min(d,r);
	d=sizeof(x);
	double** temp;
	temp= (double**)malloc(d*sizeof(double*));
	for (i = 1; i < m+1; i++)
	{
		temp[i] = (double*) calloc(d , sizeof(double));
	} 
	for(i=1;i<d+1;i++)
	{
		temp[i][i]=cos(x[i]);
	}
	return temp;
}

	double Cosinus2::get_value(double x)
{
	return cos(x)*cos(x);
}

	double Cosinus2::grad(double x)
	{
		return -2*cos(x)*sin(x);
	}

	double Cosinus2::laplacian(double x)
	{
		return -2*cos(2*x);
	}

	 double Cosinus2::get_value(std::vector<double> x)
	 {
		 double result, temp=0;
		 int i;
		 for(i=1;i<sizeof(x)+1;i++)
			temp+=x[i];
		result=cos(temp);
		 return result;
	 }
		 std::vector<double> Cosinus2::grad(std::vector<double> x)
	{
		 std::vector<double> result;
		 int i;
		 for(i=1;i<sizeof(x)+1;i++)
			result.push_back(cos(x[i]));
		 return result;
	 }
		 std::vector<std::vector<double> > Cosinus2::laplacian(std::vector<double> x)
	{
		 std::vector< std::vector<double> > result;
		 int n,i,j;
		 n=sizeof(x);
			result.resize(n);
		 for(i=1;i<n+1;i++)
		 	result[i].resize(n);
		 for(i=1;i<n+1;i++)
		 {
			 for(j=1;j<n+1;j++)
			 {
				result[i][j]=0;
			}
			result[i][i]=-cos(x[i]);
		 }
		 return result;
	 }
		 std::vector<std::vector<double> > Cosinus2::get_value( std::vector<double>  x, int d, int r)
		 {
		std::vector<std::vector<double> > result;
		 int i,j;
		 result.resize(d);
		 for(i=1;i<d+1;i++)
		 	result[i].resize(r);
		 for(i=1;i<d+1;i++)
		 {
			 for(j=1;j<r+1;j++)
				{
					result[i][j]=1;
				}
			result[i][i]=cos(x[i]);
		}
		 return result;
		 }
		 
	std::vector<double > Cosinus2::get_value( std::vector<double>  x, int d)
		 {
			 std::vector<double> result;
			 int i;
			 for(i=1;i<d+1;i++)
				result.push_back(cos(x[i]));
			 return result;
		 }
		 
	double Cosinus2::D3(double x)
	{
		return 4*sin(2*x);
	}

		double Cosinus2::D4(double x)
	{
		return 8*cos(2*x);
	}
	
double Cosinus2::get_value(double* x)
	 {
		 double result, temp=0;
		 int i;
		 for(i=1;i<sizeof(x)+1;i++)
			temp+=x[i];
		result=cos(temp);
		 return result;
	 }
	 
double * Cosinus2::get_value(double* x, int d)
	 {
		 double temp=0;
		 double * result;
		 result=(double*)calloc(d,sizeof(double));
		 int i;
		 for(i=1;i<sizeof(x)+1;i++)
			result[i]=cos(x[i]);
		 return result;
	 }
double * Cosinus2::grad(double*)
	 {
		 double* result;
		 return result;
	 }
double** Cosinus2::laplacian(double*)
	 {
		 double** result;
		 return result;
	 }
double*** Cosinus2::D3(double*)
	 {
		 double*** result;
		 return result;
	 }

