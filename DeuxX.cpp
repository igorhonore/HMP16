#include "DeuxX.h"
#include <cmath>
#include <iostream>
#include <vector>


using namespace::std;

	double DeuxX::get_value(double x)
{
	return 2*x;
}

	double DeuxX::grad(double x)
{
	return 2;
}

	double DeuxX::laplacian(double x)
{
	return 0;
}

	 double DeuxX::get_value(std::vector<double> x)
 {
	 double result=0;
	 int i,n;
	 n=sizeof(x);
	for(i=1;i<n+1;i++)
		result+=2*x[i];
	 return result;
 }
	 std::vector<double> DeuxX::grad(std::vector<double> x)
{
	 std::vector<double> result=x;
	 return result;
 }
	 std::vector<std::vector<double> > DeuxX::laplacian(std::vector<double> x)
{
	 std::vector< std::vector<double> > result;
	 return result;
 }
	
	 std::vector<std::vector<double> > DeuxX::get_value( std::vector<double>  x, int d, int r)
 {
	std::vector<std::vector<double> > result;
	int i,j,n;
	n=sizeof(x);
	result.resize(n);
	for(i=1;i<n+1;i++)
		result[i].resize(n);
	 for(i=1;i<n+1;i++)
	 {
		for(j=1;j<r+1;j++)
			result[i][j]=0;
		result[i][i]=2*x[i];
	 }
	 return result;
 }
	 
std::vector<double > DeuxX::get_value( std::vector<double>  x, int d)
 {
	 std::vector<double> result;
	 return result;
 }
