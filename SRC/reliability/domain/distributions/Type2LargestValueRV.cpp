/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.3 $
// $Date: 2001-08-01 00:25:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/distributions/Type2LargestValueRV.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) during Spring 2000
// Revised: haukaas 06/00 (core code)
//			haukaas 06/01 (made part of official OpenSees)
//

#include <Type2LargestValueRV.h>
#include <GammaRV.h>
#include <math.h>
#include <string.h>

Type2LargestValueRV::Type2LargestValueRV(int passedTag, 
		 double passedMean,
		 double passedStdv,
		 double passedStartValue)
:RandomVariable(passedTag, passedMean, passedStdv, passedStartValue)
{
	type = new char[100];
	strcpy(type,"TYPE2LARGESTVALUE");
	tag = passedTag ;
	setParameters(passedMean,passedStdv);
	startValue = passedStartValue;
}
Type2LargestValueRV::Type2LargestValueRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4,
		 double passedStartValue)
:RandomVariable(passedTag, passedParameter1, passedParameter2, passedParameter3, passedParameter4, passedStartValue)
{
	type = new char[100];
	strcpy(type,"TYPE2LARGESTVALUE");
	tag = passedTag ;
	u = passedParameter1;
	k = passedParameter2;
	startValue = passedStartValue;
}
Type2LargestValueRV::Type2LargestValueRV(int passedTag, 
		 double passedMean,
		 double passedStdv)
:RandomVariable(passedTag, passedMean, passedStdv)
{
	type = new char[100];
	strcpy(type,"TYPE2LARGESTVALUE");
	tag = passedTag ;
	setParameters(passedMean,passedStdv);
	startValue = getMean();
}
Type2LargestValueRV::Type2LargestValueRV(int passedTag, 
		 double passedParameter1,
		 double passedParameter2,
		 double passedParameter3,
		 double passedParameter4)
:RandomVariable(passedTag, passedParameter1, passedParameter2, passedParameter3, passedParameter4)
{
	type = new char[100];
	strcpy(type,"TYPE2LARGESTVALUE");
	tag = passedTag ;
	u = passedParameter1;
	k = passedParameter2;
	startValue = getMean();
}


Type2LargestValueRV::~Type2LargestValueRV()
{
  if (type != 0)
    delete [] type;
}


void
Type2LargestValueRV::Print(ostream &s, int flag)
{
}


double
Type2LargestValueRV::getPDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = k/u * pow(u/rvValue,k+1) * exp(-pow(u/rvValue,k));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
Type2LargestValueRV::getCDFvalue(double rvValue)
{
	double result;
	if ( 0.0 < rvValue ) {
		result = exp(-pow(u/rvValue,k));
	}
	else {
		result = 0.0;
	}
	return result;
}


double
Type2LargestValueRV::getInverseCDFvalue(double probValue)
{
	return u * pow((-log(probValue)),(-1.0/k));
}


char *
Type2LargestValueRV::getType()
{
	return type;
}


double 
Type2LargestValueRV::getMean()
{
	double result;
	GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	result = u * aGammaRV->gammaFunction(1.0-1.0/k);
	delete aGammaRV;
	return result;
}



double 
Type2LargestValueRV::getStdv()
{
	GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	double a = aGammaRV->gammaFunction(1.0-2.0/k);
	double b = aGammaRV->gammaFunction(1.0-1.0/k);
	delete aGammaRV;
	double result = u*sqrt(a-b*b);
	return result;
}


double 
Type2LargestValueRV::getStartValue()
{
	return startValue;
}

double Type2LargestValueRV::getParameter1()  {return u;}
double Type2LargestValueRV::getParameter2()  {return k;}
double Type2LargestValueRV::getParameter3()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}
double Type2LargestValueRV::getParameter4()  {cerr<<"No such parameter in r.v. #"<<tag<<endl; return 0.0;}



void
Type2LargestValueRV::setParameters(double mean, double stdv)
{
	double cov = stdv/mean;
	double xk = 1.5;
	function131(xk, cov,mean);
}


void
Type2LargestValueRV::function131(double xk, double cov, double mean)
{
	GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	xk = xk + 1.0;
	double x1 = 1.0 - 1.0/xk;
	double x2 = 1.0 - 2.0/xk;
	double gm1 = aGammaRV->gammaFunction(x1);
	double gm2 = aGammaRV->gammaFunction(x2);
	delete aGammaRV;
	double vy = sqrt(gm2/gm1/gm1 - 1.0);
	if (cov-vy < 0.0) 
		function131(xk,cov,mean);
	else if (cov-vy == 0.0) 
		function134(xk,cov,mean);
	else 
		function132(xk,cov,mean);
}
void
Type2LargestValueRV::function132(double xk, double cov, double mean)
{
	xk = xk - 0.1;
	double x1 = 1.0 - 1.0/xk;
	double x2 = 1.0 - 2.0/xk;
	GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	double gm1 = aGammaRV->gammaFunction(x1);
	double gm2 = aGammaRV->gammaFunction(x2);
	delete aGammaRV;
	double vy = sqrt(gm2/gm1/gm1 - 1.0);
	if (cov-vy < 0.0) 
		function133(xk,cov,mean);
	else if (cov-vy == 0.0) 
		function134(xk,cov,mean);
	else 
		function132(xk,cov,mean);
}
void
Type2LargestValueRV::function133(double xk, double cov, double mean)
{
	xk = xk + 0.01;
	double x1 = 1.0 - 1.0/xk;
	double x2 = 1.0 - 2.0/xk;
	GammaRV *aGammaRV = new GammaRV(1, 0.0, 1.0, 0.0);
	double gm1 = aGammaRV->gammaFunction(x1);
	double gm2 = aGammaRV->gammaFunction(x2);
	delete aGammaRV;
	double vy = sqrt(gm2/gm1/gm1 - 1.0);
	if (cov-vy < 0.0) 
		function133(xk,cov,mean);
	else if (cov-vy == 0.0) 
		function134(xk,cov,mean);
	else 
		function134(xk,cov,mean);


}
void
Type2LargestValueRV::function134(double xk, double gm1, double mean)
{
	u = mean/gm1;
	k = xk;
}
