/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2002-06-26 23:00:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/FedeasConcr2Material.cpp,v $
                                                                      
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasConcr2Material. FedeasConcr2Material wraps the FEDEAS
// 1d material subroutine Concr_2.

#include <FedeasConcr2Material.h>

FedeasConcr2Material::FedeasConcr2Material(int tag,
					 double fc, double ec, double fu, double eu,
					 double ratio, double ft, double Ets):
// 2 history variables and 7 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasConcrete2, 2, 7)
{
	data[0]  = fc;
	data[1]  = ec;
	data[2]  = fu;
	data[3]  = eu;
	data[4]  = ratio;
	data[5]  = ft;
	data[6]  = Ets;
}

FedeasConcr2Material::FedeasConcr2Material(int tag, const Vector &d):
// 2 history variables and 7 material parameters
FedeasMaterial(tag, MAT_TAG_FedeasConcrete2, 2, 7)
{
	if (d.Size() != numData)
		g3ErrorHandler->fatal("%s -- not enough input arguments",
		"FedeasConcr2Material::FedeasConcr2Material");

	for (int i = 0; i < numData; i++)
		data[i] = d(i);
}

FedeasConcr2Material::FedeasConcr2Material(void):
FedeasMaterial(0, MAT_TAG_FedeasConcrete2, 2, 7)
{
	// Does nothing
}

FedeasConcr2Material::~FedeasConcr2Material(void)
{
	// Does nothing
}

UniaxialMaterial*
FedeasConcr2Material::getCopy(void)
{
  Vector d(data, numData);

  FedeasConcr2Material *theCopy = new FedeasConcr2Material(this->getTag(), d);
  
  // Copy history variables
  for (int i = 0; i < 2*numHstv; i++)
    theCopy->hstv[i] = hstv[i];
  
  theCopy->epsilonP = epsilonP;
  theCopy->sigmaP   = sigmaP;
  theCopy->tangentP = tangentP;
  
  return theCopy;
}

double
FedeasConcr2Material::getInitialTangent(void)
{
	//return 2.0*fc/ec;
	return 2.0*data[0]/data[1];
}
