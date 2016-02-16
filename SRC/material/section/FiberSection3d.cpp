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
                                                                        
// $Revision: 1.31 $
// $Date: 2009-09-28 22:48:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSection3d.cpp,v $
                                                                        
// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSection2d.

#include <stdlib.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <FiberSection3d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <iostream>
#include <fstream>
using std::string;
using namespace std;


ID FiberSection3d::code(3);

// constructors:
FiberSection3d::FiberSection3d(int tag, int num, Fiber **fibers): 
  SectionForceDeformation(tag, SEC_TAG_FiberSection3d),
  numFibers(num), theMaterials(0), matData(0),
  yBar(0.0), zBar(0.0), e(4), eCommit(4), s(0), ks(0)
{
  if (numFibers != 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSection3d::FiberSection3d -- failed to allocate Material pointers\n";
      exit(-1);
    }

    matData = new double [numFibers*4];

    if (matData == 0) {
      opserr << "FiberSection3d::FiberSection3d -- failed to allocate double array for material data\n";
      exit(-1);
    }

    double Qz = 0.0;
    double Qy = 0.0;
    double A  = 0.0;
	double Heightt;
    
    for (int i = 0; i < numFibers; i++) {
      Fiber *theFiber = fibers[i];
      double yLoc, zLoc, Area;
      theFiber->getFiberLocation(yLoc, zLoc);
      Area = theFiber->getArea();
	  Heightt = theFiber->getd();

      Qz += yLoc*Area;
      Qy += zLoc*Area;
      A  += Area;

      matData[i*4] = yLoc;
      matData[i*4+1] = zLoc;
      matData[i*4+2] = Area;
	  matData[i*4+3] = Heightt;
      UniaxialMaterial *theMat = theFiber->getMaterial();
      theMaterials[i] = theMat->getCopy();

      if (theMaterials[i] == 0) {
	opserr << "FiberSection3d::FiberSection3d -- failed to get copy of a Material\n";
	exit(-1);
      }
    }

    yBar = -Qz/A;
    zBar = Qy/A;
  }

  s = new Vector(sData, 5);
  ks = new Matrix(kData, 5, 5);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  sData[3] = 0.0;
  sData[4] = 0.0;

  for (int i=0; i<25; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;

 // AddingSensitivity:BEGIN ////////////////////////////////////
  parameterID = 0;
  SHVs=0;
  // AddingSensitivity:END //////////////////////////////////////
}

// constructor for blank object that recvSelf needs to be invoked upon
FiberSection3d::FiberSection3d():
  SectionForceDeformation(0, SEC_TAG_FiberSection3d),
  numFibers(0), theMaterials(0), matData(0),
  yBar(0.0), zBar(0.0), e(4), eCommit(4), s(0), ks(0)
{
  s = new Vector(sData, 5);
  ks = new Matrix(kData, 5, 5);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  sData[3] = 0.0;
  sData[4] = 0.0;

  for (int i=0; i<25; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;

 // AddingSensitivity:BEGIN ////////////////////////////////////
  parameterID = 0;
  SHVs=0;
  // AddingSensitivity:END //////////////////////////////////////
}

int
FiberSection3d::addFiber(Fiber &newFiber)
{
  // need to create a larger array
  int newSize = numFibers+1;

  UniaxialMaterial **newArray = new UniaxialMaterial *[newSize]; 
  double *newMatData = new double [4 * newSize];
  
  if (newArray == 0 || newMatData == 0) {
    opserr << "FiberSection3d::addFiber -- failed to allocate Fiber pointers\n";
    exit(-1);
  }

  // copy the old pointers
  int i;
  for (i = 0; i < numFibers; i++) {
    newArray[i] = theMaterials[i];
    newMatData[4*i] = matData[4*i];
    newMatData[4*i+1] = matData[4*i+1];
    newMatData[4*i+2] = matData[4*i+2];
	newMatData[4*i+3] = matData[4*i+3];
  }
  // set the new pointers
  double yLoc, zLoc, Area, Height;
  newFiber.getFiberLocation(yLoc, zLoc);
  Area = newFiber.getArea();
  Height = newFiber.getd();
  newMatData[numFibers*4] = -yLoc;
  newMatData[numFibers*4+1] = zLoc;
  newMatData[numFibers*4+2] = Area;
  newMatData[numFibers*4+3] = Height;
  UniaxialMaterial *theMat = newFiber.getMaterial();
  newArray[numFibers] = theMat->getCopy();

  if (newArray[numFibers] == 0) {
    opserr << "FiberSection3d::addFiber -- failed to get copy of a Material\n";
    exit(-1);

    delete [] newArray;
    delete [] newMatData;
    return -1;
  }

  numFibers++;
  
  if (theMaterials != 0) {
    delete [] theMaterials;
    delete [] matData;
  }

  theMaterials = newArray;
  matData = newMatData;

  double Qz = 0.0;
  double Qy = 0.0;
  double A  = 0.0;

  // Recompute centroid
  for (i = 0; i < numFibers; i++) {
    yLoc = -matData[4*i];
    zLoc = matData[4*i+1];
    Area = matData[4*i+2];
	Height = matData[4*i+3];
    A  += Area;
    Qz += yLoc*Area;
    Qy += zLoc*Area;
  }

  yBar = -Qz/A;
  zBar = Qy/A;

  return 0;
}



// destructor:
FiberSection3d::~FiberSection3d()
{
  if (theMaterials != 0) {
    for (int i = 0; i < numFibers; i++)
      if (theMaterials[i] != 0)
	delete theMaterials[i];
      
    delete [] theMaterials;
  }

  if (matData != 0)
    delete [] matData;

  if (s != 0)
    delete s;

  if (ks != 0)
    delete ks;
}

int
FiberSection3d::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;
  e = deforms;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  kData[4] = 0.0; kData[5] = 0.0; kData[6] = 0.0; kData[7] = 0.0;
  kData[8] = 0.0; kData[9] = 0.0; kData[10] = 0.0; kData[11] = 0.0;
  kData[12] = 0.0; kData[13] = 0.0; kData[14] = 0.0; kData[15] = 0.0;
  kData[16] = 0.0; kData[17] = 0.0; kData[18] = 0.0; kData[19] = 0.0;
  kData[20] = 0.0; kData[21] = 0.0; kData[22] = 0.0; kData[23] = 0.0; kData[24] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0; sData[3] = 0.0; sData[4] = 0.0;

  int loc = 0;

  double d0 = deforms(0);
  double d1 = deforms(1);
  double d2 = deforms(2);
  double d3 = deforms(3);
  double d4 = deforms(4);
  double d5 = deforms(5);
  double d6 = deforms(6);
  double d7 = deforms(7);
  
  

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];
	double Height=matData[loc++];

	// calculate sectorial area
	double omig=0.0;
	if (y>0.0)
		omig = -z*(y-Height);
	else
		omig = -z* (y+Height);


    // determine material strain and set it, include second order terms
    double strain = d0 - y*d1 - z*d2 - omig*d3 + 0.5*d5*d5 + 0.5*d6*d6 + 0.5*(y*y+z*z)*d4*d4 - y*d7*d2 + z*d7*d1;
    double tangent, stress;
    res += theMat->setTrial(strain, stress, tangent);

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

	// section stiffness matrix k, refer to Alemdar

    kData[0] += value;
	kData[3] += (y*y+z*z)*value;
	kData[6] += vas1 * y;
	kData[12] += vas2 * z; 
	kData[15] += (y*y+z*z)*value;
	kData[18] += (y*y+z*z)*(y*y+z*z)*value;
	kData[24] += omig*omig*value;

    double fs0 = stress * A;

    sData[0] += fs0;
    sData[1] += -1.0 * fs0 * y;
    sData[2] += -1.0 * fs0 * z;
	sData[3] += fs0 * (y*y+z*z);
	sData[4] += -fs0 * omig;
  }

  return res;
}

const Matrix&
FiberSection3d::getInitialTangent(void)
{
  static double kInitialData[25];
  static Matrix kInitial(kInitialData, 5, 5);
  for (int i=0; i<25; i++)
	  kInitialData[i]=0.0;
  

  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];
	double Height = matData[loc++];
	// calculate sectorial area
	double omig;
	if (y>0.0)
		omig = -z*(y-Height);
	else
		omig = -z* (y+Height);

    double tangent = theMat->getInitialTangent();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kInitialData[0] += value;
	kInitialData[3] += (y*y+z*z)*value;
	kInitialData[6] += vas1 * y;
	kInitialData[12] += vas2 * z; 
	kInitialData[15] += (y*y+z*z)*value;
	kInitialData[18] += (y*y+z*z)*(y*y+z*z)*value;
	kInitialData[24] += omig*omig*value;
  }


  return kInitial;
}

const Vector&
FiberSection3d::getSectionDeformation(void)
{
  return e;
}

const Matrix&
FiberSection3d::getSectionTangent(void)
{
  return *ks;
}

const Vector&
FiberSection3d::getStressResultant(void)
{
  return *s;
}

SectionForceDeformation*
FiberSection3d::getCopy(void)
{
  FiberSection3d *theCopy = new FiberSection3d ();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new UniaxialMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr << "FiberSection3d::FiberSection3d -- failed to allocate Material pointers\n";
      exit(-1);			    
    }

    theCopy->matData = new double [numFibers*4];

    if (theCopy->matData == 0) {
      opserr << "FiberSection3d::FiberSection3d -- failed to allocate double array for material data\n";
      exit(-1);
    }
			    
    
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*4] = matData[i*4];
      theCopy->matData[i*4+1] = matData[i*4+1];
      theCopy->matData[i*4+2] = matData[i*4+2];
	  theCopy->matData[i*4+3] = matData[i*4+3];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == 0) {
	opserr << "FiberSection3d::getCopy -- failed to get copy of a Material\n";
	exit(-1);
      }
    }    
  }

  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;

  for (int i=0; i<25; i++)
    theCopy->kData[i] = kData[i];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];
  theCopy->sData[2] = sData[2];
  theCopy->sData[3] = sData[3];
  theCopy->sData[4] = sData[4];
  return theCopy;
}

const ID&
FiberSection3d::getType ()
{
  return code;
}

int
FiberSection3d::getOrder () const
{
  return 5;
}

int
FiberSection3d::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  eCommit = e;

  return err;
}

int
FiberSection3d::revertToLastCommit(void)
{
  int err = 0;

  // Last committed section deformations
  e = eCommit;


  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  kData[4] = 0.0; kData[5] = 0.0; kData[6] = 0.0; kData[7] = 0.0;
  kData[8] = 0.0; kData[9] = 0.0; kData[10] = 0.0; kData[11] = 0.0;
  kData[12] = 0.0; kData[13] = 0.0; kData[14] = 0.0; kData[15] = 0.0;
  kData[16] = 0.0; kData[17] = 0.0; kData[18] = 0.0; kData[19] = 0.0;
  kData[20] = 0.0; kData[21] = 0.0; kData[22] = 0.0; kData[23] = 0.0; kData[24] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0; sData[3] = 0.0; sData[4] = 0.0;


  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];
	double Height = matData[loc++];

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;
	double omig;
	if (y>0.0)
		omig = -z*(y-Height);
	else
		omig = -z* (y+Height);

    kData[0] += value;
	kData[3] += (y*y+z*z)*value;
	kData[6] += vas1 * y;
	kData[12] += vas2 * z; 
	kData[15] += (y*y+z*z)*value;
	kData[18] += (y*y+z*z)*(y*y+z*z)*value;
	kData[24] += omig*omig*value;
    
    double fs0 = stress * A;
    
    sData[0] += fs0;
    sData[1] += -1.0 * fs0 * y;
    sData[2] += -1.0 * fs0 * z;
	sData[3] += fs0 * (y*y+z*z);
	sData[4] += -fs0 * omig;
  }

  return err;
}

int
FiberSection3d::revertToStart(void)
{
  // revert the fibers to start    
  int err = 0;


  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  kData[4] = 0.0; kData[5] = 0.0; kData[6] = 0.0; kData[7] = 0.0;
  kData[8] = 0.0; kData[9] = 0.0; kData[10] = 0.0; kData[11] = 0.0;
  kData[12] = 0.0; kData[13] = 0.0; kData[14] = 0.0; kData[15] = 0.0;
  kData[16] = 0.0; kData[17] = 0.0; kData[18] = 0.0; kData[19] = 0.0;
  kData[20] = 0.0; kData[21] = 0.0; kData[22] = 0.0; kData[23] = 0.0; kData[24] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0; sData[3] = 0.0; sData[4] = 0.0;

  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];
	double Height = matData[loc++];
	double omig;
	if (y>0.0)
		omig = -z*(y-Height);
	else
		omig = -z* (y+Height);

    // invoke revertToStart on the material
    err += theMat->revertToStart();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kData[0] += value;
	kData[3] += (y*y+z*z)*value;
	kData[6] += vas1 * y;
	kData[12] += vas2 * z; 
	kData[15] += (y*y+z*z)*value;
	kData[18] += (y*y+z*z)*(y*y+z*z)*value;
	kData[24] += omig*omig*value;
    
    double fs0 = stress * A;

    sData[0] += fs0;
    sData[1] += -1.0 * fs0 * y;
    sData[2] += -1.0 * fs0 * z;
	sData[3] += fs0 * (y*y+z*z);
	sData[4] += -fs0 * omig;
  }


  return err;
}

int
FiberSection3d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and numFibers, 
  //     size 3 so no conflict with matData below if just 1 fiber
  static ID data(3);
  data(0) = this->getTag();
  data(1) = numFibers;
  int dbTag = this->getDbTag();
  res += theChannel.sendID(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "FiberSection2d::sendSelf - failed to send ID data\n";
    return res;
  }    

  if (numFibers != 0) {
    
    // create an id containingg classTag and dbTag for each material & send it
    ID materialData(2*numFibers);
    for (int i=0; i<numFibers; i++) {
      UniaxialMaterial *theMat = theMaterials[i];
      materialData(2*i) = theMat->getClassTag();
      int matDbTag = theMat->getDbTag();
      if (matDbTag == 0) {
	matDbTag = theChannel.getDbTag();
	if (matDbTag != 0)
	  theMat->setDbTag(matDbTag);
      }
      materialData(2*i+1) = matDbTag;
    }    
    
    res += theChannel.sendID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "FiberSection2d::sendSelf - failed to send material data\n";
     return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 4*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "FiberSection2d::sendSelf - failed to send material data\n";
     return res;
    }    

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);
  }

  return res;
}

int
FiberSection3d::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(3);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);

  if (res < 0) {
   opserr << "FiberSection2d::sendSelf - failed to recv ID data\n";
   return res;
  } 
   
  this->setTag(data(0));

  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "FiberSection2d::sendSelf - failed to send material data\n";
     return res;
    }    

    // if current arrays not of correct size, release old and resize
    if (theMaterials == 0 || numFibers != data(1)) {
      // delete old stuff if outa date
      if (theMaterials != 0) {
	for (int i=0; i<numFibers; i++)
	  delete theMaterials[i];
	delete [] theMaterials;
	if (matData != 0)
	  delete [] matData;
	matData = 0;
	theMaterials = 0;
      }

      // create memory to hold material pointers and fiber data
      numFibers = data(1);
      if (numFibers != 0) {

	theMaterials = new UniaxialMaterial *[numFibers];
	
	if (theMaterials == 0) {
	  opserr << "FiberSection2d::recvSelf -- failed to allocate Material pointers\n";
	  exit(-1);
	}

	for (int j=0; j<numFibers; j++)
	  theMaterials[j] = 0;
	
	matData = new double [numFibers*4];

	if (matData == 0) {
	  opserr << "FiberSection2d::recvSelf  -- failed to allocate double array for material data\n";
	  exit(-1);
	}
      }
    }

    Vector fiberData(matData, 4*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "FiberSection2d::sendSelf - failed to send material data\n";
     return res;
    }    
    
    int i;
    for (i=0; i<numFibers; i++) {
      int classTag = materialData(2*i);
      int dbTag = materialData(2*i+1);

      // if material pointed to is blank or not of corrcet type, 
      // release old and create a new one
      if (theMaterials[i] == 0)
	theMaterials[i] = theBroker.getNewUniaxialMaterial(classTag);
      else if (theMaterials[i]->getClassTag() != classTag) {
	delete theMaterials[i];
	theMaterials[i] = theBroker.getNewUniaxialMaterial(classTag);      
      }

      if (theMaterials[i] == 0) {
	opserr << "FiberSection2d::recvSelf -- failed to allocate double array for material data\n";
	exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    double Qz = 0.0;
    double Qy = 0.0;
    double A  = 0.0;
    double yLoc, zLoc, Area, Height;

    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
      yLoc = -matData[4*i];
      zLoc = matData[4*i+1];
      Area = matData[4*i+2];
	  Height=matData[4*i+3];
      A  += Area;
      Qz += yLoc*Area;
      Qy += zLoc*Area;
    }
    
    yBar = -Qz/A;
    zBar = Qy/A;
  }    

  return res;
}

void
FiberSection3d::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {
    for (int i = 0; i < numFibers; i++) {
      s << -matData[4*i] << " "  << matData[4*i+1] << " "  << matData[4*i+2] << " " ;
      s << theMaterials[i]->getStress() << " "  << theMaterials[i]->getStrain() << endln;
    } 
  } else {
    s << "\nFiberSection3d, tag: " << this->getTag() << endln;
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << numFibers << endln;
    s << "\tCentroid: (" << -yBar << ", " << zBar << ')' << endln;
    
    if (flag == 1) {
      for (int i = 0; i < numFibers; i++) {
	s << "\nLocation (y, z) = (" << -matData[4*i] << ", " << matData[4*i+1] << ")";
	s << "\nArea = " << matData[4*i+2] << endln;
      theMaterials[i]->Print(s, flag);
      }
    }
  }
}

Response*
FiberSection3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  const ID &type = this->getType();
  int typeSize = this->getOrder();

  Response *theResponse =0;

  output.tag("SectionOutput");
  output.attr("secType", this->getClassType());
  output.attr("secTag", this->getTag());

  // deformations
  if (strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"deformation") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","kappaZ");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","eps");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","gammaY");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","kappaY");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","gammaZ");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","theta");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new MaterialResponse(this, 1, this->getSectionDeformation());
  
  // forces
  } else if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","Mz");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","P");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","Vy");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","My");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","Vz");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","T");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new MaterialResponse(this, 2, this->getStressResultant());
  
  // force and deformation
  } else if (strcmp(argv[0],"forceAndDeformation") == 0) { 
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","kappaZ");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","eps");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","gammaY");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","kappaY");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","gammaZ");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","theta");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    for (int j=0; j<typeSize; j++) {
      int code = type(j);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","Mz");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","P");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","Vy");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","My");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","Vz");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","T");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }

    theResponse =  new MaterialResponse(this, 4, Vector(2*this->getOrder()));
  
  }  
  
  else {
    if (argc > 2 || strcmp(argv[0],"fiber") == 0) {

      int key = numFibers;
      int passarg = 2;
      
      
      if (argc <= 3)	{  // fiber number was input directly
	
	key = atoi(argv[1]);
	
      } else if (argc > 4) {         // find fiber closest to coord. with mat tag
	int matTag = atoi(argv[3]);
	double yCoord = atof(argv[1]);
	double zCoord = atof(argv[2]);
	double closestDist;
	double ySearch, zSearch, dy, dz;
	double distance;
	int j;
	
	// Find first fiber with specified material tag
	for (j = 0; j < numFibers; j++) {
	  if (matTag == theMaterials[j]->getTag()) {
	    ySearch = -matData[4*j];
	    zSearch =  matData[4*j+1];
	    dy = ySearch-yCoord;
	    dz = zSearch-zCoord;
	    closestDist = sqrt(dy*dy + dz*dz);
	    key = j;
	    break;
	  }
	}
	
	// Search the remaining fibers
	for ( ; j < numFibers; j++) {
	  if (matTag == theMaterials[j]->getTag()) {
	    ySearch = -matData[4*j];
	    zSearch =  matData[4*j+1];
	    dy = ySearch-yCoord;
	    dz = zSearch-zCoord;
	    distance = sqrt(dy*dy + dz*dz);
	    if (distance < closestDist) {
	      closestDist = distance;
	      key = j;
	    }
	  }
	}
	passarg = 4;
      }
      
      else {                  // fiber near-to coordinate specified
	double yCoord = atof(argv[1]);
	double zCoord = atof(argv[2]);
	double closestDist;
	double ySearch, zSearch, dy, dz;
	double distance;
	ySearch = -matData[0];
	zSearch =  matData[1];
	dy = ySearch-yCoord;
	dz = zSearch-zCoord;
	closestDist = sqrt(dy*dy + dz*dz);
	key = 0;
	for (int j = 1; j < numFibers; j++) {
	  ySearch = -matData[4*j];
	  zSearch =  matData[4*j+1];
	  dy = ySearch-yCoord;
	  dz = zSearch-zCoord;
	  distance = sqrt(dy*dy + dz*dz);
	  if (distance < closestDist) {
	    closestDist = distance;
	    key = j;
	  }
	}
	passarg = 3;
      }
      
      if (key < numFibers && key >= 0) {
	output.tag("FiberOutput");
	output.attr("yLoc",-matData[4*key]);
	output.attr("zLoc",matData[4*key+1]);
	output.attr("area",matData[4*key+2]);
	
	theResponse =  theMaterials[key]->setResponse(&argv[passarg], argc-passarg, output);
	
	output.endTag();
      }
    }
  }

  output.endTag();
  return theResponse;
}


int 
FiberSection3d::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}

int
FiberSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 3)
    return -1;


  int result = 0;

  // A material parameter
  if (strstr(argv[0],"material") != 0) {

    // Get the tag of the material
    int paramMatTag = atoi(argv[1]);

    // Loop over fibers to find the right material(s)
    int ok = 0;
    for (int i = 0; i < numFibers; i++)
      if (paramMatTag == theMaterials[i]->getTag()) {
	ok = theMaterials[i]->setParameter(&argv[2], argc-2, param);
	if (ok != -1)
	  result = ok;
      }
    
    return result;
  }    

  int ok = 0;
  
  // loop over every material
  for (int i = 0; i < numFibers; i++) {
    ok = theMaterials[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  return result;
}

const Vector &
FiberSection3d::getSectionDeformationSensitivity(int gradIndex)
{
	static Vector dummy(3);
	dummy.Zero();
	if (SHVs !=0) {
		dummy(0) = (*SHVs)(0,gradIndex);
		dummy(1) = (*SHVs)(1,gradIndex);
		dummy(2) = (*SHVs)(2,gradIndex);
	}
	return dummy;
}

   
const Vector &
FiberSection3d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  
  static Vector ds(3);
  
  ds.Zero();
  
  double  stressGradient;
  int loc = 0;
  
  
  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];
    stressGradient = theMaterials[i]->getStressSensitivity(gradIndex,conditional);
    stressGradient *=  A;
    ds(0) += stressGradient;
    ds(1) += stressGradient * y;
    ds(2) += stressGradient * z;
    
  } 
  
  return ds;
}

const Matrix &
FiberSection3d::getSectionTangentSensitivity(int gradIndex)
{
  static Matrix something(2,2);
  
  something.Zero();
  
  return something;
}

int
FiberSection3d::commitSensitivity(const Vector& defSens, int gradIndex, int numGrads)
{

  // here add SHVs to store the strain sensitivity.

  if (SHVs == 0) {
    SHVs = new Matrix(3,numGrads);
  }
  
  (*SHVs)(0,gradIndex) = defSens(0);
  (*SHVs)(1,gradIndex) = defSens(1);
  (*SHVs)(2,gradIndex) = defSens(2);

  int loc = 0;

  double d0 = defSens(0);
  double d1 = defSens(1);
  double d2 = defSens(2);

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
   	double y = matData[loc++] - yBar;
	double z = matData[loc++] - zBar;
	loc++;   // skip A data.

	double strainSens = d0 + y*d1 + z*d2;


    
	theMat->commitSensitivity(strainSens,gradIndex,numGrads);
  }

  return 0;
}

// AddingSensitivity:END ///////////////////////////////////


