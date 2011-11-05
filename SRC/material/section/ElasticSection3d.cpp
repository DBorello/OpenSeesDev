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
                                                                        
// $Revision: 1.10 $
// $Date: 2008-08-26 16:46:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticSection3d.cpp,v $

#include <ElasticSection3d.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <Parameter.h>
#include <stdlib.h>
#include <string.h>

#include <classTags.h>

Vector ElasticSection3d::s(4);
Matrix ElasticSection3d::ks(4,4);
ID ElasticSection3d::code(4);

ElasticSection3d::ElasticSection3d(void)
:SectionForceDeformation(0, SEC_TAG_Elastic3d),
 E(0.0), A(0.0), Iz(0.0), Iy(0.0), G(0.0), J(0.0),
 e(4), eCommit(4)
{
    if (code(0) != SECTION_RESPONSE_P)
    {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
	code(2) = SECTION_RESPONSE_MY;	// My is the third 
	code(3) = SECTION_RESPONSE_T;	// T is the fourth
    }
}

ElasticSection3d::ElasticSection3d
(int tag, double E_in, double A_in, double Iz_in, double Iy_in, double G_in, double J_in)
:SectionForceDeformation(tag, SEC_TAG_Elastic3d),
 E(E_in), A(A_in), Iz(Iz_in), Iy(Iy_in), G(G_in), J(J_in),
 e(4), eCommit(4)
{
    if (E <= 0.0)  {
      opserr << "ElasticSection3d::ElasticSection3d -- Input E <= 0.0 ... setting E to 1.0\n";
      E = 1.0;
    }
    
    if (A <= 0.0)  {
      opserr << "ElasticSection3d::ElasticSection3d -- Input A <= 0.0 ... setting A to 1.0\n";
      A = 1.0;
    }

    if (Iz <= 0.0)  {
      opserr << "ElasticSection3d::ElasticSection3d -- Input Iz <= 0.0 ... setting Iz to 1.0\n";
      Iz = 1.0;
    }
    
    if (Iy <= 0.0)  {
      opserr << "ElasticSection3d::ElasticSection3d -- Input Iy <= 0.0 ... setting Iy to 1.0\n";
      Iy = 1.0;
    }

    if (G <= 0.0)  {
      opserr << "ElasticSection3d::ElasticSection3d -- Input G <= 0.0 ... setting G to 1.0\n";
      G = 1.0;
    }
    
    if (J <= 0.0)  {
      opserr << "ElasticSection3d::ElasticSection3d -- Input J <= 0.0 ... setting J to 1.0\n";
      J = 1.0;
    }
    
    if (code(0) != SECTION_RESPONSE_P)
    {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
	code(2) = SECTION_RESPONSE_MY;	// My is the third 
	code(3) = SECTION_RESPONSE_T;	// T is the fourth
    }
}

ElasticSection3d::~ElasticSection3d(void)
{
    return;
}

int 
ElasticSection3d::commitState(void)
{
	eCommit = e;

    return 0;
}

int 
ElasticSection3d::revertToLastCommit(void)
{
	e = eCommit;

    return 0;
}

int 
ElasticSection3d::revertToStart(void)
{
	eCommit.Zero();

    return 0;
}

int
ElasticSection3d::setTrialSectionDeformation (const Vector &def)
{
    e = def;
    
	return 0;
}

const Vector &
ElasticSection3d::getSectionDeformation (void)
{
    return e;
}

const Vector &
ElasticSection3d::getStressResultant (void)
{
  s(0) = E*A*e(0);
  s(1) = E*Iz*e(1);
  s(2) = E*Iy*e(2);
  s(3) = G*J*e(3);
  
  return s;
}

const Matrix &
ElasticSection3d::getSectionTangent(void)
{
  ks(0,0) = E*A;
  ks(1,1) = E*Iz;
  ks(2,2) = E*Iy;
  ks(3,3) = G*J;
  
  return ks;
}

const Matrix &
ElasticSection3d::getInitialTangent(void)
{
  ks(0,0) = E*A;
  ks(1,1) = E*Iz;
  ks(2,2) = E*Iy;
  ks(3,3) = G*J;
  
  return ks;
}

const Matrix &
ElasticSection3d::getSectionFlexibility (void)
{
  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*Iz);
  ks(2,2) = 1.0/(E*Iy);
  ks(3,3) = 1.0/(G*J);
  
  return ks;
}

const Matrix &
ElasticSection3d::getInitialFlexibility (void)
{
  ks(0,0) = 1.0/(E*A);
  ks(1,1) = 1.0/(E*Iz);
  ks(2,2) = 1.0/(E*Iy);
  ks(3,3) = 1.0/(G*J);
  
  return ks;
}

SectionForceDeformation*
ElasticSection3d::getCopy ()
{
    // Make a copy of the hinge
    ElasticSection3d *theCopy =
	new ElasticSection3d (this->getTag(), E, A, Iz, Iy, G, J);

    theCopy->eCommit = eCommit;

    return theCopy;
}

const ID&
ElasticSection3d::getType ()
{
    return code;
}

int
ElasticSection3d::getOrder () const
{
    return 4;
}

int
ElasticSection3d::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(11);

    int dataTag = this->getDbTag();
    
	data(0) = this->getTag();
    data(1) = E;
    data(2) = A;    
    data(3) = Iz;
    data(4) = Iy;
    data(5) = G;
    data(6) = J;
    data(7) = eCommit(0);
	data(8) = eCommit(1);
	data(9) = eCommit(2);
	data(10) = eCommit(3);
    
    res += theChannel.sendVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticSection3d::sendSelf -- failed to send data\n";
      return res;
    }
    
    return res;
}

int
ElasticSection3d::recvSelf(int commitTag, Channel &theChannel,
					 FEM_ObjectBroker &theBroker)
{
    int res = 0;
    
	static Vector data(11);

    int dataTag = this->getDbTag();

    res += theChannel.recvVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticSection3d::recvSelf -- failed to receive data\n";
      return res;
    }

	this->setTag((int)data(0));
    E = data(1);
    A = data(2);    
    Iz = data(3);
    Iy = data(4);
    G = data(5);
    J = data(6);    
    eCommit(0) = data(7);
	eCommit(1) = data(8);
	eCommit(2) = data(9);
	eCommit(3) = data(10);

    return res;
}
 
void
ElasticSection3d::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {

  } else {
    s << "ElasticSection3d, tag: " << this->getTag() << endln;
    s << "\t E: " << E << endln;
    s << "\t A: " << A << endln;
    s << "\tIz: " << Iz << endln;
    s << "\tIy: " << Iy << endln;
    s << "\t G: " << G << endln;
    s << "\t J: " << J << endln;
  }
}

int
ElasticSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"E") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"A") == 0)
    return param.addObject(2, this);

  if (strcmp(argv[0],"Iz") == 0)
    return param.addObject(3, this);

  if (strcmp(argv[0],"Iy") == 0)
    return param.addObject(4, this);

  if (strcmp(argv[0],"G") == 0)
    return param.addObject(5, this);

  if (strcmp(argv[0],"J") == 0)
    return param.addObject(6, this);

  return -1;
}

int
ElasticSection3d::updateParameter(int paramID, Information &info)
{
  if (paramID == 1)
    E = info.theDouble;
  if (paramID == 2)
    A = info.theDouble;
  if (paramID == 3)
    Iz = info.theDouble;
  if (paramID == 4)
    Iy = info.theDouble;
  if (paramID == 5)
    G = info.theDouble;
  if (paramID == 6)
    J = info.theDouble;

  return 0;
}

int
ElasticSection3d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
ElasticSection3d::getStressResultantSensitivity(int gradIndex,
						bool conditional)
{
  s.Zero();

  if (parameterID == 1) { // E
    s(0) = A*e(0);
    s(1) = Iz*e(1);
    s(2) = Iy*e(2);
  }
  if (parameterID == 2) // A
    s(0) = E*e(0);
  if (parameterID == 3) // Iz
    s(1) = E*e(1);
  if (parameterID == 4) // Iy
    s(2) = E*e(2);
  if (parameterID == 5) // G
    s(3) = J*e(3);
  if (parameterID == 6) // J
    s(3) = G*e(3);

  return s;
}

const Matrix&
ElasticSection3d::getInitialTangentSensitivity(int gradIndex)
{
  ks.Zero();

  return ks;
}
