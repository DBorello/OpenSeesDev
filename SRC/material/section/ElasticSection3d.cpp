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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticSection3d.cpp,v $
                                                                        
                                                                        
///////////////////////////////////////////////////////
// File:  ~/Src/element/hinge/ElasticSection3d.cpp
//
// Written: MHS
// Date: May 2000
//
//
// Purpose:  This file contains the function definitions
// for the ElasticSection3d class.

#include <ElasticSection3d.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <stdlib.h>

#include <classTags.h>

Vector ElasticSection3d::s(4);
ID ElasticSection3d::code(4);

ElasticSection3d::ElasticSection3d(void)
:SectionForceDeformation(0, SEC_TAG_Elastic3d),
 E(0), A(0), Iz(0), Iy(0), G(0), J(0),
 k(4,4), f(4,4), e(4), eCommit(4)
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
 k(4,4), f(4,4), e(4), eCommit(4)
{
    if (E <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input E <= 0.0 ... setting E to 1.0",
			"ElasticSection3d::ElasticSection3d");
		E = 1.0;
    }
    
    if (A <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input A <= 0.0 ... setting A to 1.0",
			"ElasticSection3d::ElasticSection3d");
		A = 1.0;
    }

    if (Iz <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input Iz <= 0.0 ... setting Iz to 1.0",
			"ElasticSection3d::ElasticSection3d");
		Iz = 1.0;
    }
    
    if (Iy <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input Iy <= 0.0 ... setting Iy to 1.0",
			"ElasticSection3d::ElasticSection3d");
		Iy = 1.0;
    }

    if (G <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input G <= 0.0 ... setting G to 1.0",
			"ElasticSection3d::ElasticSection3d");
		G = 1.0;
    }
    
    if (J <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input J <= 0.0 ... setting J to 1.0",
			"ElasticSection3d::ElasticSection3d");
		J = 1.0;
    }
    
    k(0,0) = E*A;
    k(1,1) = E*Iz;
    k(2,2) = E*Iy;
    k(3,3) = G*J;
    
    f(0,0) = 1/k(0,0);
    f(1,1) = 1/k(1,1);
    f(2,2) = 1/k(2,2);
    f(3,3) = 1/k(3,3);

    if (code(0) != SECTION_RESPONSE_P)
    {
	code(0) = SECTION_RESPONSE_P;	// P is the first quantity
	code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
	code(2) = SECTION_RESPONSE_MY;	// My is the third 
	code(3) = SECTION_RESPONSE_T;	// T is the fourth
    }
}

ElasticSection3d::ElasticSection3d
(int tag, double EA_in, double EIz_in, double EIy_in, double GJ_in)
:SectionForceDeformation(tag, SEC_TAG_Elastic3d),
 E(1), A(EA_in), Iz(EIz_in), Iy(EIy_in), G(1), J(GJ_in),
 k(4,4), f(4,4), e(4), eCommit(4)
{
    if (A <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input EA <= 0.0 ... setting EA to 1.0",
			"ElasticSection3d::ElasticSection3d");
		A = 1.0;
    }

    if (Iz <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input EIz <= 0.0 ... setting EIz to 1.0",
			"ElasticSection3d::ElasticSection3d");
		Iz = 1.0;
    }
    
    if (Iy <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input EIy <= 0.0 ... setting EIy to 1.0",
			"ElasticSection3d::ElasticSection3d");
		Iy = 1.0;
    }

    if (J <= 0.0)  {
		g3ErrorHandler->warning("%s -- Input GJ <= 0.0 ... setting GJ to 1.0",
			"ElasticSection3d::ElasticSection3d");
		J = 1.0;
    }
    
    k(0,0) = A;
    k(1,1) = Iz;
    k(2,2) = Iy;
    k(3,3) = J;
    
    f(0,0) = 1/A;
    f(1,1) = 1/Iz;
    f(2,2) = 1/Iy;
    f(3,3) = 1/J;

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
    // s = k*e;
    s(0) = k(0,0)*e(0);
    s(1) = k(1,1)*e(1);
    s(2) = k(2,2)*e(2);
    s(3) = k(3,3)*e(3);

    return s;
}

const Matrix &
ElasticSection3d::getSectionTangent(void)
{
    return k;
}

const Matrix &
ElasticSection3d::getSectionFlexibility (void)
{
    return f;
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
ElasticSection3d::getType () const
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
		g3ErrorHandler->warning("%s -- failed to send data",
			"ElasticSection3d::sendSelf");
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
		g3ErrorHandler->warning("%s -- failed to receive data",
			"ElasticSection3d::recvSelf");
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

	// Set section stiffness and flexibility matrices
	k(0,0) = E*A;
	k(1,1) = E*Iz;
	k(2,2) = E*Iy;
	k(3,3) = G*J;

	f(0,0) = 1/k(0,0);
	f(1,1) = 1/k(1,1);
	f(2,2) = 1/k(2,2);
	f(3,3) = 1/k(3,3);

    return res;
}
 
void
ElasticSection3d::Print(ostream &s, int flag)
{
    s << "ElasticSection3d, tag: " << this->getTag() << endl;
    s << "\t E: " << E << endl;
	s << "\t A: " << A << endl;
	s << "\tIz: " << Iz << endl;
	s << "\tIy: " << Iy << endl;
	s << "\t G: " << G << endl;
	s << "\t J: " << J << endl;
}
