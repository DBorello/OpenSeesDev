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
// $Date: 2002-12-05 22:33:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/WilsonTheta.h,v $
                                                                        
                                                                        
#ifndef WilsonTheta_h
#define WilsonTheta_h

// File: ~/analysis/integrator/WilsonTheta.h
// 
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for WilsonTheta.
// WilsonTheta is an algorithmic class for performing a transient analysis
// using the WilsonTheta integration scheme.
//
// What: "@(#) WilsonTheta.h, revA"

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class WilsonTheta : public TransientIntegrator
{
  public:
    WilsonTheta();
    WilsonTheta(double theta);
    WilsonTheta(double theta, double alphaM, double betaK, double betaKi);    
    ~WilsonTheta();

    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);        

    int domainChanged(void);    
    int newStep(double deltaT);    
    int update(const Vector &deltaU);
    int commit(void);

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(ostream &s, int flag =0);        
    
  protected:
    
  private:
    double theta;    
    double deltaT;

    double alphaM, betaK, betaKi;  
    double c1, c2, c3;  // some constants we need to keep
    Vector *Ut, *Utdot, *Utdotdot; // response quantities at time t
    Vector *U, *Udot, *Udotdot; // response quantities at time t+deltat
};

#endif

