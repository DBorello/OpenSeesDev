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
                                                                        
// $Revision: 1.1 $
// $Date: 2001-11-26 22:55:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dPointLoad.h,v $
                                                                        
#ifndef Beam2dPointLoad_h
#define Beam2dPointLoad_h

// Written: fmk 

// Purpose: This file contains the class definition for Beam2dPointLoad.

#include <ElementalLoad.h>

class Beam2dPointLoad : public ElementalLoad
{
  public:
    Beam2dPointLoad(int tag, double P, double x, const ID &theElementTags);
    Beam2dPointLoad();    
    ~Beam2dPointLoad();

    const Vector &getData(int &type, double loadFactor);

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker);
    void Print(ostream &s, int flag =0);       

  protected:
	
  private:
    double P;     // magnitude of the load
    double x;     // relative distance (x/L) along length from end 1 of element
    static Vector data;
};

#endif

