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
                                                                        
// $Revision: 1.0 $
// $Date: 2012-09-17 10:51:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMSolver.h,v $
                                                                        
                                                                        
#ifndef PFEMSolver_h
#define PFEMSolver_h

// File: ~/system_of_eqn/linearSOE/sparseGEN/PFEMSolver.h
//
// Written: Minjie 
// Created: Sep 17 2012
//
// Description: This file contains the class definition for PFEMSolver.
// A PFEMSolver object can be constructed to solve a PFEMLinSOE
// object. It obtains the solution by making calls on the
// the PFEMSolver library developed at UC Berkeley by Prof. James Demmel, 
// Xiaoye S. Li and John R. Gilbert.
// The PFEMSolver uses Fractional Step Method to solve PFEM equations. 
//
// What: "@(#) PFEMSolver.h, revA"

#include <SparseGenColLinSolver.h>
extern "C" {
#include <cs.h>
}

class ID;


class PFEMSolver : public SparseGenColLinSolver
{
public:
    PFEMSolver();
    virtual ~PFEMSolver();

    int solve();
    int setSize();

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);  

    virtual int setDofIDs();
    virtual cs* setMatIDs();

protected:


private:
    ID* mID, *pID, *piID, *mIDall;
    ID* Mid, *Mhatid, *Gid, *Gtid, *Lid, *Qtid;
    cs* G, *Gt, *L, *Qt;
};

#endif

