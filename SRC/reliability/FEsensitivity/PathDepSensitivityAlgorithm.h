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
                                                                        
// $Revision: 1.1 $
// $Date: 2001-07-31 22:11:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/PathDepSensitivityAlgorithm.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu), July 2001
// Revised: 
//

#ifndef PathDepSensitivityAlgorithm_h
#define PathDepSensitivityAlgorithm_h

#include <ReliabilityDomain.h>
#include <EquiSolnAlgo.h>
#include <SensitivityIntegrator.h>
#include <SensitivityAlgorithm.h>

class PathDepSensitivityAlgorithm : public SensitivityAlgorithm
{
  public:
	PathDepSensitivityAlgorithm::PathDepSensitivityAlgorithm(ReliabilityDomain *passedReliabilityDomain,
										   EquiSolnAlgo *passedAlgorithm,
										   SensitivityIntegrator *passedSensitivityIntegrator);
    ~PathDepSensitivityAlgorithm();
	int computeGradients(void);
    
  protected:
    
  private:
    ReliabilityDomain *theReliabilityDomain;
	EquiSolnAlgo *theAlgorithm;
	SensitivityIntegrator *theSensitivityIntegrator;
};

#endif

