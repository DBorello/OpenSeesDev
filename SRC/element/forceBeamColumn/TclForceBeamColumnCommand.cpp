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
// $Date: 2006-09-05 23:24:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/TclForceBeamColumnCommand.cpp,v $
                                                                        
// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the implementation of the 
// TclModelBuilder_addDispBeamColumn() command. 

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <TclModelBuilder.h>

#include <ForceBeamColumn2d.h>
#include <ForceBeamColumn3d.h>

#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>
#include <UserDefinedBeamIntegration.h>

#include <HingeMidpointBeamIntegration.h>
#include <HingeEndpointBeamIntegration.h>
#include <HingeRadauBeamIntegration.h>
#include <HingeRadauTwoBeamIntegration.h>
#include <UserDefinedHingeIntegration.h>
#include <DistHingeIntegration.h>

#include <TrapezoidalBeamIntegration.h>
#include <FixedLocationBeamIntegration.h>
#include <LowOrderBeamIntegration.h>
#include <MidDistanceBeamIntegration.h>

#include <ElasticSection2d.h>
#include <ElasticSection3d.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addForceBeamColumn(ClientData clientData, Tcl_Interp *interp,  
				   int inArgc, 
				   TCL_Char **inArgv, 
				   Domain*theTclDomain,
				   TclModelBuilder *theTclBuilder)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }
  
  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();
  
  int ok = 0;
  if (ndm == 2 && ndf == 3)
    ok = 1;
  if (ndm == 3 && ndf == 6)
    ok = 1;
  
  if (ok == 0) {
    opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
	 << " not compatible with forceBeamColumn element" << endln;
    return TCL_ERROR;
  }
  
  // split possible lists present in argv
  char *List;

  List = Tcl_Merge (inArgc, inArgv);
  if (List == 0) {
    opserr << "WARNING - TclModelBuilder_addForceBeamColumn - problem merging list\n";
    return TCL_ERROR;
  }

  //  opserr << "List :" << List << endln;

  // remove braces from list
  for (int i = 0; List[i] != '\0'; i++) {
    if ((List[i] == '{')  ||  (List[i] == '}'))
      List[i] = ' ';
  }
  
  int argc;
  TCL_Char **argv;
       
  if (Tcl_SplitList(interp, List, &argc, &argv) != TCL_OK) {
    opserr <<  "WARNING - TclModelBuilder_addForceBeamColumn - problem spliting list\n";
    return TCL_ERROR;
  }
      
  Tcl_Free (List);


  if (argc < 6) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? ...\n";
    return TCL_ERROR;
  }

  int eleTag, iNode, jNode, transfTag;
  CrdTransf2d *theTransf2d = 0;
  CrdTransf3d *theTransf3d = 0;
  Element *theElement = 0;


  if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
    opserr << "WARNING invalid forceBeamColumn eleTag" << endln;
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "forceBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode\n";
    opserr << "forceBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }


  //
  // fmk UNDOCUMENTED FEATURE - 
  // all to take similar command to nonlinearBeamColumn & dispBeamColumn 
  // 

  if ((strcmp(argv[6],"Lobatto") != 0) &&
      (strcmp(argv[6],"Legendre") != 0) &&
      (strcmp(argv[6],"Radau") != 0) &&
      (strcmp(argv[6],"NewtonCotes") != 0) &&
      (strcmp(argv[6],"UserDefined") != 0) &&
      (strcmp(argv[6],"HingeMidpoint") != 0) &&
      (strcmp(argv[6],"HingeEndpoint") != 0) &&
      (strcmp(argv[6],"HingeRadau") != 0) &&
      (strcmp(argv[6],"HingeRadauTwo") != 0) &&
      (strcmp(argv[6],"UserHinge") != 0) &&
      (strcmp(argv[6],"DistHinge") != 0) &&
      (strcmp(argv[6],"Trapezoidal") != 0) &&
      (strcmp(argv[6],"FixedLocation") != 0) &&
      (strcmp(argv[6],"LowOrder") != 0) &&
      (strcmp(argv[6],"MidDistance") != 0)) {

    int nIP, secTag;

    if (Tcl_GetInt(interp, argv[5], &nIP) != TCL_OK) {
      opserr << "WARNING invalid nIP\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[6], &secTag) != TCL_OK) {
      opserr << "WARNING invalid secTag\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);
    if (theSection == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTag;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &transfTag) != TCL_OK) {
      opserr << "WARNING invalid transfTag\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    int argi = 8;
    int numIter = 0;
    double tol = 0.0;
    if (argc > argi) {
      if (strcmp(argv[argi],"-iter") == 0) {
	if (argc < argi+3) {
	  opserr << "WARNING not enough -iter args need -iter numIter? tol?\n";
	  opserr << "forceBeamColumn element: " << eleTag << endln;
	  return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[argi+1], &numIter) != TCL_OK) {
	  opserr << "WARNING invalid numIter\n";
	  opserr << "forceBeamColumn element: " << eleTag << endln;
	  return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[argi+2], &tol) != TCL_OK) {
	  opserr << "WARNING invalid tol\n";
	  opserr << "forceBeamColumn element: " << eleTag << endln;
	  return TCL_ERROR;
	}
      }
    }
	
      
    if (ndm == 2) {
      
      theTransf2d = theTclBuilder->getCrdTransf2d(transfTag);
      
      if (theTransf2d == 0) {
	opserr << "WARNING transformation not found\n";
	opserr << "transformation: " << transfTag;
	opserr << "\nforceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
    }
    
    if (ndm == 3) {
      
      theTransf3d = theTclBuilder->getCrdTransf3d(transfTag);
      
      if (theTransf3d == 0) {
	opserr << "WARNING transformation not found\n";
	opserr << "transformation: " << transfTag;
	opserr << "\nforceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
    }

    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (int i = 0; i < nIP; i++)
      sections[i] = theSection;

    LobattoBeamIntegration beamIntegr;

    if (ndm == 2) {
      if (tol == 0.0)
	theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					   beamIntegr, *theTransf2d);
      else
	theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					   beamIntegr, *theTransf2d, 0.0, numIter, tol);
    }
    else {
      if (tol == 0.0)      
	theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					   beamIntegr, *theTransf3d);
      else
	theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					   beamIntegr, *theTransf3d, 0.0, numIter, tol);
    }

    delete [] sections;    
    if (theElement == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    if (theTclDomain->addElement(theElement) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      delete theElement;
      return TCL_ERROR;
    }

    return TCL_OK;
  } 

  
  //
  // otherwise use correct format of command as found in current documentation
  //

  if (Tcl_GetInt(interp, argv[5], &transfTag) != TCL_OK) {
    opserr << "WARNING invalid transfTag\n";
    opserr << "forceBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  if (ndm == 2) {
    
    theTransf2d = theTclBuilder->getCrdTransf2d(transfTag);
    
    if (theTransf2d == 0) {
      opserr << "WARNING transformation not found\n";
      opserr << "transformation: " << transfTag;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
  }
  
  if (ndm == 3) {
    
    theTransf3d = theTclBuilder->getCrdTransf3d(transfTag);
    
    if (theTransf3d == 0) {
      opserr << "WARNING transformation not found\n";
      opserr << "transformation: " << transfTag;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
  }
  
  if (strcmp(argv[6],"Lobatto") == 0 || strcmp(argv[6],"Legendre") == 0 
      || strcmp(argv[6],"Radau") == 0 || strcmp(argv[6],"NewtonCotes") == 0
      || strcmp(argv[6],"Trapezoidal") == 0) {
    int secTag, nIP;
    
    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? Lobatto secTag? nIP?\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &secTag) != TCL_OK) {
      opserr << "WARNING invalid secTag\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[8], &nIP) != TCL_OK) {
      opserr << "WARNING invalid nIP\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);
    if (theSection == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTag;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (int i = 0; i < nIP; i++)
      sections[i] = theSection;
    
    BeamIntegration *beamIntegr = 0;
    if (strcmp(argv[6],"Lobatto") == 0)
      beamIntegr = new LobattoBeamIntegration();
    else if (strcmp(argv[6],"Legendre") == 0)
      beamIntegr = new LegendreBeamIntegration();
    else if (strcmp(argv[6],"Radau") == 0)
      beamIntegr = new RadauBeamIntegration();
    else if (strcmp(argv[6],"NewtonCotes") == 0)
      beamIntegr = new NewtonCotesBeamIntegration();
    else if (strcmp(argv[6],"Trapezoidal") == 0)
      beamIntegr = new TrapezoidalBeamIntegration();
    else {
      opserr << "ERROR: invalid integration type: " << argv[6] << endln;
      return TCL_ERROR;
    }

    if (ndm == 2)
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					 *beamIntegr, *theTransf2d);
    else
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					 *beamIntegr, *theTransf3d);
    
    delete beamIntegr;
    delete [] sections;
  }

  else if (strcmp(argv[6],"UserDefined") == 0) {

    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? UserDefined nIP? secTag1? ... pt1? ... wt1? ...\n";
      return TCL_ERROR;
    }

    int nIP;
    
    if (Tcl_GetInt(interp, argv[7], &nIP) != TCL_OK) {
      opserr << "WARNING invalid nIP\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    ID secs(nIP);
    Vector pts(nIP);
    Vector wts(nIP);

    int i, j;
    for (i = 0, j = 8; i < nIP; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
	opserr << "WARNING invalid sec\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+nIP], &pt) != TCL_OK) {
	opserr << "WARNING invalid pt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+2*nIP], &wt) != TCL_OK) {
	opserr << "WARNING invalid wt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      secs(i) = sec;
      pts(i)  = pt;
      wts(i)  = wt;
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (i = 0; i < nIP; i++) {
      SectionForceDeformation *theSection = theTclBuilder->getSection(secs(i));
      if (theSection == 0) {
	opserr << "WARNING section not found\n";
	opserr << "Section: " << secs(i);
	opserr << "\nforceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      sections[i] = theSection;
    }
    
    UserDefinedBeamIntegration beamIntegr(nIP, pts, wts);

    if (ndm == 2)
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf2d);
    else
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf3d);
    
    delete [] sections;
  }

  else if (strcmp(argv[6],"HingeMidpoint") == 0 ||
	   strcmp(argv[6],"HingeRadau") == 0 ||
	   strcmp(argv[6],"HingeRadauTwo") == 0 ||
	   strcmp(argv[6],"HingeEndpoint") == 0) {
    
    if (argc < 12) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? type secTagI? lpI? secTagJ? lpJ? secTagE?\n";
      return TCL_ERROR;
    }

    int secTagI, secTagJ, secTagE;
    double lpI, lpJ;
    
    if (Tcl_GetInt(interp, argv[7], &secTagI) != TCL_OK) {
      opserr << "WARNING invalid secTagI\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &lpI) != TCL_OK) {
      opserr << "WARNING invalid lpI\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[9], &secTagJ) != TCL_OK) {
      opserr << "WARNING invalid secTagJ\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10], &lpJ) != TCL_OK) {
      opserr << "WARNING invalid lpJ\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[11], &secTagE) != TCL_OK) {
      opserr << "WARNING invalid secTagE\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    SectionForceDeformation *sectionI = theTclBuilder->getSection(secTagI);
    if (sectionI == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTagI;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    SectionForceDeformation *sectionJ = theTclBuilder->getSection(secTagJ);
    if (sectionJ == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTagJ;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    SectionForceDeformation *sectionE = theTclBuilder->getSection(secTagE);
    if (sectionJ == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTagE;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    BeamIntegration *beamIntegr = 0;

    SectionForceDeformation *sections[6];

    int Np;
    if (strcmp(argv[6],"HingeMidpoint") == 0) {
      beamIntegr = new HingeMidpointBeamIntegration(lpI, lpJ);
      Np = 4;
      sections[0] = sectionI;
      sections[1] = sectionE;
      sections[2] = sectionE;
      sections[3] = sectionJ;
    }
    else if (strcmp(argv[6],"HingeRadau") == 0) {
      beamIntegr = new HingeRadauBeamIntegration(lpI, lpJ);
      Np = 6;
      sections[0] = sectionI;
      sections[1] = sectionE;
      sections[2] = sectionE;
      sections[3] = sectionE;
      sections[4] = sectionE;
      sections[5] = sectionJ;
    }
    else if (strcmp(argv[6],"HingeRadauTwo") == 0) {
      beamIntegr = new HingeRadauTwoBeamIntegration(lpI, lpJ);
      Np = 6;
      sections[0] = sectionI;
      sections[1] = sectionI;
      sections[2] = sectionE;
      sections[3] = sectionE;
      sections[4] = sectionJ;
      sections[5] = sectionJ;
    }
    else {
      beamIntegr = new HingeEndpointBeamIntegration(lpI, lpJ);
      Np = 4;
      sections[0] = sectionI;
      sections[1] = sectionE;
      sections[2] = sectionE;
      sections[3] = sectionJ;
    }

    if (ndm == 2)
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, Np, sections,
					 *beamIntegr, *theTransf2d);
    else
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, Np, sections,
					 *beamIntegr, *theTransf3d);

    delete beamIntegr;
  }

  else if (strcmp(argv[6],"UserHinge") == 0) {

    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? UserHinge secTagE? npL? secTagL1? ... ptL1? ... wtL1? ... npR? secTagR1? ... ptR1? ... wtR1? ...\n";
      return TCL_ERROR;
    }

    int secTagE;
    
    if (Tcl_GetInt(interp, argv[7], &secTagE) != TCL_OK) {
      opserr << "WARNING invalid secTagE\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    int argStart = 8;

    int npL, npR;
      
    if (Tcl_GetInt(interp, argv[argStart], &npL) != TCL_OK) {
      opserr << "WARNING invalid npL\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[argStart+3*npL+1], &npR) != TCL_OK) {
      opserr << "WARNING invalid npR\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    int nIP = npL+npR;

    ID secs(nIP);
    Vector ptsL(npL);
    Vector wtsL(npL);
    Vector ptsR(npR);
    Vector wtsR(npR);

    int i, j;
    for (i = 0, j = argStart+1; i < npL; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
	opserr << "WARNING invalid sec\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+npL], &pt) != TCL_OK) {
	opserr << "WARNING invalid pt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+2*npL], &wt) != TCL_OK) {
	opserr << "WARNING invalid wt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      secs(i) = sec;
      ptsL(i) = pt;
      wtsL(i) = wt;
    }
    for (i = 0, j = 1+(argStart+1)+3*npL; i < npR; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
	opserr << "WARNING invalid sec\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+npR], &pt) != TCL_OK) {
	opserr << "WARNING invalid pt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+2*npR], &wt) != TCL_OK) {
	opserr << "WARNING invalid wt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      secs(i+npL) = sec;
      ptsR(i)     = pt;
      wtsR(i)     = wt;
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP+2];

    for (i = 0; i < nIP; i++) {
      SectionForceDeformation *theSection = theTclBuilder->getSection(secs(i));
      if (theSection == 0) {
	opserr << "WARNING section not found\n";
	opserr << "Section: " << secs(i);
	opserr << "\nforceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      sections[i] = theSection;
    }

    SectionForceDeformation *sectionE = theTclBuilder->getSection(secTagE);
    if (sectionE == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTagE;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    sections[nIP]   = sectionE;
    sections[nIP+1] = sectionE;
    
    UserDefinedHingeIntegration beamIntegr(npL, ptsL, wtsL,
					   npR, ptsR, wtsR);
		
    if (ndm == 2)
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP+2, sections,
					 beamIntegr, *theTransf2d);
    else
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP+2, sections,
					 beamIntegr, *theTransf3d);
    
    delete [] sections;
  }

  else if (strcmp(argv[6],"DistHinge") == 0) {
    
    if (argc < 16) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? type distType nIP? secTagI? lpI? secTagJ? lpJ? secTagE?\n";
      return TCL_ERROR;
    }

    int secTagI, secTagJ, secTagE;
    double lpI, lpJ;
    int nIP;

    BeamIntegration *otherBeamInt = 0;
    if (strcmp(argv[7],"Lobatto") == 0)
      otherBeamInt = new LobattoBeamIntegration();
    else if (strcmp(argv[7],"Legendre") == 0)
      otherBeamInt = new LegendreBeamIntegration();
    else if (strcmp(argv[7],"Radau") == 0)
      otherBeamInt = new RadauBeamIntegration();
    else if (strcmp(argv[7],"NewtonCotes") == 0)
      otherBeamInt = new NewtonCotesBeamIntegration();
    else if (strcmp(argv[7],"Trapezoidal") == 0)
      otherBeamInt = new TrapezoidalBeamIntegration();
    else {
      opserr << "ERROR: invalid integration type: " << argv[7] << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[8], &nIP) != TCL_OK) {
      opserr << "WARNING invalid nIP\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[9], &secTagI) != TCL_OK) {
      opserr << "WARNING invalid secTagI\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10], &lpI) != TCL_OK) {
      opserr << "WARNING invalid lpI\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[11], &secTagJ) != TCL_OK) {
      opserr << "WARNING invalid secTagJ\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[12], &lpJ) != TCL_OK) {
      opserr << "WARNING invalid lpJ\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[13], &secTagE) != TCL_OK) {
      opserr << "WARNING invalid secTagE\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    SectionForceDeformation *sectionI = theTclBuilder->getSection(secTagI);
    if (sectionI == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTagI;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    SectionForceDeformation *sectionJ = theTclBuilder->getSection(secTagJ);
    if (sectionJ == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTagJ;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    SectionForceDeformation *sectionE = theTclBuilder->getSection(secTagE);
    if (sectionJ == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTagE;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    int numSections = 2*nIP;
    SectionForceDeformation **sections = new SectionForceDeformation *[numSections];
    for (int i = 0; i < nIP; i++) {
      sections[i] = sectionI;
      sections[numSections-1-i] = sectionJ;
    }

    sections[numSections]   = sectionE;
    sections[numSections+1] = sectionE;
    
    DistHingeIntegration beamIntegr(lpI, lpJ, *otherBeamInt);

    if (ndm == 2)
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, numSections+2,
					 sections, beamIntegr, *theTransf2d);
    else
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, numSections+2,
					 sections, beamIntegr, *theTransf3d);

    if (otherBeamInt != 0)
      delete otherBeamInt;
    delete [] sections;
  }

  else if (strcmp(argv[6],"FixedLocation") == 0) {

    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? FixedLocation nIP? secTag1? ... pt1? ... \n";
      return TCL_ERROR;
    }

    int nIP;
    
    if (Tcl_GetInt(interp, argv[7], &nIP) != TCL_OK) {
      opserr << "WARNING invalid nIP\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    ID secs(nIP);
    Vector pts(nIP);
    Vector wts(nIP);

    int i, j;
    for (i = 0, j = 8; i < nIP; i++, j++) {
      int sec;
      double pt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
	opserr << "WARNING invalid sec\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+nIP], &pt) != TCL_OK) {
	opserr << "WARNING invalid pt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      secs(i) = sec;
      pts(i)  = pt;
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (i = 0; i < nIP; i++) {
      SectionForceDeformation *theSection = theTclBuilder->getSection(secs(i));
      if (theSection == 0) {
	opserr << "WARNING section not found\n";
	opserr << "Section: " << secs(i);
	opserr << "\nforceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      sections[i] = theSection;
    }
    
    FixedLocationBeamIntegration beamIntegr(nIP, pts);

    if (ndm == 2)
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf2d);
    else
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf3d);
    
    delete [] sections;
  }

  else if (strcmp(argv[6],"LowOrder") == 0) {

    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? LowOrder nIP? secTag1? ... pt1? ... wt1? ...\n";
      return TCL_ERROR;
    }

    int nIP;
    
    if (Tcl_GetInt(interp, argv[7], &nIP) != TCL_OK) {
      opserr << "WARNING invalid nIP\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    ID secs(nIP);
    Vector pts(nIP);
    Vector wts(nIP);

    int i, j;
    int nc = 0;
    for (i = 0, j = 8; i < nIP; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
	opserr << "WARNING invalid sec\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      secs(i) = sec;

      if (Tcl_GetDouble(interp, argv[j+nIP], &pt) != TCL_OK) {
	opserr << "WARNING invalid pt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      pts(i)  = pt;

      if (j+2*nIP < argc) {
	if (Tcl_GetDouble(interp, argv[j+2*nIP], &wt) != TCL_OK) {
	  opserr << "WARNING invalid wt\n";
	  opserr << "forceBeamColumn element: " << eleTag << endln;
	  return TCL_ERROR;
	}
	else {
	  wts(i)  = wt;
	  nc++;
	}
      }
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (i = 0; i < nIP; i++) {
      SectionForceDeformation *theSection = theTclBuilder->getSection(secs(i));
      if (theSection == 0) {
	opserr << "WARNING section not found\n";
	opserr << "Section: " << secs(i);
	opserr << "\nforceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      sections[i] = theSection;
    }
    
    LowOrderBeamIntegration beamIntegr(nIP, pts, nc, wts);

    if (ndm == 2)
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf2d);
    else
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf3d);
    
    delete [] sections;
  }

  else if (strcmp(argv[6],"MidDistance") == 0) {

    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? MidDistance nIP? secTag1? ... pt1? ... \n";
      return TCL_ERROR;
    }

    int nIP;
    
    if (Tcl_GetInt(interp, argv[7], &nIP) != TCL_OK) {
      opserr << "WARNING invalid nIP\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    ID secs(nIP);
    Vector pts(nIP);
    Vector wts(nIP);

    int i, j;
    for (i = 0, j = 8; i < nIP; i++, j++) {
      int sec;
      double pt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
	opserr << "WARNING invalid sec\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+nIP], &pt) != TCL_OK) {
	opserr << "WARNING invalid pt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      secs(i) = sec;
      pts(i)  = pt;
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (i = 0; i < nIP; i++) {
      SectionForceDeformation *theSection = theTclBuilder->getSection(secs(i));
      if (theSection == 0) {
	opserr << "WARNING section not found\n";
	opserr << "Section: " << secs(i);
	opserr << "\nforceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      sections[i] = theSection;
    }
    
    MidDistanceBeamIntegration beamIntegr(nIP, pts);

    if (ndm == 2)
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf2d);
    else
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf3d);
    
    delete [] sections;
  }

  else {
    opserr << "Unknown integration type: " << argv[6] << endln;
    opserr << "forceBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "forceBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  if (theTclDomain->addElement(theElement) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "forceBeamColumn element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }
  
  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}
