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
                                                                        
// $Revision: 1.11 $
// $Date: 2007-04-25 23:44:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/partitioned/PartitionedDomain.cpp,v $
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for PartitionedDomain.
// PartitionedDomain is an abstract class. The class is responsible for holding
// and providing access to the Elements, Nodes, LoadCases, SP_Constraints 
// and MP_Constraints just like a normal domain. In addition the domain provides
// a method to partition the domain into Subdomains.
//
// ModelBuilder. There are no partitions in a PartitionedDomain.
//
// What: "@(#) PartitionedDomain.C, revA"

#include <PartitionedDomain.h>
#include <stdlib.h>

#include <Matrix.h>

#include <DomainPartitioner.h>
#include <Element.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>
#include <ArrayOfTaggedObjects.h>
#include <ArrayOfTaggedObjectsIter.h>
#include <Subdomain.h>
#include <DomainPartitioner.h>
#include <PartitionedDomain.h>
#include <PartitionedDomainEleIter.h>
#include <PartitionedDomainSubIter.h>
#include <SingleDomEleIter.h>
#include <Vertex.h>
#include <Graph.h>
#include <LoadPattern.h>
#include <NodalLoad.h>
#include <ElementalLoad.h>
#include <SP_Constraint.h>
#include <Recorder.h>
#include <Parameter.h>

PartitionedDomain::PartitionedDomain()
:Domain(),
 theSubdomains(0),theDomainPartitioner(0),
 theSubdomainIter(0), mySubdomainGraph(0)
{
    elements = new ArrayOfTaggedObjects(1024);    
    theSubdomains = new ArrayOfTaggedObjects(32);
    theSubdomainIter = new PartitionedDomainSubIter(theSubdomains);

    mainEleIter = new SingleDomEleIter(elements);    
    theEleIter = new PartitionedDomainEleIter(this);
    
    if (theSubdomains == 0 || elements == 0 ||
	theSubdomainIter == 0 || 
	theEleIter == 0 || mainEleIter == 0) {
	
	opserr << "FATAL: PartitionedDomain::PartitionedDomain ";
	opserr << "  - ran out of memory\n";
	exit(-1);
    }
}


PartitionedDomain::PartitionedDomain(DomainPartitioner &thePartitioner)
:Domain(),
 theSubdomains(0),theDomainPartitioner(&thePartitioner),
 theSubdomainIter(0), mySubdomainGraph(0)
{
    elements = new ArrayOfTaggedObjects(1024);    
    theSubdomains = new ArrayOfTaggedObjects(32);
    theSubdomainIter = new PartitionedDomainSubIter(theSubdomains);

    mainEleIter = new SingleDomEleIter(elements);    
    theEleIter = new PartitionedDomainEleIter(this);
    
    if (theSubdomains == 0 || elements == 0 ||
	theSubdomainIter == 0 || theDomainPartitioner == 0 ||
	theEleIter == 0 || mainEleIter == 0) {
	
	opserr << "FATAL: PartitionedDomain::PartitionedDomain ";
	opserr << "  - ran out of memory\n";
	exit(-1);
    }
}


PartitionedDomain::PartitionedDomain(int numNodes, int numElements,
				     int numSPs, int numMPs, int numLoadPatterns,
				     int numSubdomains,
				     DomainPartitioner &thePartitioner)

:Domain(numNodes,0,numSPs,numMPs,numLoadPatterns),
 theSubdomains(0),theDomainPartitioner(&thePartitioner),
 theSubdomainIter(0), mySubdomainGraph(0)
{
    elements = new ArrayOfTaggedObjects(numElements);    
    theSubdomains = new ArrayOfTaggedObjects(numSubdomains);
    theSubdomainIter = new PartitionedDomainSubIter(theSubdomains);

    mainEleIter = new SingleDomEleIter(elements);    
    theEleIter = new PartitionedDomainEleIter(this);
    
    if (theSubdomains == 0 || elements == 0 ||
	theSubdomainIter == 0 || 
	theEleIter == 0 || mainEleIter == 0) {
	
	opserr << "FATAL: PartitionedDomain::PartitionedDomain(int ..) ";
	opserr << "  - ran out of memory\n";
	exit(-1);
    }
}




PartitionedDomain::~PartitionedDomain()
{
  this->clearAll();

  if (elements != 0)
    delete elements;
  
  if (theSubdomains != 0)
    delete theSubdomains;
  
  if (theSubdomainIter != 0)
    delete theSubdomainIter;
  
  if (theEleIter != 0)
    delete theEleIter;
}

void
PartitionedDomain::clearAll(void)
{
  this->Domain::clearAll();
  elements->clearAll();

  SubdomainIter &mySubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = mySubdomains()) != 0) 
    theSub->clearAll();

  theSubdomains->clearAll();
}
    


bool 
PartitionedDomain::addElement(Element *elePtr)
{
  if (elePtr->isSubdomain() == true)
    return this->addSubdomain((Subdomain *)elePtr);

  int eleTag = elePtr->getTag();
#ifdef _DEBUG      
    if (check == true) {

	// check ele Tag >= 0
	if (eleTag < 0) {
	    opserr << "PartitionedDomain::addElement - Element " << eleTag;
	    opserr << " tag must be >= 0\n";
	    return false;
	}      
	
	// check its not in this or any of the subdomains
	// MISSING CODE	
	
	// check all the elements nodes exist in the domain
	const ID &nodes = elePtr->getExternalNodes();
	for (int i=0; i<nodes.Size(); i++) {
	    int nodeTag = nodes(i);
	    Node *nodePtr = this->getNode(nodeTag);
	    if (nodePtr == 0) {
		opserr << "PartitionedDomain::addElement - In element " << eleTag;
		opserr << " no node " << nodeTag << " exists in the domain\n";
		return false;
	    }      	
	}
	
    }
#endif
    
    TaggedObject *other = elements->getComponentPtr(eleTag);
    if (other != 0)
	return false;
    
    bool result = elements->addComponent(elePtr);
    if (result == true) {
	elePtr->setDomain(this);
	elePtr->update();
	this->domainChange();
    }
    
    return result;
}    




bool 
PartitionedDomain::addNode(Node *nodePtr)
{
#ifdef _DEBUG    
    if (check == true) {
	// check its not in this or any of the subdomains

	// MISSING CODE	
    }
#endif
    return (this->Domain::addNode(nodePtr));    
}



bool
PartitionedDomain::addSP_Constraint(SP_Constraint *load)
{
  int nodeTag = load->getNodeTag();
  
  // check the Node exists in the Domain or one of Subdomains

  // if in Domain add it as external .. ignore Subdomains
  Node *nodePtr = this->getNode(nodeTag);
  if (nodePtr != 0) {
    return (this->Domain::addSP_Constraint(load));    
  }

  // find subdomain with node and add it .. break if find as internal node
  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = theSubdomains()) != 0) {
    bool res = theSub->hasNode(nodeTag);
    if (res == true) 
      return theSub->addSP_Constraint(load);
  }

  // if no subdomain .. node not in model .. error message and return failure
  opserr << "PartitionedDomain::addSP_Constraint - cannot add as node with tag" <<
    nodeTag << "does not exist in model\n"; 

  return false;
}



int
PartitionedDomain::addSP_Constraint(int startTag, int axisDirn, double axisValue, 
				   const ID &fixityCodes, double tol)
{
  int spTag = startTag;

  spTag = this->Domain::addSP_Constraint(spTag, axisDirn, axisValue, fixityCodes, tol);

  // find subdomain with node and add it .. break if find as internal node
  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = theSubdomains()) != 0) {
    spTag = theSub->addSP_Constraint(spTag, axisDirn, axisValue, fixityCodes, tol);
  }

  return spTag;
}


bool
PartitionedDomain::addSP_Constraint(SP_Constraint *load, int pattern)
{
  int nodeTag = load->getNodeTag();
  
  // check the Node exists in the Domain or one of Subdomains

  // if in Domain add it as external .. ignore Subdomains
  Node *nodePtr = this->getNode(nodeTag);
  if (nodePtr != 0) {
    return (this->Domain::addSP_Constraint(load, pattern));    
  }

  // find subdomain with node and add it .. break if find as internal node
  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = theSubdomains()) != 0) {
    bool res = theSub->hasNode(nodeTag);
    if (res == true) 
      return theSub->addSP_Constraint(load, pattern);
  }

    
  // if no subdomain .. node not in model .. error message and return failure
  opserr << "PartitionedDomain::addSP_Constraint - cannot add as node with tag" <<
    nodeTag << "does not exist in model\n"; 

  return false;
}


bool
PartitionedDomain::addMP_Constraint(MP_Constraint *load)
{
  bool res = true;
  bool getRetained = false;
  bool addedMain = false;

  // to every domain with the constrained node we must
  // 1. add the retained node if not already there & any sp constraints on that node
  // 2. add the constraint.

  int retainedNodeTag = load->getNodeRetained();
  int constrainedNodeTag = load->getNodeConstrained();

  //
  // first we check the main domain
  // if has the constrained but not retained we mark as needing retained
  //

  Node *retainedNodePtr = this->Domain::getNode(retainedNodeTag);
  Node *constrainedNodePtr = this->Domain::getNode(constrainedNodeTag);
  if (constrainedNodePtr != 0) {
    if (retainedNodePtr != 0) {

      res = this->Domain::addMP_Constraint(load);    
      if (res == false) {
	opserr << "PartitionedDomain::addMP_Constraint - problems adding to main domain\n";
	return res;
      }
      addedMain = true;
    } else {
      getRetained = true;
    }
  }

  //
  // now we check all subdomains
  // if a subdomain has both nodes we add the constraint, if only the
  // constrained node we mark as needing retained
  //

  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = theSubdomains()) != 0) {

    bool hasConstrained = theSub->hasNode(constrainedNodeTag);

    if (hasConstrained == true) {

      bool hasRestrained = theSub->hasNode(retainedNodeTag);

      if (hasRestrained == true) {

	res = theSub->addMP_Constraint(load);

	if (res == false) {
	  opserr << "PartitionedDomain::addMP_Constraint - problems adding to subdomain with retained\n";
	  return res;
	}
      } else
	getRetained = true;
    }
  }

  //
  // if getRetained is true, a subdomain or main domain has the constrained node
  // but no retained node .. we must go get it && SP constraints as well
  // 1. we first go get it
  // 2. then we add to main domain
  // 3. then we add to any subdomain
  // 

  if (getRetained == true) {

    // we get a copy of the retained Node, set mass equal 0 (don't want to double it)
    if (retainedNodePtr == 0) {
      SubdomainIter &theSubdomains2 = this->getSubdomains();
      while ((theSub = theSubdomains2()) != 0 && retainedNodePtr == 0) {
	
	bool hasRestrained = theSub->hasNode(retainedNodeTag);

	if (hasRestrained == true) {
	  retainedNodePtr = theSub->getNode(retainedNodeTag);

	  Matrix mass(retainedNodePtr->getNumberDOF(), retainedNodePtr->getNumberDOF());
	  mass.Zero();
	  retainedNodePtr->setMass(mass);
	}
      }

    } else {
      // get a copy & zero the mass
      retainedNodePtr = new Node(*retainedNodePtr, false);
    }

    if (retainedNodePtr == 0) {
      opserr << "partitionedDomain::addMP_Constraint - can't find retained node anywhere!\n";
      return false;
    }

    //
    // if main has it we add the retained to main & constraint
    //

    if (constrainedNodePtr != 0 && addedMain == false) {
      res = this->Domain::addNode(retainedNodePtr);

      if (res == false) {
	opserr << "PartitionedDomain::addMP_Constraint - problems adding retained to main domain\n";
	return res;
      } 
      res = this->Domain::addMP_Constraint(load);

      if (res == false) {
	opserr << "PartitionedDomain::addMP_Constraint - problems adding constraint to main domain after adding node\n";
	return res;
      } 
    }

    //
    // to subdmains that have the constrained but no retained
    // 1. we add a copy of retained
    // 2. we add the constraint
    //

    SubdomainIter &theSubdomains3 = this->getSubdomains();
    while ((theSub = theSubdomains3()) != 0 && retainedNodePtr == 0) {
      bool hasConstrained = theSub->hasNode(constrainedNodeTag);
      if (hasConstrained == true) {
	bool hasRestrained = theSub->hasNode(retainedNodeTag);
	if (hasRestrained == false) {
	  res = theSub->addNode(retainedNodePtr);

	  if (res == false) {
	    opserr << "PartitionedDomain::addMP_Constraint - problems adding retained to subdomain\n";
	    return res;
	  } 
	  res = theSub->addMP_Constraint(load);

	  if (res == false) {
	    opserr << "PartitionedDomain::addMP_Constraint - problems adding constraint to subdomain after adding node\n";
	    return res;
	  } 
	}
      }
    }

    // clean up memory

    if (constrainedNodePtr != 0 && addedMain == false) 
      ; 
    else
      delete retainedNodePtr;
    
  }

  return res;

}


bool 
PartitionedDomain::addLoadPattern(LoadPattern *loadPattern)
{
  bool result = true;

  int tag = loadPattern->getTag();
  if (this->getLoadPattern(tag) != 0) {
    opserr << "PartitionedDomain::addLoadPattern - cannot add as LoadPattern with tag" <<
      tag << "already exists in model\n";             
    return false;
  }

  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = theSubdomains()) != 0) {
    bool res = theSub->addLoadPattern(loadPattern);
    if (res != true) {
      opserr << "PartitionedDomain::addLoadPattern - cannot add as LoadPattern with tag: " <<
	tag << " to subdomain\n";             
      result = res;
    }
  }

  this->Domain::addLoadPattern(loadPattern);

  return result;
}    


bool 
PartitionedDomain::addNodalLoad(NodalLoad *load, int pattern)
{
  int nodeTag = load->getNodeTag();
  
  // check the Node exists in the Domain or one of Subdomains

  // if in Domain add it as external .. ignore Subdomains
  Node *nodePtr = this->getNode(nodeTag);
  if (nodePtr != 0) {
    return (this->Domain::addNodalLoad(load, pattern));    
  }


  // find subdomain with node and add it .. break if find as internal node
  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = theSubdomains()) != 0) {
    bool res = theSub->hasNode(nodeTag);
    if (res == true) {
      // opserr << "PartitionedDomain::addLoadPattern(LoadPattern *loadPattern) SUB " << theSub->getTag() << *load;
      return theSub->addNodalLoad(load, pattern);
    }
  }

  // if no subdomain .. node not in model
  opserr << "PartitionedDomain::addNodalLoad - cannot add as node with tag" <<
    nodeTag << "does not exist in model\n"; 
  return false;
}    


bool 
PartitionedDomain::addElementalLoad(ElementalLoad *load, int pattern)
{
  opserr << "PartitionedDomain::addElementalLoad - not yet implemented\n";
  return false;
}


Element *
PartitionedDomain::removeElement(int tag)
{
    // we first see if its in the original domain
    TaggedObject *res = elements->removeComponent(tag);
    Element *result = 0;
    if (res != 0) {
	result = (Element *)res;
	this->domainChange();
	return result;
    }

    // if not there we must check all the other subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    result = theSub->removeElement(tag);
	    if (result != 0) {
		return result;
	    }
	}
    }

    // its not there
    return 0;
}    


Node *
PartitionedDomain::removeNode(int tag)    
{
    // we first remove it form the original domain (in case on boundary)
    Node *result = this->Domain::removeNode(tag);

    // we must also try removing from the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    Node *res = theSub->removeNode(tag);
	    if (res != 0) 
		result = res;
	}
    }
    
    if (result != 0) 
	this->domainChange();
    
    return result;
}    

SP_Constraint *
PartitionedDomain::removeSP_Constraint(int tag)
{
    // we first see if its in the original domain
    SP_Constraint *result = this->Domain::removeSP_Constraint(tag);
    if (result != 0) {
	this->domainChange();
	return result;
    }
	

    // if not there we must check all the other subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    result = theSub->removeSP_Constraint(tag);
	    if (result != 0) {
		return result;
	    }
	}
    }

    // its not there
    return 0;
}


MP_Constraint *
PartitionedDomain::removeMP_Constraint(int tag)
{
    // we first see if its in the original domain
    MP_Constraint *result = this->Domain::removeMP_Constraint(tag);
    if (result != 0) {
	this->domainChange();
	return result;
    }

    // if not there we must check all the other subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);		
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    result = theSub->removeMP_Constraint(tag);
	    if (result != 0) {
		return result;
	    }
	}
    }    

    // its not there
    return 0;
}


LoadPattern * 
PartitionedDomain::removeLoadPattern(int tag)
{
    // we first see if its in the original domain
    LoadPattern *result = this->Domain::removeLoadPattern(tag);

    // we must also try removing from the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    LoadPattern *res = theSub->removeLoadPattern(tag);
	    if (res != 0) 
		result = res;
	}
    }
    
    if (result != 0) 
	this->domainChange();
    
    return result;
}    

// public member functions which have to be modified
ElementIter       &
PartitionedDomain::getElements()
{
    theEleIter->reset();
    return *theEleIter;
}    


Element  *
PartitionedDomain::getElement(int tag) 
{
    // we first see if its in the original domain
    TaggedObject *res = elements->getComponentPtr(tag);
    Element *result =0;
    if (res != 0) {
	result = (Element *)res;
	return result;
    }

    /*
    // go through the other subdomains until we find it or we run out of subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    result = theSub->getElement(tag);
	    if (result != 0) 
		return result;
	}
    }
    */

    // its not there
    return 0;
}


int 		
PartitionedDomain::getNumElements(void) const
{
    int result = elements->getNumComponents();

    // add the number of subdomains
    result +=  theSubdomains->getNumComponents();
    return result;
}

void
PartitionedDomain::applyLoad(double timeStep)
{
    this->Domain::applyLoad(timeStep);

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->applyLoad(timeStep);
	}
    }
}


void
PartitionedDomain::setCommitTag(int newTag)
{
    this->Domain::setCommitTag(newTag);

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->setCommitTag(newTag);
	}
    }
}



void
PartitionedDomain::setCurrentTime(double newTime)
{
    this->Domain::setCurrentTime(newTime);

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->setCurrentTime(newTime);
	}
    }
}


void
PartitionedDomain::setCommittedTime(double newTime)
{
    this->Domain::setCommittedTime(newTime);

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->setCommittedTime(newTime);
	}
    }
}


void
PartitionedDomain::setLoadConstant(void)
{
    this->Domain::setLoadConstant();

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->setLoadConstant();
	}
    }
}


int
PartitionedDomain::setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc)
{
    this->Domain::setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
	}
    }
    return 0;
}


int
PartitionedDomain::update(void)
{
  int res = this->Domain::update();

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      theSub->computeNodalResponse();
      theSub->update();
    }
  }

#ifdef _PARALLEL_PROCESSING
  return this->barrierCheck(res);
#endif
  return 0;
}




#ifdef _PARALLEL_PROCESSING
int
PartitionedDomain::barrierCheck(int res)
{
  int result = res;

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      int subResult = theSub->barrierCheckIN();
      if (subResult != 0)
	result = subResult;
    }

    ArrayOfTaggedObjectsIter theSubsIter1(*theSubdomains);	
    while ((theObject = theSubsIter1()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      theSub->barrierCheckOUT(result);
    }
  }

  return result;
}
#endif

int
PartitionedDomain::update(double newTime, double dT)
{
  this->applyLoad(newTime);
  int res = this->Domain::update();

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      theSub->computeNodalResponse();
      theSub->update(newTime, dT);
    }
  }

#ifdef _PARALLEL_PROCESSING
  return this->barrierCheck(res);
#endif
  return 0;

  /*

  opserr << "PartitionedDomain::update(double newTime, double dT) -1\n";
  int result = 0;


  opserr << "PartitionedDomain::update(double newTime, double dT) -2\n";
  this->update();
  opserr << "PartitionedDomain::update(double newTime, double dT) -2a\n";

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      theSub->update(newTime, dT);
    }
    this->barrierCheck(result);
  }
  opserr << "PartitionedDomain::update(double newTime, double dT) -3\n";
  return result;

*/

}


int
PartitionedDomain::hasDomainChanged(void)
{
  return this->Domain::hasDomainChanged();
}

int
PartitionedDomain::newStep(double dT)
{
  // first we need to see if any subdomain has changed & mark the change in domain
  bool domainChangedAnySubdomain = this->getDomainChangeFlag();
  if (domainChangedAnySubdomain == false) {
    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while (((theObject = theSubsIter()) != 0) && (domainChangedAnySubdomain == false)) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    domainChangedAnySubdomain = theSub->getDomainChangeFlag();
	}
    }
  }

  if (domainChangedAnySubdomain == true) {
    this->Domain::domainChange();
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while (((theObject = theSubsIter()) != 0) && (domainChangedAnySubdomain == false)) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->domainChange();
	}
    }
  }

  this->Domain::newStep(dT);

    int res = 0;
    // do the same for all the subdomains
    if (theSubdomains != 0) {
      ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
      TaggedObject *theObject;
      while ((theObject = theSubsIter()) != 0) {
	Subdomain *theSub = (Subdomain *)theObject;	    
	res += theSub->newStep(dT);
	if (res != 0) 
	  opserr << "PartitionedDomain::step - subdomain " << theSub->getTag() << " failed in step\n";
      }
    }
    return res;
}




int
PartitionedDomain::commit(void)
{
  int result = this->Domain::commit();
  if (result < 0) {
    opserr << "PartitionedDomain::commit(void) - failed in Domain::commit()\n";
    return result;
  }

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      int res = theSub->commit();
      if (res < 0) {
	opserr << "PartitionedDomain::commit(void)";
	opserr << " - failed in Subdomain::commit()\n";
	return res;
      }	    
    }
  }

  // now we load balance if we have subdomains and a partitioner
  int numSubdomains = this->getNumSubdomains();
  if (numSubdomains != 0 && theDomainPartitioner != 0)  {
    Graph &theSubGraphs = this->getSubdomainGraph();
    theDomainPartitioner->balance(theSubGraphs);
  }

  return 0;
}


int
PartitionedDomain::revertToLastCommit(void)
{
    int result = this->Domain::revertToLastCommit();
    if (result < 0) {
	opserr << "PartitionedDomain::revertToLastCommit(void) - failed in Domain::revertToLastCommit()\n";
	return result;
    }

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    int res = theSub->revertToLastCommit();
	    if (res < 0) {
		opserr << "PartitionedDomain::revertToLastCommit(void)";
		opserr << " - failed in Subdomain::revertToLastCommit()\n";
		return res;
	    }	    
	}
    }

    return 0;
}

int
PartitionedDomain::revertToStart(void)
{
    int result = this->Domain::revertToStart();
    if (result < 0) {
	opserr << "PartitionedDomain::revertToLastCommit(void) - failed in Domain::revertToLastCommit()\n";
	return result;
    }

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    int res = theSub->revertToStart();
	    if (res < 0) {
		opserr << "PartitionedDomain::revertToLastCommit(void)";
		opserr << " - failed in Subdomain::revertToLastCommit()\n";
		return res;
	    }	    
	}
    }

    return 0;
}


int  
PartitionedDomain::addRecorder(Recorder &theRecorder)
{
  if (this->Domain::addRecorder(theRecorder) < 0)
    return -1;

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      int res = theSub->addRecorder(theRecorder);
      if (res < 0) {
	opserr << "PartitionedDomain::revertToLastCommit(void)";
	opserr << " - failed in Subdomain::revertToLastCommit()\n";
	return res;
      }	    
    }
  }
  return 0;
}

int  
PartitionedDomain::removeRecorders(void)
{
  if (this->Domain::removeRecorders() < 0)
    return -1;

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      int res = theSub->removeRecorders();
      if (res < 0) {
	opserr << "PartitionedDomain::revertToLastCommit(void)";
	opserr << " - failed in Subdomain::revertToLastCommit()\n";
	return res;
      }	    
    }
  }
  return 0;
}

void 
PartitionedDomain::Print(OPS_Stream &s, int flag)
{
  this->Domain::Print(s, flag);

  s << "\nELEMENT DATA: NumEle: " << elements->getNumComponents() << "\n";
  elements->Print(s);
	
  // print all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      theObject->Print(s, flag);
    }
  }
}


void 
PartitionedDomain::Print(OPS_Stream &s, ID *nodeTags, ID *eleTags, int flag)
{
  if (nodeTags != 0)
    this->Domain::Print(s, nodeTags, 0, flag);

  if (eleTags != 0) {
    int numEle = eleTags->Size();
    for (int i=0; i<numEle; i++) {
      int eleTag = (*eleTags)(i);
      TaggedObject *theEle = elements->getComponentPtr(eleTag);
      if (theEle != 0)
	theEle->Print(s, flag);
    }
  }
  

  // print all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      theSub->Print(s, nodeTags, eleTags, flag);
    }
  }
}



int 
PartitionedDomain::setPartitioner(DomainPartitioner *thePartitioner)
{
  theDomainPartitioner = thePartitioner;
  return 0;
}


int 
PartitionedDomain::partition(int numPartitions, bool usingMain, int mainPartitionID)
{
  int result = 0;
    // need to create element graph before create new subdomains
    // DO NOT REMOVE THIS LINE __ EVEN IF COMPILER WARNING ABOUT UNUSED VARIABLE
    Graph &theEleGraph = this->getElementGraph();
    
    // now we call partition on the domainPartitioner which does the partitioning
    DomainPartitioner *thePartitioner = this->getPartitioner();
    if (thePartitioner != 0) {
      thePartitioner->setPartitionedDomain(*this);
      result =  thePartitioner->partition(numPartitions, usingMain, mainPartitionID);
    } else {
      opserr << "PartitionedDomain::partition(int numPartitions) - no associated partitioner\n";
      return -1;
    }

    //
    // add recorder objects
    //

    // do the same for all the subdomains
    if (theSubdomains != 0) {
      ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
      TaggedObject *theObject;
      while ((theObject = theSubsIter()) != 0) {
	Subdomain *theSub = (Subdomain *)theObject;	    
	for (int i=0; i<numRecorders; i++) {
	  int res = theSub->addRecorder(*theRecorders[i]);
	  if (res != 0) {
	    opserr << "PartitionedDomain::revertToLastCommit(void)";
	    opserr << " - failed in Subdomain::revertToLastCommit()\n";
	    return res;
	  }	  
	}  
      }
    }

    
    return result;
}


bool 
PartitionedDomain::addSubdomain(Subdomain *theSubdomain)
{
    int eleTag = theSubdomain->getTag();
    TaggedObject *other = theSubdomains->getComponentPtr(eleTag);
    if (other != 0)
	return false;
  
    bool result = theSubdomains->addComponent(theSubdomain);
    if (result == true) {
	theSubdomain->setDomain(this);
	this->domainChange();
    }
  
  return result;

}

int 
PartitionedDomain::getNumSubdomains(void)
{
    return theSubdomains->getNumComponents();
}
    
Subdomain *
PartitionedDomain::getSubdomainPtr(int tag)
{
    TaggedObject *mc = theSubdomains->getComponentPtr(tag);
    if (mc == 0) return 0;
    Subdomain *result = (Subdomain *)mc;
    return result;
}    

SubdomainIter &
PartitionedDomain::getSubdomains(void)
{
    theSubdomainIter->reset();
    return *theSubdomainIter;
}



DomainPartitioner *
PartitionedDomain::getPartitioner(void) const
{
    return theDomainPartitioner;
}
	



int 
PartitionedDomain::buildEleGraph(Graph *theEleGraph)
{
    int numVertex = elements->getNumComponents();

    // see if quick return

    if (numVertex == 0) 
	return 0;
    
    // create another vertices array which aids in adding edges
    
    int *theElementTagVertices = 0;
    int maxEleNum = 0;
    
    TaggedObject *tagdObjPtr;
    TaggedObjectIter &theEles = elements->getComponents();
    while ((tagdObjPtr = theEles()) != 0)
	if (tagdObjPtr->getTag() > maxEleNum)
	    maxEleNum = tagdObjPtr->getTag();

    theElementTagVertices = new int[maxEleNum+1];

    if (theElementTagVertices == 0) {
	opserr << "WARNING Domain::buildEleGraph ";
	opserr << " - Not Enough Memory for ElementTagVertices\n";
	return -1;
    }

    for (int j=0; j<=maxEleNum; j++) theElementTagVertices[j] = -1;

    // now create the vertices with a reference equal to the element number.
    // and a tag which ranges from 0 through numVertex-1
    
    TaggedObjectIter &theEles2 = elements->getComponents();
    
    int count = START_VERTEX_NUM;
    while ((tagdObjPtr = theEles2()) != 0) {
	int ElementTag = tagdObjPtr->getTag();
	Vertex *vertexPtr = new Vertex(count,ElementTag);

	if (vertexPtr == 0) {
	    opserr << "WARNING Domain::buildEleGraph";
	    opserr << " - Not Enough Memory to create ";
	    opserr << count << "th Vertex\n";
	    delete [] theElementTagVertices;
	    return -1;
	}

	theEleGraph->addVertex(vertexPtr);
	theElementTagVertices[ElementTag] = count++;
	
    }

    // We now need to determine which elements are asssociated with each node.
    // As this info is not in the Node interface we must build it; which we
    // do using vertices for each node, when we addVertex at thes nodes we
    // will not be adding vertices but element tags.

    Vertex **theNodeTagVertices = 0;
    int maxNodNum = 0;
    Node *nodPtr;
    NodeIter &nodeIter = this->getNodes();
    while ((nodPtr = nodeIter()) != 0)
	if (nodPtr->getTag() > maxNodNum)
	    maxNodNum = nodPtr->getTag();

    theNodeTagVertices = new Vertex *[maxNodNum+1];

    if (theNodeTagVertices == 0) {
	opserr << "WARNING Domain::buildEleGraph ";
	opserr << " - Not Enough Memory for NodeTagVertices\n";
	return -1;
    }

    for (int l=0; l<=maxNodNum; l++) theNodeTagVertices[l] = 0;

    // now create the vertices with a reference equal to the node number.
    // and a tag which ranges from 0 through numVertex-1 and placed in
    // theNodeTagVertices at a position equal to the node's tag.

    NodeIter &nodeIter2 = this->getNodes();
    count = START_VERTEX_NUM;
    while ((nodPtr = nodeIter2()) != 0) {
	int nodeTag = nodPtr->getTag();
	Vertex *vertexPtr = new Vertex(count++,nodeTag);
	theNodeTagVertices[nodeTag] = vertexPtr;

	if (vertexPtr == 0) {
	    opserr << "WARNING Domain::buildEleGraph";
	    opserr << " - Not Enough Memory to create ";
	    opserr << count << "th Node Vertex\n";
	    delete [] theNodeTagVertices;
	    return -1;
	}
    }

    // now add the the Elements to the nodes
    Element *elePtr;
    TaggedObjectIter &theEles3 = elements->getComponents();
    
    while((tagdObjPtr = theEles3()) != 0) {
	elePtr = (Element *)tagdObjPtr;
	int eleTag = elePtr->getTag();
	const ID &id = elePtr->getExternalNodes();

	int size = id.Size();
	for (int i=0; i<size; i++) 
	    theNodeTagVertices[id(i)]->addEdge(eleTag);
    }

    // now add the edges to the vertices of our element graph;
    // this is done by looping over the Node vertices, getting their 
    // Adjacenecy and adding edges between elements with common nodes


    Vertex *vertexPtr;
    for (int k=0; k<=maxNodNum; k++)
	if ((vertexPtr = theNodeTagVertices[k]) != 0) {

	    const ID &id = vertexPtr->getAdjacency();

	    int size = id.Size();
	    for (int i=0; i<size; i++) {
		int Element1 = id(i);

		int vertexTag1 = theElementTagVertices[Element1];

		for (int j=0; j<size; j++) 
		    if (i != j) {

			int Element2 = id(j);
			int vertexTag2 = theElementTagVertices[Element2];

			// addEdge() adds for both vertices - do only once
			if (vertexTag1 > vertexTag2) 
			    theEleGraph->addEdge(vertexTag1,vertexTag2);
			    theEleGraph->addEdge(vertexTag2,vertexTag1);			
		    }
	    }
	}

    // done now delete theElementTagVertices, the node Vertices and
    // theNodeTagVertices
   
    delete [] theElementTagVertices;    
    
    for (int i=0; i<=maxNodNum; i++)
	if ((vertexPtr = theNodeTagVertices[i]) != 0) 
	    delete vertexPtr;
	    
    delete [] theNodeTagVertices;

    return 0;
    
}



// a method which will only remove a node from the partitioned domain
// it does not touch the subdomains .. can be dangerous to use.
Node *
PartitionedDomain::removeExternalNode(int tag)
{
    return (this->Domain::removeNode(tag));        
}

Graph &
PartitionedDomain::getSubdomainGraph(void)
{
    // delete the old always - only object that will 
    // use this is a DomainBalancer & it is always looking for latest
    if (mySubdomainGraph != 0) {
	delete mySubdomainGraph;
	mySubdomainGraph = 0;
    }

    // create a new graph
    if (mySubdomainGraph == 0)
        mySubdomainGraph = new Graph(this->getNumSubdomains()+START_VERTEX_NUM);

    if (mySubdomainGraph == 0) // if still 0 try a smaller one
        mySubdomainGraph = new Graph();    

    int numVertex = theSubdomains->getNumComponents();

    // see if quick return

    if (numVertex == 0) 
	return *mySubdomainGraph;
    
    // create another vertices array which aids in adding edges
    
    int *theElementTagVertices = 0;
    int maxEleNum = 0;

    TaggedObject *tagdObjPtr;
    TaggedObjectIter &theEles = theSubdomains->getComponents();
    while ((tagdObjPtr = theEles()) != 0)
	if (tagdObjPtr->getTag() > maxEleNum)
	    maxEleNum = tagdObjPtr->getTag();

    theElementTagVertices = new int[maxEleNum+1];

    if (theElementTagVertices == 0) {
	opserr << "WARNING PartitionedDomain::buildEleGraph ";
	opserr << " - Not Enough Memory for ElementTagVertices\n";
	exit(-1);
    }

    for (int j=0; j<=maxEleNum; j++) theElementTagVertices[j] = -1;

    // now create the vertices with a reference equal to the subdomain number.
    // and a tag equal to the subdomain number and a weighed according to 
    // the subdomain cost 
    
    TaggedObjectIter &theEles2 = theSubdomains->getComponents();

    while ((tagdObjPtr = theEles2()) != 0) {
	Subdomain *theSub = (Subdomain *)tagdObjPtr; // upward cast ok as
	                                     // only subdomais can be added
	int ElementTag = tagdObjPtr->getTag();

	Vertex *vertexPtr = new Vertex(ElementTag, ElementTag, theSub->getCost()); 
	if (vertexPtr == 0) {
	    opserr << "WARNING Domain::buildEleGraph";
	    opserr << " - Not Enough Memory to create ";
	    opserr << ElementTag << "th Vertex\n";
	    delete [] theElementTagVertices;
	    exit(-1);
	}

	mySubdomainGraph->addVertex(vertexPtr);

	theElementTagVertices[ElementTag] = ElementTag;
    }

    // We now need to determine which theSubdomains are asssociated with each node.
    // As this info is not in the Node interface we must build it; which we
    // do using vertices for each node, when we addVertex at thes nodes we
    // will not be adding vertices but element tags.

    Vertex **theNodeTagVertices = 0;
    int maxNodNum = 0;
    Node *nodPtr;
    NodeIter &nodeIter = this->getNodes();
    while ((nodPtr = nodeIter()) != 0)
	if (nodPtr->getTag() > maxNodNum)
	    maxNodNum = nodPtr->getTag();

    theNodeTagVertices = new Vertex *[maxNodNum+1];

    if (theNodeTagVertices == 0) {
	opserr << "WARNING Domain::buildEleGraph ";
	opserr << " - Not Enough Memory for NodeTagVertices\n";
	exit(-1);
    }

    for (int l=0; l<=maxNodNum; l++) theNodeTagVertices[l] = 0;

    // now create the vertices with a reference equal to the node number.
    // and a tag which ranges from 0 through numVertex-1 and placed in
    // theNodeTagVertices at a position equal to the node's tag.

    NodeIter &nodeIter2 = this->getNodes();
    int count = START_VERTEX_NUM;
    while ((nodPtr = nodeIter2()) != 0) {
	int nodeTag = nodPtr->getTag();
	Vertex *vertexPtr = new Vertex(count++,nodeTag);
	theNodeTagVertices[nodeTag] = vertexPtr;

	if (vertexPtr == 0) {
	    opserr << "WARNING Domain::buildEleGraph";
	    opserr << " - Not Enough Memory to create "; opserr << count << "th Node Vertex\n";
	    delete [] theNodeTagVertices;
	    exit(-1);
	}
    }

    // now add the the TheSubdomains to the nodes
    Element *elePtr;
    TaggedObjectIter &theEles3 = theSubdomains->getComponents();
    
    while((tagdObjPtr = theEles3()) != 0) {
	elePtr = (Element *)tagdObjPtr;
	int eleTag = elePtr->getTag();
	const ID &id = elePtr->getExternalNodes();

	int size = id.Size();
	for (int i=0; i<size; i++) 
	    theNodeTagVertices[id(i)]->addEdge(eleTag);
    }

    // now add the edges to the vertices of our element graph;
    // this is done by looping over the Node vertices, getting their 
    // Adjacenecy and adding edges between theSubdomains with common nodes

    Vertex *vertexPtr;
    for (int k=0; k<=maxNodNum; k++)
	if ((vertexPtr = theNodeTagVertices[k]) != 0) {

	    const ID &id = vertexPtr->getAdjacency();

	    int size = id.Size();
	    for (int i=0; i<size; i++) {
		int Element1 = id(i);

		int vertexTag1 = theElementTagVertices[Element1];

		for (int j=0; j<size; j++) 
		    if (i != j) {

			int Element2 = id(j);
			int vertexTag2 = theElementTagVertices[Element2];

			// addEdge() adds for both vertices - do only once
			if (vertexTag1 > vertexTag2) 
			    mySubdomainGraph->addEdge(vertexTag1,vertexTag2);
			    mySubdomainGraph->addEdge(vertexTag2,vertexTag1);			
		    }
	    }
	}

    // done now delete theElementTagVertices, the node Vertices and
    // theNodeTagVertices
   
    delete [] theElementTagVertices;    

    for (int i=0; i<=maxNodNum; i++)
	if ((vertexPtr = theNodeTagVertices[i]) != 0) 
	    delete vertexPtr;

    delete [] theNodeTagVertices;

    return *mySubdomainGraph;
}


double
PartitionedDomain::getNodeDisp(int nodeTag, int dof, int &errorFlag)
{
  double result = this->Domain::getNodeDisp(nodeTag, dof, errorFlag);

  if (errorFlag != 0) {

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0 && errorFlag != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    result = theSub->getNodeDisp(nodeTag, dof, errorFlag);
	    if (errorFlag == 0)
	      return result;
	}	    
    }
  }
  
  return result;
}


int
PartitionedDomain::setMass(const Matrix &mass, int nodeTag)
{
  int result = this->Domain::setMass(mass, nodeTag);

  if (result != 0) {

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0 && result != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    result = theSub->setMass(mass, nodeTag);
	}	    
    }
  }
  
  return result;
}

const Vector *
PartitionedDomain::getNodeResponse(int nodeTag, NodeResponseType response)
{
  const Vector *res = this->Domain::getNodeResponse(nodeTag, response); 
  if (res != 0)
    return res;

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      const Vector *result = theSub->getNodeResponse(nodeTag, response); 
      if (result != 0)
	return result;
    }	    
  }

  return NULL;
}

int 
PartitionedDomain::calculateNodalReactions(bool inclInertia)
{
  int res = this->Domain::calculateNodalReactions(inclInertia); 

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      res += theSub->calculateNodalReactions(inclInertia); 
    }
  }
  return res;
}

bool 
PartitionedDomain::addParameter(Parameter *param)
{
    bool res = this->Domain::addParameter(param);

    // do the same for all the subdomains
    if (theSubdomains != 0) {
      ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
      TaggedObject *theObject;
      while ((theObject = theSubsIter()) != 0) {
	Subdomain *theSub = (Subdomain *)theObject;	    
	theSub->addParameter(param);
      }
    }

    return res;
}

Parameter *
PartitionedDomain::removeParameter(int tag)
{
  Parameter *res = this->Domain::removeParameter(tag);

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      theSub->removeParameter(tag);
    }
  }

  return res;
}


int 
PartitionedDomain::updateParameter(int tag, int value)
{
  int res = this->Domain::updateParameter(tag, value);

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {

      Subdomain *theSub = (Subdomain *)theObject;	    
      theSub->updateParameter(tag, value);
    }
  }

  return res;
}


int 
PartitionedDomain::updateParameter(int tag, double value)
{
  int res = this->Domain::updateParameter(tag, value);

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;	    
      theSub->updateParameter(tag, value);
    }
  }

  return res;
}

