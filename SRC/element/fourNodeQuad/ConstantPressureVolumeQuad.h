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
                                                                        
// $Revision: 1.7 $
// $Date: 2002-06-07 00:29:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/ConstantPressureVolumeQuad.h,v $

// Ed "C++" Love
//
// Constant Presssure/Volume Four Node Quadrilateral
// Plane Strain (NOT PLANE STRESS)

#include <iostream.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>

class ConstantPressureVolumeQuad : public Element 
{

  public :
    
    //null constructor
    ConstantPressureVolumeQuad( ) ;
  
    //full constructor
    ConstantPressureVolumeQuad( int tag, 
			        int node1,
			        int node2,
			        int node3,
			        int node4,
			        NDMaterial &theMaterial ) ;

    //destructor 
    virtual ~ConstantPressureVolumeQuad( ) ;

    //set domain
    void setDomain( Domain *theDomain ) ;

    //get the number of external nodes
    int getNumExternalNodes( ) const ;
 
    //return connected external nodes
    const ID &getExternalNodes( ) ;

    //return number of dofs
    int getNumDOF( ) ;

    //commit state
    int commitState( ) ;
    
    //revert to last commit 
    int revertToLastCommit( ) ;
    
    //revert to start 
    int revertToStart( ) ;

    //print out element data
    void Print( ostream &s, int flag ) ;
	
    //return stiffness matrix 
    const Matrix &getTangentStiff( ) ;
    
    //return secant matrix 
    const Matrix &getSecantStiff( ) ;
    
    //return damping matrix
    const Matrix &getDamp( ) ;
    
    //return mass matrix
    const Matrix &getMass( ) ;

    void zeroLoad( ) ;
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    //get residual
    const Vector &getResistingForce( ) ;
    
    //get residual with inertia terms
    const Vector &getResistingForceIncInertia( ) ;

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);

    //plotting 
    int displaySelf(Renderer &theViewer, int displayMode, float fact);
  
  private : 

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;
    
    //volume-pressure constants
    static double one3 ;
    static double two3 ;
    static double four3 ;
    static double one9 ;
    
    //quadrature data
    static double root3 ;
    static double one_over_root3 ;    
    static double sg[4] ;
    static double tg[4] ;
    static double wg[4] ;

    
    //node information
    ID connectedExternalNodes ;  //four node numbers
    Node *nodePointers[4] ;      //pointers to four nodes
					
    //material information
    NDMaterial *materialPointers[4] ; //pointers to four materials
					  
    //nodal coordinates, two coordinates for each of four nodes
    double xl[2][4] ; 
    
    //form residual and tangent					  
    void formResidAndTangent( int tang_flag ) ;

    //inertia terms
    void formInertiaTerms( int tangFlag ) ;

    Matrix transpose( int dim1, int dim2, const Matrix &M ) ;
   
    //shape function routine for four node quads
    void shape2d( double ss, double tt, 
		  const double x[2][4], 
		  double shp[3][4], 
		  double &xsj, 
		  Matrix &sx ) ;

    Vector *load;
} ; 
