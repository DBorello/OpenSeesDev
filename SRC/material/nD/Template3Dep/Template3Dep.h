/*
################################################################################
# COPYRIGHT (C):     :-))                                                      #
# PROJECT:           Object Oriented Finite Element Program                    #
# PURPOSE:           General platform for elaso-plastic constitutive model     #
#                    implementation                                            #
# CLASS:             Template3Dep (the base class for all material point)      #
#                                                                              #
# VERSION:                                                                     #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
# TARGET OS:         DOS || UNIX || . . .                                      #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                               #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                               #
#                                                                              #
#                                                                              #
# DATE:              08-03-2000                                                #
# UPDATE HISTORY:    09-12-2000                                                #
#                                                                              #
#                                                                              #
# SHORT EXPLANATION: The Template3Dep class is used to hold specific           #
#                    yield surface, potential surface, Evolution law(s)        #
#                    and EPState of a 3D elasto-plastic material model for one #
#                    gauss point!! It is worthwhile noting that one model may  #
#                    have multiple evolution law.  Each evlotuion law is       #
#                    used to evolve one internal var.                          #
#                                                                              #
#                                                                              #
################################################################################
*/

#ifndef Template3Dep_H
#define Template3Dep_H

#include <NDMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <stresst.h>
#include <straint.h>
#include <BJtensor.h>

#include <MD_YS.h>
#include <MD_PS.h>
#include <EL_S.h>
#include <EL_T.h>
#include <EPState.h>
//#include <CDriver.h>

class Template3Dep : public NDMaterial
{
 public:

  public:
    // constructor
    Template3Dep(  int tag               ,
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
              	   EPState          *EPS_,
	       	   EvolutionLaw_S   *ELS1_ , 
	       	   EvolutionLaw_S   *ELS2_ , 
	       	   EvolutionLaw_S   *ELS3_ , 
	       	   EvolutionLaw_S   *ELS4_ , 
	       	   EvolutionLaw_T   *ELT1_ ,
	       	   EvolutionLaw_T   *ELT2_ ,
	       	   EvolutionLaw_T   *ELT3_ ,
	       	   EvolutionLaw_T   *ELT4_  );
    
    // Constructor0
    // If no evolution law is provided, then there will be no hardening or softening!
    Template3Dep(  int tag               ,
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
              	   EPState          *EPS_);

    // Constructor1
    // Only one scalar evolution law is provided!
    Template3Dep(  int tag               ,
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
              	   EPState          *EPS_,
	       	   EvolutionLaw_S   *ELS1_ );

    // Constructor2
    // Only one tensorial evolution law is provided!
    Template3Dep(  int tag               ,
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
              	   EPState          *EPS_,
	       	   EvolutionLaw_T   *ELT1_ );

    // Constructor 3
    // One scalar evolution law and one tensorial evolution law are provided!
    Template3Dep(  int tag               , 
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
              	   EPState          *EPS_,
	       	   EvolutionLaw_S   *ELS1_, 
	       	   EvolutionLaw_T   *ELT1_ );
    
    // Constructor 4
    // Two scalar evolution laws and one tensorial evolution law are provided!
    Template3Dep(  int tag               ,
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
              	   EPState          *EPS_,
	       	   EvolutionLaw_S   *ELS1_, 
	       	   EvolutionLaw_S   *ELS2_, 
	       	   EvolutionLaw_T   *ELT1_ );

    // Constructor 5
    // Two scalar evolution laws and two tensorial evolution laws are provided!
    Template3Dep(  int tag               ,
                   YieldSurface     *YS_ ,        
                   PotentialSurface *PS_ ,
              	   EPState          *EPS_,
	       	   EvolutionLaw_S   *ELS1_, 
	       	   EvolutionLaw_S   *ELS2_, 
	       	   EvolutionLaw_T   *ELT1_,
	       	   EvolutionLaw_T   *ELT2_ );

    // For parallel processing
    Template3Dep(void);
    virtual ~Template3Dep(void);

    // methods to set state and retrieve state using Matrix and Vector classes
    int setTrialStrain(const Vector &v);
    int setTrialStrain(const Vector &v, const Vector &r);
    int setTrialStrainIncr(const Vector &v) ;
    int setTrialStrainIncr(const Vector &v, const Vector &r) ;
    const Matrix &getTangent(void) ;
    const Vector &getStress(void) ;
    const Vector &getStrain(void) ;

    // methods to set and retrieve state using the Tensor class    
    int setTrialStrain(const Tensor &v) ;
    int setTrialStrain(const Tensor &v, const Tensor &r) ;    
    int setTrialStrainIncr(const Tensor &v) ;
    int setTrialStrainIncr(const Tensor &v, const Tensor &r) ;
    const Tensor &getTangentTensor(void) ;
    const Tensor &getStressTensor(void) ;
    const Tensor &getStrainTensor(void) ;

    int commitState(void) ;
    int revertToLastCommit(void) ;
    int revertToStart(void) ;
    
    NDMaterial *getCopy(void) ;
    NDMaterial *getCopy(const char *code) ;

    const char *getType(void) const ;
    int getOrder(void) const ;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel,  
    FEM_ObjectBroker &theBroker);    

    void Print(ostream &s, int flag =0);

    tensor ElasticComplianceTensor(void) const;
    tensor ElasticStiffnessTensor(void) const;

    // Utility method
    YieldSurface *getYS() const;
    PotentialSurface *getPS() const;
    EPState *getEPS() const;
    
    //get scalar evolution laws
    EvolutionLaw_S *getELS1() const;
    EvolutionLaw_S *getELS2() const;
    EvolutionLaw_S *getELS3() const;
    EvolutionLaw_S *getELS4() const;

    //get tensorial evolution laws
    EvolutionLaw_T *getELT1() const;
    EvolutionLaw_T *getELT2() const;
    EvolutionLaw_T *getELT3() const;
    EvolutionLaw_T *getELT4() const;
    
    void setEPS( EPState &eps);

    //================================================================================
    // Overloaded Insertion Operator
    //================================================================================
    friend ostream& operator<< (ostream& os, const Template3Dep & MP);

   private:

    YieldSurface *YS;

    PotentialSurface *PS;

    EPState *EPS;    
                               
    //Scalar variable evolution laws (currently at most 4  allowed)
    EvolutionLaw_S *ELS1; 
    EvolutionLaw_S *ELS2; 
    EvolutionLaw_S *ELS3; 
    EvolutionLaw_S *ELS4; 
    
    //Tensorial variable evolution laws (currently at most 4  allowed)
    EvolutionLaw_T *ELT1; 
    EvolutionLaw_T *ELT2; 
    EvolutionLaw_T *ELT3; 
    EvolutionLaw_T *ELT4; 

};


#endif

