// $Revision: 1.19 $
// $Date: 2002-06-21 00:28:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PressureIndependMultiYield.cpp,v $
                                                                        
// Written: ZHY
// Created: August 2000

//
// PressureIndependMultiYield.cpp
// -------------------
//
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>
#include <PressureIndependMultiYield.h>
#include <Information.h>
#include <ID.h>

int PressureIndependMultiYield::loadStage = 0;
Matrix PressureIndependMultiYield::theTangent(6,6);
T2Vector PressureIndependMultiYield::subStrainRate;

PressureIndependMultiYield::PressureIndependMultiYield (int tag, int nd, 
							double r, double refShearModul,
							double refBulkModul, 
							double cohesi, double peakShearStra, 
							double frictionAng, double refPress, double pressDependCoe,
							int numberOfYieldSurf)
 : NDMaterial(tag,ND_TAG_PressureIndependMultiYield), currentStress(),
   trialStress(), currentStrain(), strainRate()
{
  if (nd !=2 && nd !=3) {
    cerr << "FATAL:PressureIndependMultiYield:: dimension error" << endl;
    cerr << "Dimension has to be 2 or 3, you give nd= " << nd << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (refShearModul <= 0) {
    cerr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: refShearModulus <= 0" << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (refBulkModul <= 0) {
    cerr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: refBulkModulus <= 0" << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (frictionAng < 0.) {
    cerr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: frictionAngle < 0" << endl;
    cerr << "Will reset frictionAngle to zero." << endl;
    frictionAng = 0.;
  }
  if (frictionAng == 0. && cohesi <= 0. ) {
    cerr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: frictionAngle && cohesion <= 0." << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (cohesi <= 0) {
    cerr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: cohesion <= 0" << endl;
    cerr << "Will reset cohesion to zero." << endl;
    cohesi = 0.;
  }
  if (peakShearStra <= 0) {
    cerr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: peakShearStra <= 0" << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (refPress <= 0) {
    cerr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: refPress <= 0" << endl;
    g3ErrorHandler->fatal(" ");
  }
  if (pressDependCoe < 0) {
    cerr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: pressDependCoe < 0" << endl;
    cerr << "Will reset pressDependCoe to zero." << endl;
    pressDependCoe = 0.;
  }
  if (numberOfYieldSurf <= 0) {
    cerr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: numberOfSurfaces <= 0" << endl;
    cerr << "Will use 10 yield surfaces." << endl;
    numberOfYieldSurf = 10;
  }
  if (numberOfYieldSurf > 40) {
    cerr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: numberOfSurfaces > 40" << endl;
    cerr << "Will use 40 yield surfaces." << endl;
    numberOfYieldSurf = 40;
  }
  if (r < 0) {
    cerr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: mass density < 0" << endl;
    cerr << "Will use rho = 0." << endl;
    r = 0.;
  }

  ndm = nd;
  loadStage = 0;   //default
  refShearModulus = refShearModul;
  refBulkModulus = refBulkModul;
  frictionAngle = frictionAng;
  peakShearStrain = peakShearStra;
  refPressure = -refPress;  //compression is negative
  cohesion = cohesi;
  pressDependCoeff = pressDependCoe;
  numOfSurfaces = numberOfYieldSurf;
  rho = r;

  e2p = 0;
  theSurfaces = new MultiYieldSurface[numOfSurfaces+1]; //first surface not used
  committedSurfaces = new MultiYieldSurface[numOfSurfaces+1]; 
  activeSurfaceNum = committedActiveSurf = 0; 

  setUpSurfaces();  // residualPress is calculated inside.
}
   

PressureIndependMultiYield::PressureIndependMultiYield () 
 : NDMaterial(0,ND_TAG_PressureIndependMultiYield), 
   currentStress(), trialStress(), currentStrain(), 
  strainRate(), theSurfaces(0), committedSurfaces(0)
{
  ndm = 3;
  loadStage = 0;   
  refShearModulus = 0.;
  refBulkModulus = 0.;
  frictionAngle = 0.;
  peakShearStrain = 0.;
  refPressure = 0.;  //compression is negative
  cohesion = 0.;
  pressDependCoeff = 0.;
  numOfSurfaces = 1;
  residualPress = 0.;
  rho = 0.;

  e2p = 0;
  activeSurfaceNum = committedActiveSurf = 0; 
}


PressureIndependMultiYield::PressureIndependMultiYield (const PressureIndependMultiYield & a)
 : NDMaterial(a.getTag(),ND_TAG_PressureIndependMultiYield), 
   currentStress(a.currentStress), trialStress(a.trialStress), 
  currentStrain(a.currentStrain), strainRate(a.strainRate)
{
  ndm = a.ndm;
  loadStage = a.loadStage;  
  refShearModulus = a.refShearModulus;
  refBulkModulus = a.refBulkModulus;
  frictionAngle = a.frictionAngle;
  peakShearStrain = a.peakShearStrain;
  refPressure = a.refPressure;
  cohesion = a.cohesion;
  pressDependCoeff = a.pressDependCoeff;
  numOfSurfaces = a.numOfSurfaces;
  residualPress = a.residualPress;
  rho = a.rho;

  e2p = a.e2p;
  committedActiveSurf = a.committedActiveSurf;
  activeSurfaceNum = a.activeSurfaceNum; 

  theSurfaces = new MultiYieldSurface[numOfSurfaces+1];  //first surface not used
  committedSurfaces = new MultiYieldSurface[numOfSurfaces+1];  
  for(int i=1; i<=numOfSurfaces; i++) {
    committedSurfaces[i] = a.committedSurfaces[i];  
    theSurfaces[i] = a.theSurfaces[i];  
  }
}


PressureIndependMultiYield::~PressureIndependMultiYield ()
{
  if (theSurfaces != 0) delete [] theSurfaces;
  if (committedSurfaces != 0) delete [] committedSurfaces;
}


void PressureIndependMultiYield::elast2Plast(void)
{
  if (loadStage == 0 || e2p == 1) return;
  e2p = 1;

  if (currentStress.volume() > 0. && frictionAngle > 0.) {
    //cerr << "WARNING:PressureIndependMultiYield::elast2Plast(): material in tension." << endl;
    currentStress.setData(currentStress.deviator(),0);
  }

  paramScaling();  // scale surface parameters corresponding to initial confinement

  // Active surface is 0, return
  if (currentStress.deviatorLength() == 0.) return;

  // Find active surface
  while (yieldFunc(currentStress, committedSurfaces, ++committedActiveSurf) > 0) {
    if (committedActiveSurf == numOfSurfaces) {
      //cerr <<"WARNING:PressureIndependMultiYield::elast2Plast(): stress out of failure surface"<<endl;
      deviatorScaling(currentStress, committedSurfaces, numOfSurfaces);
      initSurfaceUpdate();
      return;
    }
  } 
  committedActiveSurf--;
  initSurfaceUpdate();
}


int PressureIndependMultiYield::setTrialStrain (const Vector &strain)
{
  static Vector temp(6);
  if (ndm==3 && strain.Size()==6) 
    temp = strain;
  else if (ndm==2 && strain.Size()==3) {
    temp[0] = strain[0];
    temp[1] = strain[1];
    temp[2] = 0.0;
    temp[3] = strain[2];
    temp[4] = 0.0;
    temp[5] = 0.0;
  }
  else {
    cerr << "Fatal:D2PressDepMYS:: Material dimension is: " << ndm << endl;
    cerr << "But strain vector size is: " << strain.Size() << endl;
    g3ErrorHandler->fatal(" ");
  }
	
  //strainRate.setData(temp-currentStrain.t2Vector(1),1);
  temp -= currentStrain.t2Vector(1);
  strainRate.setData(temp, 1);
	
  return 0;
}


int PressureIndependMultiYield::setTrialStrain (const Vector &strain, const Vector &rate)
{
  return setTrialStrain (strain);
}


int PressureIndependMultiYield::setTrialStrainIncr (const Vector &strain)
{
  static Vector temp(6);
  if (ndm==3 && strain.Size()==6) 
    temp = strain;
  else if (ndm==2 && strain.Size()==3) {
    temp[0] = strain[0];
    temp[1] = strain[1];
    temp[3] = strain[2];
  }
  else {
    cerr << "Fatal:D2PressDepMYS:: Material dimension is: " << ndm << endl;
    cerr << "But strain vector size is: " << strain.Size() << endl;
    g3ErrorHandler->fatal(" ");
  }
  
  strainRate.setData(temp,1);
  return 0;
}


int PressureIndependMultiYield::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return setTrialStrainIncr(strain);
}


const Matrix & PressureIndependMultiYield::getTangent (void)
{
  if (loadStage != 0 && e2p == 0) elast2Plast();

  if (loadStage==0) {  //linear elastic
    for (int i=0;i<6;i++) 
      for (int j=0;j<6;j++) {
	theTangent(i,j) = 0.;
	if (i==j) theTangent(i,j) += refShearModulus;
	if (i<3 && j<3 && i==j) theTangent(i,j) += refShearModulus;
	if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
      }
  }
  else {
    double coeff;
    static Vector devia(6);
    
    if (committedActiveSurf > 0) {
      //devia = currentStress.deviator()-committedSurfaces[committedActiveSurf].center();
      devia = currentStress.deviator();
      devia -= committedSurfaces[committedActiveSurf].center();
	
      double size = committedSurfaces[committedActiveSurf].size();
      double plastModul = committedSurfaces[committedActiveSurf].modulus();
      coeff = 6.*refShearModulus*refShearModulus/(2.*refShearModulus+plastModul)/size/size;
    }
    else coeff = 0.;
    
    for (int i=0;i<6;i++) 
      for (int j=0;j<6;j++) {
	theTangent(i,j) = - coeff*devia[i]*devia[j];
        if (i==j) theTangent(i,j) += refShearModulus;
        if (i<3 && j<3 && i==j) theTangent(i,j) += refShearModulus;
	if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
      }
  }

  if (ndm==3) 
    return theTangent;
  else {
  	static Matrix workM(3,3);
    workM(0,0) = theTangent(0,0);
    workM(0,1) = theTangent(0,1);
    workM(0,2) = theTangent(0,3);
    workM(1,0) = theTangent(1,0);
    workM(1,1) = theTangent(1,1);
    workM(1,2) = theTangent(1,3);
    workM(2,0) = theTangent(3,0);
    workM(2,1) = theTangent(3,1);
    workM(2,2) = theTangent(3,3);
    return workM;
  }
}


const Vector & PressureIndependMultiYield::getStress (void)
{
  int i;
  if (loadStage != 0 && e2p == 0) elast2Plast();

  if (loadStage==0) {  //linear elastic
    static Vector trialStrain(6);
    //trialStrain = currentStrain.t2Vector(1) + strainRate.t2Vector(1);
    trialStrain = currentStrain.t2Vector(1);
    trialStrain += strainRate.t2Vector(1);
    
    getTangent();
    static Vector a(6);
    a.addMatrixVector(0.0, theTangent, trialStrain, 1.0);
    trialStress.setData(a);
  }

  else {
    for (i=1; i<=numOfSurfaces; i++) theSurfaces[i] = committedSurfaces[i];
    activeSurfaceNum = committedActiveSurf;
    subStrainRate = strainRate;
    setTrialStress(currentStress);
    if (isLoadReversal()) {
      updateInnerSurface();
      activeSurfaceNum = 0;
    }
    int numSubIncre = setSubStrainRate();
    
    for (i=0; i<numSubIncre; i++) {
      if (i==0)  
	setTrialStress(currentStress);
      else 
	setTrialStress(trialStress);
      if (activeSurfaceNum==0 && !isCrossingNextSurface()) continue;
      if (activeSurfaceNum==0) activeSurfaceNum++;
      stressCorrection(0);
      updateActiveSurface();
    }
    //volume stress change
    double volum = refBulkModulus*(strainRate.volume()*3.);
    volum += currentStress.volume();
    if (volum > 0) volum = 0.;
    trialStress.setData(trialStress.deviator(),volum);
  }

  if (ndm==3)
    return trialStress.t2Vector();
  else {
    static Vector workV(3);
    workV[0] = trialStress.t2Vector()[0];
    workV[1] = trialStress.t2Vector()[1];
    workV[2] = trialStress.t2Vector()[3];
    return workV;
  }
}


const Vector & PressureIndependMultiYield::getStrain (void)
{
  return getCommittedStrain();
}


int PressureIndependMultiYield::commitState (void)
{
  currentStress = trialStress;
  
  //currentStrain = T2Vector(currentStrain.t2Vector() + strainRate.t2Vector());
  static Vector temp(6);
  temp = currentStrain.t2Vector();
  temp += strainRate.t2Vector();
  currentStrain.setData(temp);
  temp.Zero();
  strainRate.setData(temp);
  
  if (loadStage) {
    committedActiveSurf = activeSurfaceNum;
    for (int i=1; i<=numOfSurfaces; i++) committedSurfaces[i] = theSurfaces[i];
  }

  return 0;
}
 

int PressureIndependMultiYield::revertToLastCommit (void)
{
  return 0;
}


NDMaterial * PressureIndependMultiYield::getCopy (void)
{
  PressureIndependMultiYield * copy = new PressureIndependMultiYield(*this);
  return copy;
}


NDMaterial * PressureIndependMultiYield::getCopy (const char *code)
{
  if (strcmp(code,"PressureIndependMultiYield") == 0 || strcmp(code,"PlaneStrain") == 0
      || strcmp(code,"ThreeDimensional") == 0) {
    PressureIndependMultiYield * copy = new PressureIndependMultiYield(*this);
    return copy;
  }

  return 0;
}


const char * PressureIndependMultiYield::getType (void) const
{
  return (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
}


int PressureIndependMultiYield::getOrder (void) const
{
  return (ndm == 2) ? 3 : 6;
}


int PressureIndependMultiYield::updateParameter(int responseID, Information &info)
{
  loadStage = responseID;
  return 0;
}


int PressureIndependMultiYield::sendSelf(int commitTag, Channel &theChannel)
{
  int i, res = 0;

  static ID idData(4);
  idData(0) = this->getTag();
  idData(1) = numOfSurfaces;
  idData(2) = loadStage;
  idData(3) = ndm;

  res += theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    g3ErrorHandler->warning("%s -- could not send Vector",
			    "PressureDependMultiYield::sendSelf");
    return res;
  }

  Vector data(23+numOfSurfaces*8);
  static Vector temp(6);
  data(0) = rho;
  data(1) = refShearModulus;
  data(2) = refBulkModulus;
  data(3) = frictionAngle;
  data(4) = peakShearStrain;
  data(5) = refPressure;
  data(6) = cohesion;
  data(7) = pressDependCoeff;
  data(8) = residualPress;
  data(9) = e2p;
  data(10) = committedActiveSurf;
	
  temp = currentStress.t2Vector();
  for(i = 0; i < 6; i++) data(i+11) = temp[i];
  
  temp = currentStrain.t2Vector();
  for(i = 0; i < 6; i++) data(i+17) = temp[i];
  
  for(i = 0; i < numOfSurfaces; i++) {
    int k = 23 + i*8;
    data(k) = committedSurfaces[i+1].size();
    data(k+1) = committedSurfaces[i+1].modulus();
    temp = committedSurfaces[i+1].center();
    data(k+2) = temp(0);
    data(k+3) = temp(1);
    data(k+4) = temp(2);
    data(k+5) = temp(3);
    data(k+6) = temp(4);
    data(k+7) = temp(5);
  }

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    g3ErrorHandler->warning("%s -- could not send Vector",
			    "FluidSolidPorousMaterial::sendSelf");
    return res;
  }
  
  return res;
}


int PressureIndependMultiYield::recvSelf(int commitTag, Channel &theChannel, 
					 FEM_ObjectBroker &theBroker)    
{
  int i, res = 0;

  static ID idData(4);
  idData(0) = this->getTag();
  idData(1) = numOfSurfaces;
  idData(2) = loadStage;
  idData(3) = ndm;

  res += theChannel.recvID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    g3ErrorHandler->warning("%s -- could not send Vector",
			    "PressureDependMultiYield::sendSelf");
    return res;
  }

  this->setTag((int)idData(0));
  loadStage = idData(2);
  ndm = idData(3);

  Vector data(23+idData(1)*8);
  static Vector temp(6);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    g3ErrorHandler->warning("%s -- could not receive Vector",
			    "FluidSolidPorousMaterial::recvSelf");
    return res;
  }
    
  rho = data(0);
  refShearModulus = data(1);
  refBulkModulus = data(2);
  frictionAngle = data(3);
  peakShearStrain = data(4);
  refPressure = data(5);
  cohesion = data(6);
  pressDependCoeff = data(7);
  residualPress = data(8);
  e2p = data(9);
  committedActiveSurf = data(10);
  
  for(i = 0; i < 6; i++) temp[i] = data(i+11);
  currentStress.setData(temp);
  
  for(i = 0; i < 6; i++) temp[i] = data(i+17);
  currentStrain.setData(temp);

  if (numOfSurfaces != idData(1)) {
    if (committedSurfaces != 0) {
      delete [] committedSurfaces;
      delete [] theSurfaces;
    }
    numOfSurfaces = idData(1);
    theSurfaces = new MultiYieldSurface[numOfSurfaces+1]; //first surface not used
    committedSurfaces = new MultiYieldSurface[numOfSurfaces+1]; 

    for (int i=1; i<=numOfSurfaces; i++) 
      committedSurfaces[i] = MultiYieldSurface();    
  }
  
  for(i = 0; i < numOfSurfaces; i++) {
    int k = 23 + i*8;
    temp(0) = data(k+2);
    temp(1) = data(k+3);
    temp(2) = data(k+4);
    temp(3) = data(k+5);
    temp(4) = data(k+6);
    temp(5) = data(k+7);
    committedSurfaces[i+1].setData(temp, data(k), data(k+1));
  }
  
  return res;
}


int PressureIndependMultiYield::getResponse (int responseID, Information &matInfo)
{
	switch (responseID) {
		case -1:
			return -1;
		case 1:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getCommittedStress();
			return 0;
		case 2:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getCommittedStrain();
			return 0;
		case 3:
			if (matInfo.theMatrix != 0)
				*(matInfo.theMatrix) = getTangent();
			return 0;
		default:
			return -1;
	}
}


void PressureIndependMultiYield::Print(ostream &s, int flag )
{
  s << "PressureIndependMultiYield" << endl;
}


const Vector & PressureIndependMultiYield::getCommittedStress (void)
{
	double scale = sqrt(3./2.)*currentStress.deviatorLength()/committedSurfaces[numOfSurfaces].size();
  if (ndm==3) {
		static Vector temp7(7), temp6(6);
		temp6 = currentStress.t2Vector();
    temp7[0] = temp6[0];
    temp7[1] = temp6[1];
    temp7[2] = temp6[2];
    temp7[3] = temp6[3];
    temp7[4] = temp6[4];
    temp7[5] = temp6[5];
    temp7[6] = scale;
		return temp7;
	}
  else {
    static Vector temp4(4), temp6(6);
		temp6 = currentStress.t2Vector();
    temp4[0] = temp6[0];
    temp4[1] = temp6[1];
    temp4[2] = temp6[3];
    temp4[3] = scale;
    return temp4;
  }
}


const Vector & PressureIndependMultiYield::getCommittedStrain (void)
{	
  if (ndm==3)
    return currentStrain.t2Vector(1);
  else {
    static Vector workV(3), temp6(6);
		temp6 = currentStrain.t2Vector(1);
    workV[0] = temp6[0];
    workV[1] = temp6[1];
    workV[2] = temp6[3];
    return workV;
  }
}


// NOTE: surfaces[0] is not used 
void PressureIndependMultiYield::setUpSurfaces (void)
{ 
	double pi = 3.14159265358979;
	double refStrain, peakShear, coneHeight;

	if (frictionAngle > 0) {
		double sinPhi = sin(frictionAngle * pi/180.);
		double Mnys = 6.*sinPhi/(3.-sinPhi);
		residualPress = 3.* cohesion / (sqrt(2.) * Mnys);
		coneHeight = - (refPressure - residualPress);
		peakShear = sqrt(2.) * coneHeight * Mnys / 3.; 
		refStrain = (peakShearStrain * peakShear) 
			          / (refShearModulus * peakShearStrain - peakShear);
	}

	else if (frictionAngle == 0.) { // cohesion = peakShearStrength
		peakShear = cohesion;
		refStrain = (peakShearStrain * peakShear) 
			          / (refShearModulus * peakShearStrain - peakShear);
		residualPress = 0.;
	}

	double  stress1, stress2, strain1, strain2, size, elasto_plast_modul, plast_modul;	
	double stressInc = peakShear / numOfSurfaces;

	for (int ii=1; ii<=numOfSurfaces; ii++){
        stress1 = ii * stressInc; 
				stress2 = stress1 + stressInc;
        strain1 = stress1 * refStrain / (refShearModulus * refStrain - stress1);
        strain2 = stress2 * refStrain / (refShearModulus * refStrain - stress2);
        if (frictionAngle > 0.) size = 3. * stress1 / sqrt(2.) / coneHeight;
        else if (frictionAngle == 0.) size = 3. * stress1 / sqrt(2.);
 
        elasto_plast_modul = 2.*(stress2 - stress1)/(strain2 - strain1);

        if ( (2.*refShearModulus - elasto_plast_modul) <= 0) 
					plast_modul = UP_LIMIT;
        else 
					plast_modul = (2.*refShearModulus * elasto_plast_modul)/
                        (2.*refShearModulus - elasto_plast_modul);
        if (plast_modul < 0) plast_modul = 0;
        if (plast_modul > UP_LIMIT) plast_modul = UP_LIMIT;
        if (ii==numOfSurfaces) plast_modul = 0;

		static Vector temp(6);
        committedSurfaces[ii] = MultiYieldSurface(temp,size,plast_modul);
  }  // ii   
}


double PressureIndependMultiYield::yieldFunc(const T2Vector & stress, 
											 const MultiYieldSurface * surfaces, int surfaceNum)
{
	static Vector temp(6);
	//temp = stress.deviator() - surfaces[surfaceNum].center();
	temp = stress.deviator();
	temp -= surfaces[surfaceNum].center();

	double sz = surfaces[surfaceNum].size();
	return 3./2.*(temp && temp) - sz * sz;
}


void PressureIndependMultiYield::deviatorScaling(T2Vector & stress, const MultiYieldSurface * surfaces, 
																			int surfaceNum)
{
	double diff = yieldFunc(stress, surfaces, surfaceNum);

	if ( surfaceNum < numOfSurfaces && diff < 0. ) {
		double sz = surfaces[surfaceNum].size();
		double deviaSz = sqrt(sz*sz + diff);
		static Vector devia(6);
		devia = stress.deviator(); 
		static Vector temp(6);
		temp = devia - surfaces[surfaceNum].center();
		double coeff = (sz-deviaSz) / deviaSz;
		if (coeff < 1.e-13) coeff = 1.e-13;
		devia.addVector(1.0, temp, coeff);
		stress.setData(devia, stress.volume());

		deviatorScaling(stress, surfaces, surfaceNum);  // recursive call
	}

	if (surfaceNum==numOfSurfaces && fabs(diff) > LOW_LIMIT) {
		double sz = surfaces[surfaceNum].size();
		static Vector newDevia(6);
		newDevia.addVector(0.0, stress.deviator(), sz/sqrt(diff+sz*sz));
		stress.setData(newDevia, stress.volume());
	}
}


void PressureIndependMultiYield::initSurfaceUpdate()
{
	if (activeSurfaceNum == 0) return; 

	static Vector devia(6);
	devia = currentStress.deviator();
	double Ms = sqrt(3./2.*(devia && devia));
	static Vector newCenter(6);

	if (activeSurfaceNum < numOfSurfaces) { // failure surface can't move
		//newCenter = devia * (1. - committedSurfaces[activeSurfaceNum].size() / Ms); 
		newCenter.addVector(0.0, devia, 1.0-committedSurfaces[activeSurfaceNum].size()/Ms);
		committedSurfaces[activeSurfaceNum].setCenter(newCenter);
	}

	for (int i=1; i<activeSurfaceNum; i++) {
	  newCenter = devia * (1. - committedSurfaces[i].size() / Ms);
	  committedSurfaces[i].setCenter(newCenter); 
	}
}


void PressureIndependMultiYield::paramScaling(void)
{
	if (frictionAngle == 0.) return;

	double conHeig = - (currentStress.volume() - residualPress);
	double scale = -conHeig / (refPressure-residualPress);
           
	scale = pow(scale, pressDependCoeff); 
	refShearModulus *= scale;
   
	double plastModul, size;
	static Vector temp(6);
	for (int i=1; i<=numOfSurfaces; i++) {
	  plastModul = committedSurfaces[i].modulus() * scale;
	  size = committedSurfaces[i].size() * conHeig;
	  committedSurfaces[i] =  MultiYieldSurface(temp,size,plastModul);
	}
}


void PressureIndependMultiYield::setTrialStress(T2Vector & stress)
{
  static Vector devia(6);
  //devia = stress.deviator() + subStrainRate.deviator()*2.*refShearModulus;
  devia = stress.deviator();
  devia.addVector(1.0, subStrainRate.deviator(), 2.0*refShearModulus);
  
  trialStress.setData(devia, stress.volume());
}


int PressureIndependMultiYield::setSubStrainRate(void)
{
	if (activeSurfaceNum==numOfSurfaces) return 1;

	//if (strainRate==T2Vector()) return 0;
	if (strainRate.isZero()) return 0;

	double elast_plast_modulus;
	if (activeSurfaceNum==0) 
	  elast_plast_modulus = 2*refShearModulus;
	else {
	  double plast_modulus = theSurfaces[activeSurfaceNum].modulus();
	  elast_plast_modulus = 2*refShearModulus*plast_modulus 
	    / (2*refShearModulus+plast_modulus);
	}
	static Vector incre(6);
	//incre = strainRate.deviator()*elast_plast_modulus;
	incre.addVector(0.0, strainRate.deviator(),elast_plast_modulus);

	static T2Vector increStress;
	increStress.setData(incre, 0);
	double singleCross = theSurfaces[numOfSurfaces].size() / numOfSurfaces;
	double totalCross = 3.*increStress.octahedralShear() / sqrt(2.);
	int numOfSub = totalCross/singleCross + 1;
	if (numOfSub > numOfSurfaces) numOfSub = numOfSurfaces;
	//incre = strainRate.t2Vector() / numOfSub;
	incre = strainRate.t2Vector();
	incre /= numOfSub;
	subStrainRate.setData(incre);

	return numOfSub;
}


void
PressureIndependMultiYield::getContactStress(T2Vector &contactStress)
{
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center(); 
	static Vector devia(6);
	//devia = trialStress.deviator() - center;
	devia = trialStress.deviator();
	devia -= center;

	double Ms = sqrt(3./2.*(devia && devia));
	//devia = devia * theSurfaces[activeSurfaceNum].size() / Ms + center;
	devia *= theSurfaces[activeSurfaceNum].size() / Ms;
	devia += center;

	contactStress.setData(devia,trialStress.volume()); 
}


int PressureIndependMultiYield::isLoadReversal(void)
{
  if(activeSurfaceNum == 0) return 0;

  static Vector surfaceNormal(6);
  getSurfaceNormal(currentStress, surfaceNormal);
 
  //(((trialStress.deviator() - currentStress.deviator()) && surfaceNormal) < 0) 
  // return 1;
  static Vector a(6);
  a = trialStress.deviator();
  a-= currentStress.deviator();
  if((a && surfaceNormal) < 0) 
    return 1;

  return 0;   
}


void
PressureIndependMultiYield::getSurfaceNormal(const T2Vector & stress, Vector &surfaceNormal)
{
  //Q = stress.deviator() - theSurfaces[activeSurfaceNum].center();
  // return Q / sqrt(Q && Q);

  surfaceNormal = stress.deviator();
  surfaceNormal -= theSurfaces[activeSurfaceNum].center();
  surfaceNormal /= sqrt(surfaceNormal && surfaceNormal);
}


double PressureIndependMultiYield::getLoadingFunc(const T2Vector & contactStress, 
																			 const Vector & surfaceNormal,
																			 int crossedSurface)
{
	double loadingFunc;
  double temp1 = 2. * refShearModulus ;
  double temp2 = theSurfaces[activeSurfaceNum].modulus();

  //for crossing first surface
  double temp = temp1 + temp2;
  //loadingFunc = (surfaceNormal && (trialStress.deviator()-contactStress.deviator()))/temp;
  static Vector tmp(6);
  tmp =trialStress.deviator();
  tmp -= contactStress.deviator();
  loadingFunc = (surfaceNormal && tmp)/temp;
   //for crossing more than one surface
  if(crossedSurface) {
    double temp3 = theSurfaces[activeSurfaceNum-1].modulus();
    loadingFunc *= (temp3 - temp2)/temp3;
  }

  return loadingFunc;
}


void PressureIndependMultiYield::stressCorrection(int crossedSurface)
{
	static T2Vector contactStress;
	this->getContactStress(contactStress);
	static Vector surfaceNormal(6);
	this->getSurfaceNormal(contactStress, surfaceNormal);
	double loadingFunc = getLoadingFunc(contactStress, surfaceNormal, crossedSurface);
	static Vector devia(6);

	//devia = trialStress.deviator() - surfaceNormal * 2 * refShearModulus * loadingFunc;
	devia.addVector(0.0, surfaceNormal, -2*refShearModulus*loadingFunc);
	devia += trialStress.deviator();

	trialStress.setData(devia, trialStress.volume());
	deviatorScaling(trialStress, theSurfaces, activeSurfaceNum);

	if (isCrossingNextSurface()) {
		activeSurfaceNum++;
		stressCorrection(1);  //recursive call
	}
}


void PressureIndependMultiYield::updateActiveSurface(void)
{
  if (activeSurfaceNum == numOfSurfaces) return;

	double A, B, C, X;
	static T2Vector direction;
	static Vector t1(6);
	static Vector t2(6);
	static Vector temp(6);
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	static Vector outcenter(6);
	outcenter= theSurfaces[activeSurfaceNum+1].center();
	double outsize = theSurfaces[activeSurfaceNum+1].size();


	//t1 = trialStress.deviator() - center;
	//t2 = center - outcenter;
	t1 = trialStress.deviator();
	t1 -= center;
	t2 = center;
	t2 -= outcenter;

	A = t1 && t1;
	B = 2. * (t1 && t2);
	C = (t2 && t2) - 2./3.* outsize * outsize;
	X = secondOrderEqn(A,B,C,0);
	if ( fabs(X-1.) < LOW_LIMIT ) X = 1.;
	if (X < 1.){
	  cerr << "FATAL:PressureIndependMultiYield::updateActiveSurface(): error in Direction of surface motion." 
	       << endl; 
	  g3ErrorHandler->fatal(" ");
	}

	//temp = (t1 * X + center) * (1. - size / outsize) - (center - outcenter * size / outsize);
	temp = center;
	temp.addVector(1.0, t1, X);
	temp *= (1.0 - size/outsize);
	t2 = center;
	t2.addVector(1.0, outcenter, -size/outsize);
	temp -= t2;

	direction.setData(temp);

	if (direction.deviatorLength() < LOW_LIMIT) return;

	temp = direction.deviator();  
	A = temp && temp;
	B = - 2 * (t1 && temp);
	if (fabs(B) < LOW_LIMIT) B = 0.; 
	C = (t1 && t1) - 2./3.* size * size;
	if ( fabs(C) < LOW_LIMIT || fabs(C)/(t1 && t1) < LOW_LIMIT ) return;

	if (B > 0. || C < 0.) {
	  cerr << "FATAL:PressureIndependMultiYield::updateActiveSurface(): error in surface motion.\n" 
	       << "A= " <<A <<" B= " <<B <<" C= "<<C <<" (t1&&t1)= "<<(t1&&t1) <<endl; 
	  g3ErrorHandler->fatal(" ");
	}
	X = secondOrderEqn(A,B,C,1);  

	//center += temp * X;
	center.addVector(1.0, temp, X);
	theSurfaces[activeSurfaceNum].setCenter(center);
}      


void PressureIndependMultiYield::updateInnerSurface(void)
{
	if (activeSurfaceNum <= 1) return;

	static Vector devia(6);
	devia = currentStress.deviator();
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	static Vector newcenter(6);

	for (int i=1; i<activeSurfaceNum; i++) {
		//newcenter = devia - (devia - center) * theSurfaces[i].size() / size;
		newcenter = center; 
		newcenter -= devia;
		newcenter *= theSurfaces[i].size()/size;
		newcenter += devia;

		theSurfaces[i].setCenter(newcenter);
	}
}


int PressureIndependMultiYield:: isCrossingNextSurface(void)
{
  if (activeSurfaceNum == numOfSurfaces) return 0;  

  if(yieldFunc(trialStress, theSurfaces, activeSurfaceNum+1) > 0) return 1;
  
  return 0;
}
 
