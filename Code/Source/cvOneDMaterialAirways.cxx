# include <cstring>
# include <cstdlib>
# include <cmath>
# include <iostream>

# include "cvOneDMaterialAirways.h"
# include "cvOneDUtility.h"

const double PI = 4.0 * atan(1.0);

using namespace std;

cvOneDMaterialAirways::cvOneDMaterialAirways(){
  p1_ = 1333.2237 * 85.0;
  zstar = 70.0;
  alpha = 9.29e-3/1333.2237;
  pcrit = 33250.0;

  //default material values
  K1_ = 2.00e7;
  K2_ = -22.5267;
 // K3_ = 8.65e5;

  rigid = 0;
}

cvOneDMaterialAirways::~cvOneDMaterialAirways(){
}

cvOneDMaterialAirways::cvOneDMaterialAirways (const cvOneDMaterialAirways &rhs){
 // this seems not used
  cvOneDMaterial::operator=(rhs);
  double k1 = 0,k2 = 0,k3 = 0;
  double pref=0;
  rhs.GetParams( &k1, &k2, &k3, &pref);
  K1_ = k1;
  K2_ = k2;
  K3_ = k3;
  p1_= pref;
}

cvOneDMaterialAirways& cvOneDMaterialAirways::operator= (const cvOneDMaterialAirways &that){
  if (this != &that) {
    cvOneDMaterial::operator=(that);
    double k1 = 0,k2 = 0,k3 = 0;

    double pref=0;
    that.GetParams( &k1, &k2, &k3, &pref);
    K1_ = k1;
    K2_ = k2;
    K3_ = k3;
    p1_= pref;
  }
  return *this;
}

void cvOneDMaterialAirways::SetMaterialType(double *mType,double Pref){
  K1_ = mType[0];
  K2_ = mType[1];
  K3_ = mType[2];
  PP1_=Pref;
  cout << "Setting material K's " << K1_ << " " << K2_ << " " << K3_ << " ..." << endl;
  cout << "Setting reference Pressure " << PP1_ << endl;
}

void cvOneDMaterialAirways::SetPeriod(double per){
  Period=per;
}

double cvOneDMaterialAirways::GetProperty(char* what)const{
  if(strcmp(what, "density") == 0)
    return density;
  else if( strcmp(what, "dynamic viscosity") == 0)
    return dynamicViscosity;
  else if( strcmp(what, "kinematic viscosity") == 0)
    return kinematicViscosity;
  else if( strcmp(what, "constant n") == 0)
    return profile_exponent;
  else if( strcmp(what, "delta") == 0)
    return delta;
  else if( strcmp(what, "N") == 0)
    return N;
  else{
    abort();
    return 0.0;
  }
}

double cvOneDMaterialAirways::GetEHR(double z)const{
  double ro = Getr1(z);
  double ans = 4./3.*(K1_*ro + K3_) ;//dyne/cm^2=g/cm/s^2
  if (fabs(ans) < 1e-12) {
    cerr << "ERROR: EHR ~ 0 at z=" << z << endl;
    abort();
  }
  return ans;
}

void cvOneDMaterialAirways::SetAreas_and_length_and_in_out(double S_top,double S_bottom,double z,double z_in,double z_out){
  Stop = S_top;//inlet
  Sbot = S_bottom;//outlet
  len = z;
  Zin = z_in;
  Zout = z_out;
  lt = len/6.0;
  double eps = 1.0e-12;
  double dS = fabs(Sbot - Stop);
  //if (dS < 2.0e-5) {
  //  cerr << "ERROR: Sbot ≈ Stop → atanh invalid" << endl;
  //  abort();
  //}
  if (dS < 1e-6) {
    // ===== Tube droit =====
    chiS = 0.0;
  } else {
    double eps = 1.0e-12;
    chiS = 2.0 * atanh(1.0 - (2.0e-5 / dS)) / lt + eps;
  }
}

double cvOneDMaterialAirways::GetS1(double z)const{
  double area;
  area = (Sbot - Stop) * ((1.0 - tanh(chiS * (z - Zin - 0.5*lt)))/2.0) + Stop;
  if (area <= 0.0) {
    cerr << "ERROR: So_ <= 0 at z=" << z << endl;
    abort();
  }
  return area;  // slightly increased pressure/decreased area
}

double cvOneDMaterialAirways::Getr1(double z)const{
  // linearly interpolated r
  double area = GetS1(z);
  double r = sqrt(area / PI);
  return r;
}

double cvOneDMaterialAirways::GetDS1Dz(double z)const{
  double dsdz = 0.5*(Sbot - Stop)*chiS*(1 - pow(tanh(chiS*(z - Zin - 0.5*lt)),2));
  return dsdz;
}

//this is in the reference state dr1dz
double cvOneDMaterialAirways::GetDr1Dz(double z)const{
  // linearly vary radius
  double dsdz = GetDS1Dz(z);
  double So_ = GetS1(z);
  double drdz = -(sqrt(PI)*dsdz)/(2.0*So_*sqrt(So_));
  return drdz;
}

double cvOneDMaterialAirways::GetArea(double pressure, double z)const{
  // NOTE: o "So_" is the LSA under pressure p1_.
  //         This property comes from the subdomain
  //
  //       o Po is the the zero transmural pressure
  double pres = pressure;
  double So_  = GetS1(z);
  double EHR  = GetEHR(z);
  double u = ((pres - p1_) / EHR) + 1.0;
  if (u <= 0.0) {
    cerr << "ERROR: negative argument in GetArea, u=" << u << endl;
    abort();
  }

  double area = So_*pow(((pres - p1_)/EHR) + 1,0.25);
  return area;
}

double cvOneDMaterialAirways::GetPressure(double S, double z)const{
  double So_ = GetS1(z);
  double EHR = GetEHR(z);
  double press = p1_ + EHR*(pow(S/So_,4) - 1);
  return press;
}

double cvOneDMaterialAirways::GetDpDS(double S, double z)const{
  double EHR = GetEHR(z);
  double So_ = GetS1(z);
  double dpds= 4.0 * EHR * pow(S,3) / pow(So_,4);
  return dpds;
}

double cvOneDMaterialAirways::GetD2pDS2(double S, double z)const{
  double EHR = GetEHR(z);
  double So_ = GetS1(z);
  return 12 * EHR * pow(S,2) / pow(So_,4);
}

double cvOneDMaterialAirways::GetOutflowFunction(double pressure, double z)const{
  return 0.0; // This is not used in our model
}

double cvOneDMaterialAirways::GetDOutflowDp(double pressure, double z)const{
  return 0.0; // Nor is this.
}

//used for viscosity term in matrix outlet flux term
double cvOneDMaterialAirways::GetDD2PDzDS(double area, double z)const{
  double EHR   = GetEHR(z);
  double ro    = Getr1(z);
  double drodz = GetDr1Dz(z);
  double dEHRdro = 4./3.*K1_;
  double derP = drodz*sqrt(PI)*0.5/area/sqrt(area)*(dEHRdro*ro+EHR);
  return derP;
}

double cvOneDMaterialAirways::GetIntegralpD2S(double area, double z)const{
  double EHR   = GetEHR(z);
  double DS1Dz = GetDS1Dz(z);
  double So_   = GetS1(z);
  double DroDz = GetDr1Dz(z);
  double DEHRinvDr = 4./3.*K1_; //should be /(EHR)^2 but included in IntegralpD2S
  double DEHRinvDz = DEHRinvDr*DroDz;
  double IntegralpD2S = 0.8*DS1Dz*EHR*(pow(area/So_,5) - 1) + So_*DEHRinvDz*(0.2*pow(area/So_,5) - area/So_ + 0.8);
  return IntegralpD2S;
}

double cvOneDMaterialAirways::GetIntegralpS(double area, double z)const{
  double EHR   = GetEHR(z);
  double So_   = GetS1(z);
  double IntegralpS = 0.80 * So_ * EHR * (pow(area/So_,5) - 1.0);
  return IntegralpS;
}

// Careful!! this D2p(S,z)Dz first derivative, 2nd variable
double cvOneDMaterialAirways::GetDpDz(double S, double z)const{
  double So_   = GetS1(z);
  double DS1Dz = GetDS1Dz(z);
  double EHR   = GetEHR(z);
  double drodz  = GetDr1Dz(z); // if straight tube ->0.0
  double dEHRdr = 4./3.*K1_; // dfdr
  double dEHRdz = dEHRdr * drodz;
  double dpdz = dEHRdz*(pow(S/So_,4) - 1) - 4*EHR*(pow(S,4)/pow(So_,5))*DS1Dz;
  return dpdz; //careful this is D2P(S,z)Dz
}

double cvOneDMaterialAirways::GetN(double S)const{
  double R    = sqrt(S/PI);
  double T = Period;
  double del  = sqrt(kinematicViscosity*T/2.0/PI); // from Lighthill, olufsen's thesis, pg56
  double altN = -2.0*PI*kinematicViscosity*R/del; // == -1.45838*r

  return N;

}

double cvOneDMaterialAirways::GetLinCompliance(double z)const{
  double EHR   = GetEHR(z);
  double S_o   = GetS1(z);
  double LinCompliance = 2.0*S_o/EHR; // =dSdP at po,S_o
  return LinCompliance;
}

double cvOneDMaterialAirways::GetnonLinCompliance(double area, double z)const{
  double EHR   = GetEHR(z);
  double S_o   = GetS1(z);
  double ratioA= sqrt(area/S_o);
  double nonLinCompliance = 2.0*S_o/EHR*ratioA*ratioA*ratioA;
  return nonLinCompliance;
}

