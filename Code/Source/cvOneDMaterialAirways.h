#ifndef CVONEDMATERIALAIRWAYS_H
#define CVONEDMATERIALAIRWAYS_H

# include "cvOneDMaterial.h"
# include "cvOneDEnums.h"

class cvOneDMaterialAirways:public cvOneDMaterial{

  public:

    cvOneDMaterialAirways();
    ~cvOneDMaterialAirways();
    cvOneDMaterialAirways (const cvOneDMaterialAirways &rhs);
    cvOneDMaterialAirways& operator= ( const cvOneDMaterialAirways &that);

    void   SetAreas_and_length_and_in_out(double S_top, double S_bottom, double z, double Z_in, double z_out);
    void   SetStop(double S){Stop = S;}
    void   SetSbottom(double S){Sbot = S;}
    void   SetLength(double length){len = length;}
    void   SetMaterialType(double*, double);
    double GetProperty( char* what) const;
    double GetArea( double pressure, double z) const;
    double GetPressure( double S, double z) const;
    double GetDpDS( double area, double z) const;
    double GetD2pDS2( double area, double z) const;
	  double GetDD2PDzDS( double area, double z) const; // for viscosity term
    double GetOutflowFunction( double pressure, double z) const;
    double GetDOutflowDp( double pressure, double z) const;
    double GetIntegralpD2S ( double area, double z) const; 
	  double GetIntegralpS ( double area, double z) const;
	  double GetDpDz( double area, double z) const;
    double GetTopArea() const {return Stop;}
    double GetBotArea() const {return Sbot;}
    double GetEHR(double z) const;
    double GetMette2(double area,double z) const;
 	  double GetN(double S) const;

    int IsRigid() {return rigid;}

    double GetLinCompliance(double z) const;
	  double GetnonLinCompliance(double area,double z) const;
	  void SetPeriod(double period);

    void GetParams(double *K1, double *K2, double *K3, double *Pref) const {*K1 = K1_; *K2 = K2_; *K3 = K3_; *Pref=PP1_;};

  private:

    double Stop;
    double Sbot;
    double len;
    double Zin;
    double Zout;
    double lt;
    double chiS;

    double Period;

    double zstar;
    double alpha;
    double pcrit;

    double K1_;
    double K2_;
    double K3_;
    double PP1_;

    int rigid;

    double GetS1( double z) const;
    double Getr1( double z) const;
    double GetDS1Dz( double z) const;
    double GetDr1Dz(double z) const;
    
};

#endif  // CVONEDMATERIALAIRWAYS_H