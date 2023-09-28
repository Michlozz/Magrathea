#ifndef EOSLIST_H_
#define EOSLIST_H_

#include "EOS.h"

extern EOS *Fe_liquid, *Fe_liquid2, *Fe_fcc, *Fe_bcc, *Fe_hcp, *Fe_hcp2, *Fe_Seager, *Fe_hcp3, *Fe_7Si, *Fe_15Si, *Fe_Dummy, *Si_Pv_Shim, *Si_Pv, *Si_PPv, *Si_PPv_Sakai, *Si_PREM, *Si_BM2fit, *Si_Seager, *Si_Dummy, *Si_liquid, *Si_Liquid_Wolf, *Fo, *Wds, *Rwd, *Akm, *Pv_Doro, *PPv_Doro, *Fo_Sotin, *En, *Mw, *Ice_Seager, *Water_ExoPlex, *Water, *Water_sc_dummy, *IceIh, *IceIh_ExoPlex, *IceVI_ExoPlex, *IceVI_Bezacier, *IceVII_ExoPlex, *IceVII_Bezacier, *IceVII, *IceVIIp, *IceVII_FFH2004, *IceVII_FFH2004fit, *IceVII_FFH2004BM, *IceVII_Fei, *IceVII_FFH2004T, *IceX_HS, *IceX, *IceZeng2013FFH, *IceZeng2013FMNR, *Ice_Dummy, *Gas, *Gas_iso, *watervapor, *Gold, *Plat, *Si_QEOS;

double dTdP_Si_Dummy (double P, double T);
// A temperature gradient that equals to the melting curve. Guarantee the temperature won't drop below the melting curve. 
double dTdP_mix(double P, double T);
//double S_Mix(double rho, double T);
//double S_AQUA(double rho, double T);
//double S_QEOS(double rho, double T);
//void print_S(double P, double T, double S, double grad); 

struct Ps {
    double values[21];// Polynomial coefficients.

    //in a form of 
    //f(x,y)=p00 + p10.*x + p01.*y + p20.*x.^2 + p11.*x.*y + p02.*y.^2 + p30.*x.^3 + p21.*x.^2.*y  ... 
     // + p12.*x.*y.^2 + p03.*y.^3 + p40.*x.^4 + p31.*x.^3.*y + p22.*x.^2.*y.^2 ...
     //+ p13.*x.*y.^3 + p04.*y.^4 + p50.*x.^5 + p41.*x.^4.*y + p32.*x.^3.*y.^2 ...
     //+ p23.*x.^2.*y.^3 + p14.*x.*y.^4 + p05.*y.^5;

    // Overload the * operator to multiply a Ps by a scalar (double)
    Ps operator*(double scalar) const {
        Ps result;
        for (int i = 0; i < 21; ++i) {
            result.values[i] = values[i] * scalar;
        }
        return result;
    }

    // Overload the - operator to subtract two Ps
    Ps operator-(const Ps& other) const {
        Ps result;
        for (int i = 0; i < 21; ++i) {
            result.values[i] = values[i] - other.values[i];
        }
        return result;
    }
     Ps operator+(const Ps& other) const {
        Ps result;
        for (int i = 0; i < 21; ++i) {
            result.values[i] = values[i] + other.values[i];
        }
        return result;
    }  
     void print() const {
        for (int i = 0; i < 21; ++i) {
            cout << "values[" << i << "] = " << values[i] << " ";
        }
        cout << endl;
    } 
};


Ps FindPsAqua(double P, double T);
Ps FindPsRock(double Tlin, double Plin);

#endif
