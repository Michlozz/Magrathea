#include "EOSlist.h"

/*
  phasetype is the name of a phase. The comment about the EOS used for the phase should be in the parentheses separated by a space.
0.	eqntype. 8-12 for RTpress style.
1.	V0 in cm^3 / mol.
	For ice, 1 \AA^3 = N_A / (2*10^24) cm^3/mol = 0.3011 cm^3/mol
2.	K0 in GPa
3.	K0p
4.	K0pp in GPa ^-1
5.	mmol in g / mol, or mean molecular weight of gas, or in g / mol for RTpress style
6.	P0 (GPa) the minimum pressure.  The pressure correspond to V0
7.	Theta0 (K), a fitting parameter to Einstein temperature or Debye temperature
8.	gamma0, a fitting parameter of Grueneisen parameter
9.	beta, a fitting parameter of Grueneisen parameter.  In RTpress style, it represents the "m" which stands for the power-law exponent in the thermal deviation term.  Theoretically expected value: 0.6.
10.	gammainf, a fitting parameter of Grueneisen parameter
11.	gamma0p, volume derivative of the Grueneisen parameter
12.	e0 (10^-6 K^-1), electronic contribution to Helmholtz free energy
13.	g, is an electronic analogue of the Grueneisen parameter
14.	n is the number of atoms in the chemical formula of the compound.  Should have n*NA atoms within V.  The n of ideal gas is the number of atoms per molecule for the purpose of adiabatic index.  NOTE: n=2 for collinear molecules e.g. carbon dioxide!  Isothermal atmosphere can be achieved by setting n=0.
15.     Z is the atomic number (number of electron)
16.	T0, the reference temperature for the thermal pressure
17.     alpha0, the zeroth order coefficient of thermal expansion at a reference pressure P0 in 10^-6 K^-1
18.     alpha1, the first order coefficient of thermal expansion at a reference pressure P0 in 10^-6 K^-2
19.	xi, a power law index to describe the pressure effect of the coefficient of thermal expansion 
20.	cp_a in 10^7 erg/g/K Heat capacity per mass at constant pressure
21.	cp_b, fitting coefficient for specific heat capacity, in 10^7 erg/g/K^2
22.	cp_c, cp = cp_a + cp_b*T - cp_c/T^2. cp in 10^7 erg/g/K, cp_c in 10^7 erg*K/g
23.	Debye_approx, whether use Debye approximation or Einstein approximation. Debye approximation is slower but more accurate at temperature lower than Debye/Einstein temperature.  Positive number for Debye, otherwise Einstein.
24.     thermal_type, indicates the thermal type of the phase.  0 indicates no temperature profile available, 1 indicates entropy method, 2 indicates the temperature gradient method.  The only method to set the gradient is using the modify_dTdP function, 3 indicates ideal gas, 4 indicates the EOS is fitted along the isentrope, type 8 indicates RTpress style.
25-32.  at1-at4 & ap1 - ap4


For RTpress style of EOS, also need a _b array. They are fitted polynomial parameters of the thermal coefficients b(V) in erg/mol.  Convert eV/atom to erg/mol need to multiply eV_erg*n*NA. For example, for MgSiO3, 0.9821 eV/atom = 4.824E12 *0.9821 erg/mol = 4.738E12 erg/mol.
*/


// ==========  Iron  ================

// ---------------------------------
// Liquid Iron, Dorogokupets et al. 2017, Scientific Reports.
// DEFAULT

double Fe_liquid_array[][2] = {{0,2}, {1,7.957}, {2,83.7}, {3,5.97}, {5,mFe}, {6, 1E-4}, {7, 263}, {8, 2.033}, {9, 1.168}, {10,0}, {12,198}, {13, 0.884}, {14, 1}, {15,26}, {16, 1811}};

EOS *Fe_liquid = new EOS("Fe liquid (Dorogokupets)", Fe_liquid_array, sizeof(Fe_liquid_array)/2/sizeof(Fe_liquid_array[0][0]));

// -----------------------------------
// Liquid Iron, Anderson & Ahrens, 1994 JGR

double Fe_liquid2_array[][2] = {{0,1}, {1,7.95626}, {2,109.7}, {3,4.66}, {4,-0.043}, {5,mFe}, {14,1}, {15,26}};

EOS *Fe_liquid2 = new EOS("Fe liquid (Anderson)", Fe_liquid2_array, sizeof(Fe_liquid2_array)/2/sizeof(Fe_liquid2_array[0][0]));

// -----------------------------------
// Alpha Iron (bcc), Dorogokupets et al. 2017, Scientific Reports 

double Fe_bcc_array[][2] = {{0,0}, {1,7.092}, {2,164.0}, {3,5.50}, {5,mFe}, {8,303}, {9,1.736}, {10,1.125}, {11,0}, {13,198}, {14,1.0}, {17,1043}};

EOS *Fe_bcc = new EOS("Fe bcc (Dorogokupets)", Fe_bcc_array, sizeof(Fe_bcc_array)/2/sizeof(Fe_bcc_array[0][0]));

// -----------------------------------
// Gamma Iron (fcc), Dorogokupets et al. 2017, Scientific Reports

double Fe_fcc_array[][2] = {{0,0}, {1,6.9285}, {2,146.2}, {3,4.67}, {5,mFe}, {8,222.5}, {9,2.203}, {10,0.01}, {11,0}, {13,198}, {14,0.5}};

EOS *Fe_fcc = new EOS("Fe fcc (Dorogokupets)", Fe_fcc_array, sizeof(Fe_fcc_array)/2/sizeof(Fe_fcc_array[0][0]));

// -----------------------------------
// Epsilon Iron (hcp), Smith et al. 2018, Nature Astronomy. (Gruneisen determined from fitting Fig. 3b)
// DEFAULT

double Fe_hcp_array[][2] = {{0,2}, {1,mFe/8.43}, {2,177.7}, {3,5.64}, {5,mFe}, {7,322}, {8,2.09}, {9,1.01}, {10,0.0500}, {14,1}, {15,26}};

EOS *Fe_hcp = new EOS("Fe hcp (Smith)", Fe_hcp_array, sizeof(Fe_hcp_array)/2/sizeof(Fe_hcp_array[0][0]));

// -----------------------------------
// Epsilon Iron (hcp), Bouchet et al. 2013, PRB 87, 094102

double Fe_hcp2_array[][2] = {{0,3}, {1,6.29}, {2,253.844}, {3,4.719}, {5,mFe}, {7,44.574}, {8,1.408}, {9,0.826}, {10,0.827}, {12,212.1}, {13,1.891}, {14,1}, {15,26}};

EOS *Fe_hcp2 = new EOS("Fe hcp (Bouchet)", Fe_hcp2_array, sizeof(Fe_hcp2_array)/2/sizeof(Fe_hcp2_array[0][0]));

// -----------------------------------
// Epsilon Iron (hcp), Dorogokupets et al. 2017, Scientific Reports.

double Fe_hcp3_array[][2] = {{0,2}, {1,6.8175}, {2,148.0}, {3,5.86}, {5,mFe}, {7, 227}, {8, 2.2}, {9, 0.01}, {10,0}, {12,126}, {13,-0.83}, {14,1}, {15,26}, {16, 298.15}};

EOS *Fe_hcp3 = new EOS("Fe hcp (Dorogokupets)", Fe_hcp3_array, sizeof(Fe_hcp3_array)/2/sizeof(Fe_hcp3_array[0][0]));

// -----------------------------------
// Iron-Silicate Alloy, 7 wt% silicate, Wicks et al. 2018, Science Advances.

double Fe_7Si_array[][2] = {{0,2}, {1,7.02}, {2,136.2}, {3,5.97}, {5,mFe*0.93+mSi*0.07}};

EOS *Fe_7Si = new EOS("Fe-7Si (Wicks)", Fe_7Si_array, sizeof(Fe_7Si_array)/2/sizeof(Fe_7Si_array[0][0]));

// -----------------------------------
// Iron-Silicate Alloy, 15 wt% silicate, Wicks et al. 2018, Science Advances.

double Fe_15Si_array[][2] = {{0,2}, {1,6.784}, {2,227.9}, {3,4.74}, {5,mFe*0.85+mSi*0.15}};

EOS *Fe_15Si = new EOS("Fe-15Si (Wicks)", Fe_15Si_array, sizeof(Fe_15Si_array)/2/sizeof(Fe_15Si_array[0][0]));

// -----------------------------------
// Iron, Seager et al. 2007 ApJ 669:1279, tabulate EOS
EOS *Fe_Seager = new EOS("Fe (Seager)", "./tabulated/iron.txt");

// -----------------------------------
// Iron Dummy, Used to fill in phase space that no EOS provided.

EOS *Fe_Dummy = new EOS("Fe Dummy", Fe_hcp3_array, sizeof(Fe_hcp3_array)/2/sizeof(Fe_hcp3_array[0][0]));

// ==========  Silicate  ================

// ---------------------------------
// Liquid Magnesium Silicate, MgSiO3, Mosenfelder et al. 2009, Journal of Geophysical Research: Solid Earth, Table 5, BM4LC38.202E-5 m^3/kg = 64.1 AA^3

double Si_liquid_array[][2] = {{0,1}, {1,64.1}, {2,24.7}, {3,9.2}, {4,-1.87}, {5,mMg+mSi+3*mO}, {8,0.37}, {9,-1.71}, {10, 0}, {14,5}, {16, 1673}, {24, 4}};

EOS *Si_liquid = new EOS("Si liquid (Mosenfelder)", Si_liquid_array, sizeof(Si_liquid_array)/2/sizeof(Si_liquid_array[0][0]));

// -----------------------------------
// Liquid Magnesium Silicate, MgSiO3, Wolf & Bower 2018, Table 1 S11, Using RTpress structure
// DEFAULT

double Si_Liquid_Wolf_array[][2] = {{0, 10}, {1, 38.99}, {2, 13.2}, {3, 8.238}, {5, mMg+mSi+3*mO}, {8, 0.1899}, {9, 0.6}, {11, -1.94}, {14, 5},  {16, 3000}};
double Si_Liquid_Wolf_b[] = {4.738E12, 2.97E12, 6.32E12, -1.4E13, -2.0E13};

EOS *Si_Liquid_Wolf = new EOS("Si liquid (Wolf)", Si_Liquid_Wolf_array, Si_Liquid_Wolf_b, sizeof(Si_Liquid_Wolf_array)/2/sizeof(Si_Liquid_Wolf_array[0][0]), sizeof(Si_Liquid_Wolf_b)/sizeof(Si_Liquid_Wolf_b[0]));

//----------------------------------------
// Forsterite/Olivine, Mg2Si04, Dorogokupets et al. 2015, Russ. Geol. Geophys.
double Fo_array[][2] = {{0,0}, {1,43.67}, {2,127.4}, {3,4.3}, {5,2*mMg+mSi+4*mO}, {7,949}, {8,1.066}, {9,2.225}, {10,0.0}, {14,7}};

EOS *Fo = new EOS("Fo/Ol (Dorogokupets)", Fo_array, sizeof(Fo_array)/2/sizeof(Fo_array[0][0]));

//----------------------------------------
// Wadsleyite, Mg2Si04, Dorogokupets et al. 2015, Russ. Geol. Geophys.

double Wds_array[][2] = {{0,0}, {1,40.54}, {2,169.0}, {3,4.14}, {5,2*mMg+mSi+4*mO}, {7,921}, {8,1.185}, {9,2.10}, {10,0.0}, {14,7}};

EOS *Wds = new EOS("Wds (Dorogokupets)", Wds_array, sizeof(Wds_array)/2/sizeof(Wds_array[0][0]));

//----------------------------------------
//Ringwoodite, Mg2Si04, Dorogokupets et al. 2015, Russ. Geol. Geophys.

double Rwd_array[][2] = {{0,0}, {1,39.5}, {2,187.4}, {3,3.98}, {5,2*mMg+mSi+4*mO}, {7,929}, {8,1.21}, {9,1.35}, {10,0.0}, {14,7}};

EOS *Rwd = new EOS("Rwd (Dorogokupets)", Rwd_array, sizeof(Rwd_array)/2/sizeof(Rwd_array[0][0]));

//----------------------------------------
// Akimotoite, MgSi03, Dorogokupets et al. 2015, Russ. Geol. Geophys.

double Akm_array[][2] = {{0,0}, {1,26.35}, {2,215.3}, {3,4.91}, {5,108.27}, {7,995}, {8,2.00}, {9,1.41}, {10,0.0}, {14,5}};

EOS *Akm = new EOS("Akm (Dorogokupets et al.)", Akm_array, sizeof(Akm_array)/2/sizeof(Akm_array[0][0]));

//--------------------------------------
// Forsterite/Olivine, Mg2SiO4, Sotin et al. 2007, Icarus, Using Duffy et al. 1995 & Bouhifd et al. 1996

double Fo_Sotin_array[][2] = {{0,0}, {1,(2*mMg+mSi+4*mO)/3.22}, {2,128}, {3,4.3}, {5,2*mMg+mSi+4*mO}, {14,7}, {16, 300}, {17,28.32}, {18, 0.00758}, {20,0.840}};

EOS *Fo_Sotin = new EOS("Fo/Ol (Sotin)", Fo_Sotin_array, sizeof(Fo_Sotin_array)/2/sizeof(Fo_Sotin_array[0][0]));

//---------------------------------------
// Enstatite/Orthopyroxene, MgSiO3, Sotin et al. 2007, Icarus, Using Vacher et al. 1998 & Anderson et al. 1991

double En_array[][2] = {{0,0}, {1,(mMg+mSi+3*mO)/3.215}, {2,111}, {3,7}, {5,mMg+mSi+3*mO}, {14,5}, {16, 300}, {17,28.6}, {18, 0.0072}, {20,0.840}};

EOS *En = new EOS("En/Opx (Sotin)", En_array, sizeof(En_array)/2/sizeof(En_array[0][0]));

//--------------------------------------
// Magnesiowustite, MgO, Sotin et al. 2007, Icarus,

double Mw_array[][2] = {{0,0}, {1,(mMg+mO)/3.584}, {2,157}, {3,4.4}, {5,mMg+mO}, {7,430}, {8, 1.45}, {9,3}, {14,2}};

EOS *Mw = new EOS("Magnesiowustite (Sotin)", Mw_array, sizeof(Mw_array)/2/sizeof(Mw_array[0][0]));

// ---------------------------------
// Bridgmanite/Perovskite, MgSiO3, Oganov & Ono 2004, Nature, GGA
// DEFAULT

double Si_Pv_array[][2] = {{0,2}, {1,25.206}, {2,230.05}, {3,4.142}, {5,mMg+mSi+3*mO}, {6, -11.2}, {7, 1000}, {8,1.506}, {9,7.02469}, {10,1.14821}, {14,5}};

EOS *Si_Pv = new EOS("Brg (Oganov)", Si_Pv_array, sizeof(Si_Pv_array)/2/sizeof(Si_Pv_array[0][0]));

// ---------------------------------
// Bridgmanite/Perovskite, MgSiO3, Shim & Duffy 2000, American Mineralogist

double Si_Pv_Shim_array[][2] = {{0,0}, {1,24.43}, {2,261}, {3,4}, {5,mMg+mSi+3*mO}, {7,1000}, {8,1.42}, {9,2}, {14,5}};

EOS *Si_Pv_Shim = new EOS("Brg (Shim)", Si_Pv_Shim_array, sizeof(Si_Pv_Shim_array)/2/sizeof(Si_Pv_Shim_array[0][0]));

//----------------------------------------
// Bridgmanite/Perovskite, MgSi03, Dorogokupets et al. 2015, Russ. Geol. Geophys.

double Pv_Doro_array[][2] = {{0,0}, {1,24.45}, {2,252.0}, {3,4.38}, {5, mMg+mSi+3*mO}, {7,943}, {8,1.70}, {9,3.00}, {10,0.0}, {14,5}};

EOS *Pv_Doro = new EOS("Pv (Dorogokupets)", Pv_Doro_array, sizeof(Pv_Doro_array)/2/sizeof(Pv_Doro_array[0][0]));

// ---------------------------------
// Post-Perovskite, MgSiO3, Sakai, Dekura, & Hirao, 2016, Scientific Reports
// DEFAULT

double Si_PPv_Sakai_array[][2] = {{0,4}, {1,24.73}, {2,203}, {3,5.35}, {5,mMg+mSi+3*mO}, {7,848}, {8,1.47}, {9,2.7}, {10,0.93}, {14,5}};

EOS *Si_PPv_Sakai = new EOS("Si PPv (Sakai)", Si_PPv_Sakai_array, sizeof(Si_PPv_Sakai_array)/2/sizeof(Si_PPv_Sakai_array[0][0]));

// ---------------------------------
// Post-Perovskite, MgSiO3, Oganov & Ono 2004, Nature, GGA

double Si_PPv_array[][2] = {{0,2}, {1,25.239}, {2,199.96}, {3,4.541}, {5,mMg+mSi+3*mO}, {7, 1500}, {8,1.553}, {9,4.731}, {10,1.114}, {14,5}};

EOS *Si_PPv = new EOS("Si PPv (Oganov)", Si_PPv_array, sizeof(Si_PPv_array)/2/sizeof(Si_PPv_array[0][0]));

//----------------------------------------
// Postperovskite, MgSi03, Dorogokupets et al. 2015, Russ. Geol. Geophys.

double PPv_Doro_array[][2] = {{0,0}, {1,24.2}, {2,253.7}, {3,4.03}, {5,mMg+mSi+3*mO}, {7,943}, {8,1.67}, {9,2.22}, {10,0.0}, {14,5}};

EOS *PPv_Doro = new EOS("PPv (Dorogokupets)", PPv_Doro_array, sizeof(PPv_Doro_array)/2/sizeof(PPv_Doro_array[0][0]));

// -----------------------------------
// Silicate PREM mantle EOS in Appendix F.1 of Stacey & Davis 2008, used in Zeng 2016
EOS *Si_PREM = new EOS("Si (PREM)", "./tabulated/SiPREM.txt");

// -----------------------------------
// Silicate PREM BM2 extrapolation used in Zeng 2016
double Si_BM2fit_array[][2] = {{0,0}, {1,25.223}, {2,206},{3,4},{5,mMg+mSi+3*mO}};
EOS *Si_BM2fit = new EOS("Si (PREM, Zeng)", Si_BM2fit_array, sizeof(Si_BM2fit_array)/2/sizeof(Si_BM2fit_array[0][0]));

// -----------------------------------

double Si_QEOS_array[][2] = {{0,5}, {1, 27.36}, {5, mSi+2*mO}};
EOS *Si_QEOS = new EOS("Si (QEOS)", Si_QEOS_array, sizeof(Si_QEOS_array)/2/sizeof(Si_QEOS_array[0][0]));


// -----------------------------------
// Silicate, Seager et al. 2007 ApJ 669:1279, tabulate EOS
EOS *Si_Seager = new EOS("Si (Seager)", "./tabulated/silicate.txt");

// ---------------------------------
// Si Dummy, Used to fill in phase space that no EOS provided.

EOS *Si_Dummy = new EOS("Si Dummy", Si_BM2fit_array, sizeof(Si_BM2fit_array)/2/sizeof(Si_BM2fit_array[0][0]));

double dTdP_Si_Dummy (double P, double T)
// A temperature gradient that equals to the melting curve. Guarantee the temperature won't drop below the melting curve. 
{
  P /= 1E10;
  if (T > 1830*pow(1+P/4.6, 0.33))
    return 131.2826087*pow(1+P/4.6, -0.67)/1E10;
  else
  {
    cout<<"Error: The pressure "<<P<<" GPa and temperature "<<T<<" K are inconsistent with liquid silicate."<<endl;
    return numeric_limits<double>::quiet_NaN();
  }
}

// ==========  Water  ================

// -----------------------------------
// Liquid water, ExoPlex, unkown source
double Water_ExoPlex_array[][2] = {{0,1}, {1,18.797}, {2,2.06}, {3,6.29}, {4,-0.9175}, {5,18.01528}, {20,4.184}};

EOS *Water_ExoPlex = new EOS("Water (ExoPlex)", Water_ExoPlex_array, sizeof(Water_ExoPlex_array)/2/sizeof(Water_ExoPlex_array[0][0]));

// -----------------------------------
// Liquid water, Valencia et al. 2007, ApJ, 656:545
// DEFAULT
double Water_array[][2] = {{0,0}, {1,18.047}, {2,2.18}, {5,18.01528}};

EOS *Water = new EOS("Water (Valencia)", Water_array, sizeof(Water_array)/2/sizeof(Water_array[0][0]));

// -----------------------------------
// Dummy for supercritical water. 
// DEFAULT
EOS *Water_sc_dummy = new EOS("Water supercritical Dummy", Water_array, sizeof(Water_array)/2/sizeof(Water_array[0][0]));

// -----------------------------------
// Ice Ih, Feistel & Wagner 2006, Acuna et al. 2021
// DEFAULT

double IceIh_array[][2] = {{0,0}, {1,19.56}, {2,9.5}, {3,5.3}, {5,18.01528}, {20,1.913}};

EOS *IceIh = new EOS("Ice Ih", IceIh_array, sizeof(IceIh_array)/2/sizeof(IceIh_array[0][0]));

// -----------------------------------
// Ice Ih, ExoPlex, unkown source

double IceIh_ExoPlex_array[][2] = {{0,0}, {1,19.65}, {2,9.2}, {3,5.5}, {5,18.01528}, {20,4.184}};

EOS *IceIh_ExoPlex = new EOS("Ice Ih (ExoPlex)", IceIh_ExoPlex_array, sizeof(IceIh_ExoPlex_array)/2/sizeof(IceIh_ExoPlex_array[0][0]));

// -----------------------------------
// Ice VI, ExoPlex, Bezacier et al. 2014

double IceVI_ExoPlex_array[][2] = {{0,0}, {1,14.17}, {2,14.01}, {3,4}, {5,18.01528}};

EOS *IceVI_ExoPlex = new EOS("Ice VI (ExoPlex)", IceVI_ExoPlex_array, sizeof(IceVI_ExoPlex_array)/2/sizeof(IceVI_ExoPlex_array[0][0]));

// -----------------------------------
// Ice VI, Bezacier et al. 2014 & Tchijov et al. 2004
// DEFAULT

double IceVI_Bezacier_array[][2] = {{0,0}, {1,14.17}, {2,14.05}, {3,4}, {5,18.01528}, {16, 300}, {17, 146.}, {20, 2.6}};

EOS *IceVI_Bezacier = new EOS("Ice VI (Bezacier)", IceVI_Bezacier_array, sizeof(IceVI_Bezacier_array)/2/sizeof(IceVI_Bezacier_array[0][0]));

// -----------------------------------
// Ice VII, Bezacier et al. 2014 & Tchijov et al. 2004
// DEFAULT

double IceVII_Bezacier_array[][2] = {{0,0}, {1,12.49}, {2,20.15}, {3,4}, {5,18.01528}, {16, 300}, {17, 115.8}, {20, 2.3}};

EOS *IceVII_Bezacier = new EOS("Ice VII (Bezacier)", IceVII_Bezacier_array, sizeof(IceVII_Bezacier_array)/2/sizeof(IceVII_Bezacier_array[0][0]));

// -----------------------------------
// Ice VII Isothermal, Bezacier et al. 2014

double IceVII_ExoPlex_array[][2] = {{0,0}, {1,12.49}, {2,20.15}, {3,4}, {5,18.01528}};

EOS *IceVII_ExoPlex = new EOS("Ice VII (ExoPlex)", IceVII_ExoPlex_array, sizeof(IceVII_ExoPlex_array)/2/sizeof(IceVII_ExoPlex_array[0][0]));

// -----------------------------------
// Ice VII, Zachary Grande

double IceVII_array[][2] = {{0,2}, {1,12.80}, {2,18.47}, {3,2.51}, {5,18.01528}};

EOS *IceVII = new EOS("Ice VII (Grande)", IceVII_array, sizeof(IceVII_array)/2/sizeof(IceVII_array[0][0]));

// -----------------------------------
// Ice VII', Zachary Grande

double IceVIIp_array[][2] = {{0,2}, {1,12.38}, {2,20.76}, {3,4.49}, {5,18.01528}};

EOS *IceVIIp = new EOS("Ice VII' (Grande)", IceVIIp_array, sizeof(IceVIIp_array)/2/sizeof(IceVIIp_array[0][0]));

// -----------------------------------
// Ice VII, Frank, Fei & Hu 2004 300K, apply 3rd Birch Murnaghan fit result in the paper to Vinet EOS.

double IceVII_FFH2004_array[][2] = {{0,2}, {1,12.4}, {2,21.1}, {3,4.4}, {5,18.01528}};

EOS *IceVII_FFH2004 = new EOS("Ice VII (FFH2004, Vinet)", IceVII_FFH2004_array, sizeof(IceVII_FFH2004_array)/2/sizeof(IceVII_FFH2004_array[0][0]));

// -----------------------------------
// Ice VII, Frank, Fei & Hu 2004 300K, fitting with 3rd order Vinet

double IceVII_FFH2004fit_array[][2] = {{0,2}, {1,12.42}, {2,19.84}, {3,4.99}, {5,18.01528}};

EOS *IceVII_FFH2004fit = new EOS("Ice VII (FFH2004fit, Vinet fit)", IceVII_FFH2004fit_array, sizeof(IceVII_FFH2004fit_array)/2/sizeof(IceVII_FFH2004fit_array[0][0]));


// -----------------------------------
// Ice VII, Frank, Fei & Hu 2004 300K, fitting with 3rd order Birch Murnaghan.

double IceVII_FFH2004BM_array[][2] = {{0,0}, {1,12.4}, {2,21.1}, {3,4.4}, {5,18.01528}};

EOS *IceVII_FFH2004BM = new EOS("Ice VII (FFH2004, BM)", IceVII_FFH2004BM_array, sizeof(IceVII_FFH2004BM_array)/2/sizeof(IceVII_FFH2004BM_array[0][0]));

// -----------------------------------
// Ice VII, Frank, Fei & Hu 2004 original with temperature using thermal expansion representation and the heat capacity from Myint et al. 2017.

double IceVII_FFH2004T_array[][2] = {{0,0}, {1,12.4}, {2,21.1}, {3,4.4}, {5,18.01528}, {16,300}, {17, -420}, {18, 1.56}, {19, 1.1}, {20, 0}, {21, 4E-3}, {24,9}};

EOS *IceVII_FFH2004T = new EOS("Ice VII (FFH2004, thermal)", IceVII_FFH2004T_array, sizeof(IceVII_FFH2004T_array)/2/sizeof(IceVII_FFH2004T_array[0][0]));

// -----------------------------------
// Ice VII, Fei et al. 1993 modified by Sotin et al. 2007

double IceVII_Fei_array[][2] = {{0,0}, {1,12.3}, {2,23.9}, {3,4.2}, {5,18.01528}, {7, 1470}, {8, 1.2}, {9, 1}, {10, 0}, {14,3}};

EOS *IceVII_Fei = new EOS("Ice VII (Fei)", IceVII_Fei_array, sizeof(IceVII_Fei_array)/2/sizeof(IceVII_Fei_array[0][0]));

// -----------------------------------
// Ice X, Zachary Grande
// DEFAULT

double IceX_array[][2] = {{0,2}, {1,10.18}, {2,50.52}, {3,4.5}, {5,18.01528}};

EOS *IceX = new EOS("Ice X (Grande)", IceX_array, sizeof(IceX_array)/2/sizeof(IceX_array[0][0]));

// -----------------------------------
// Ice X, Hermann & Schwerdtfeger 2011, PRL 106, 187403

double IceX_HS_array[][2] = {{0,0}, {1,7.588}, {2,194.73}, {3,4}, {5,18.01528}};

EOS *IceX_HS = new EOS("Ice X (Hermann)", IceX_HS_array, sizeof(IceX_HS_array)/2/sizeof(IceX_HS_array[0][0]));

// -----------------------------------
// Ice, Seager et al. 2007 ApJ 669:1279, tabulate EOS
EOS *Ice_Seager = new EOS("Ice (Seager)", "./tabulated/water.txt");

// -----------------------------------
// Ice Dummy,  Used to fill in phase space that no EOS provided.

EOS *Ice_Dummy = new EOS("Ice Dummy", IceVI_ExoPlex_array, sizeof(IceVI_ExoPlex_array)/2/sizeof(IceVI_ExoPlex_array[0][0]));

// -----------------------------------
// Modified EOSs to match Zeng 2013
double FMNR[2][13] = {{2.45, 2.5, 3.25, 3.5, 3.75, 4, 5, 6, 7, 9, 11, 13, 15}, {37.3, 43.6, 140, 183, 234, 290, 587, 966, 1440, 2703, 4405, 6416, 8893}};

double Zeng2013FFH(double P, double T)
{
  P/=1E10;			// convert to GPa

  return 18.01528 * (0.0805 + 0.0229*(1-exp(-0.0743*P)) + 0.1573*(1-exp(-0.0061*P)));
}
  
EOS *IceZeng2013FFH = new EOS("Ice (FFH 2004)", Zeng2013FFH);

EOS *IceZeng2013FMNR = new EOS("Ice (FMNR 2009)", FMNR[1], FMNR[0], 13);

// =========  Atmosphere  ================

// -----------------------------------
// Adiabtic Ideal Gas, Parameter 5 can be changed for mean molecular weight of gas. 3 g/mol = mix of H/He
// DEFAULT

double Gas_array[3][2] = {{0,6}, {5,3}, {14,2}};

EOS *Gas = new EOS("Ideal Gas", Gas_array, 3);

// -----------------------------------
// Isothermal Ideal Gas
// DEFAULT

double Gas_iso_array[3][2] = {{0,6}, {5,3}, {14,0}};

EOS *Gas_iso = new EOS("Isothermal Ideal Gas", Gas_iso_array, 3);

// -----------------------------------
// Water Vapor Ideal Gas, same as adiabatic but 5 and 14 changed for water

double watervapor_array[3][2] = {{0,6}, {5,18}, {14,3}};

EOS *watervapor = new EOS("Water vapor", watervapor_array, 3);

// ==========  OTHER  ================

// -----------------------------------
// Gold, Heinz & Jeanloz 1983, J. Appl. Phys. (included for a Hitchiker's-related joke)

double Gold_array[][2] = {{0,0}, {1,10.215}, {2,166.65}, {3,5.4823}, {5,196.96657}, {7,170}, {8,2.95}, {11,1.7}, {14,1}};

EOS *Gold = new EOS("Gold", Gold_array, sizeof(Gold_array)/2/sizeof(Gold_array[0][0]));

// -----------------------------------
// Platinum, Matsui et al. 2009, J. Appl. Phys. (included for a Hitchiker's-related joke)

double Plat_array[][2] = {{0,2}, {1,10.03}, {2,273}, {3,5.20}, {5,195.084}, {7,230}, {8,2.7}, {9,1.1}, {14,1}};

EOS *Plat = new EOS("Plat", Plat_array, sizeof(Plat_array)/2/sizeof(Plat_array[0][0]));


// ============== An example on the format of dTdP function ==============
double dTdP_gas(double P, double T)
// return the adiabatic temperature gradient at given pressure and temperature point for ideal gas.
{
  if (P != 0)
    return 2.*T / (7.*P);
  else
  {
    cout<<"Error: Can't get adiabatic temperature gradient for diatomic gas at P=0."<<endl;
    return numeric_limits<double>::quiet_NaN();
  }
}

double dTdP_QEOS(double P, double T)
// return the adiabatic temperature gradient at given pressure and temperature point for ideal gas.
{
  if (P != 0){
    //cout <<"dTdP_S="<<dTdP_S(P,T);
    //cin.ignore();
    return 2.*T / (7.*P);
    }
  else
  {
    cout<<"Error: Can't get adiabatic temperature gradient at P=0."<<endl;
    return numeric_limits<double>::quiet_NaN();
  }
}

double S_Mix(double rho, double T)
{
double Xr=0.5; //mass function of rock
double  n=3.0;
//double 
// double tmp1=S_QEOS(rho, T);
// double tmp2=S_AQUA(rho, T);

// cout<<"im in S_mix"<<endl;
// cout<<"T="<<T<<endl;
// cout<<"rho="<<rho<<endl;
// cout<<"S_QEOS="<<tmp1<<endl;
// cout<<"S_AQUA(1.5 , 700)="<<S_AQUA(1.5 , 700)<<endl;
// cout<<"S_QEOS(1.5 , 700)="<<S_QEOS(1.5 , 700)<<endl;

  return (Xr*S_QEOS(rho, T) +(1-Xr)*S_AQUA(rho, T) )/(n*R);

}
double S_QEOS(double rho, double T){
    //rock, based on QEOS
    T= log10(T/(1.1605e+7)); // to logKev 
    rho=log10(rho); //to log scale

    T= round(T * 10000.0) / 10000.0; //round a bit
     rho= round(rho * 10000.0) / 10000.0; //round a bit

    double s=0;
    double x = T;
    double y = rho;

    double p00 = 0;
    double p10 = 0;
    double p01 = 0;
    double p20 = 0;
    double p11 = 0;
    double p02 = 0;
    double p30 = 0;
    double p21 = 0;
    double p12 = 0;
    double p03 = 0;
    double p40 = 0;
    double p31 = 0;
    double p22 = 0;
    double p13 = 0;
    double p04 = 0;
    double p50 = 0;
    double p41 = 0;
    double p32 = 0;
    double p23 = 0;
    double p14 = 0;
    double p05 = 0;

    if (T < -4.2 && rho <= -10) { // reg 1
        p00 = 0.004069;
        p10 = 0.002223;
        p01 = -0.004093;
        p20 = 0.001253;
        p11 = -0.0008889;
        p02 = 0.004368;
        p21 = -0.000793;
        p12 = 0.001997;
        p03 = -0.005752;
        p22 = 0.0005726;
        p13 = -0.001047;
        p04 = 0.002362;

        x = (x + 5.15) / 0.5193;
        y = (y + 11.3) / 0.7797;
    }

    if (T < -4.2 && rho > -10 && rho <= -1.3) { // reg 2
        p00 = -0.5919;
        p10 = -0.5767;
        p01 = -0.0001139;
        p20 = -0.1951;
        p11 = -1.754e-05;
        p02 = -2.455e-05;
        p30 = -0.028;
        p21 = -3.813e-06;
        p12 = 6.834e-08;
        p03 = -2.544e-06;
        p40 = -0.001459;
        p31 = -5.085e-07;
        p22 = 4.466e-07;
        p13 = -3.179e-07;
    }

    if (T < -3.4 && rho > -1.3 && rho < 0.4) { // reg 3
        p00 = 3.583;
        p10 = 3.943;
        p01 = -1.223;
        p20 = 1.744;
        p11 = -1.059;
        p02 = 0.2075;
        p30 = 0.3843;
        p21 = -0.3403;
        p12 = 0.1345;
        p03 = -0.02197;
        p40 = 0.042;
        p31 = -0.04812;
        p22 = 0.02866;
        p13 = -0.009447;
        p04 = 0.001556;
        p50 = 0.001817;
        p41 = -0.002528;
        p32 = 0.00201;
        p23 = -0.0009945;
        p14 = 0.000328;
        p05 = -6.824e-05;
    }

    if (T < -3.4 && rho >= 0.4) { // reg 4
        p00 = -0.7245;
        p10 = -0.665;
        p01 = 1.163;
        p20 = -0.2116;
        p11 = 0.8952;
        p02 = -0.3606;
        p30 = -0.0286;
        p21 = 0.2429;
        p12 = -0.2175;
        p03 = 0.03663;
        p40 = -0.001406;
        p31 = 0.02806;
        p22 = -0.04054;
        p13 = 0.01742;
        p04 = 0.0003217;
        p41 = 0.001177;
        p32 = -0.002402;
        p23 = 0.001771;
        p14 = -0.0002133;
        p05 = -0.0001547;


        x = (x + 3.975) / 0.1941;
        y = (y + 6.632) / 2.954;
    }

    if (T >= -3.4) { // reg 5
        if (rho <= 6.0833 * T + 8.0833) { // reg 5a
            p00 = 0.1124;
            p10 = -0.6784;
            p01 = -0.08861;
            p20 = -0.6971;
            p11 = 0.05124;
            p02 = -0.0005395;
            p30 = -0.1636;
            p21 = 0.02927;
            p12 = -0.0009521;
        } else { // reg 5b
            p00 = 0.6412;
            p10 = 0.507;
            p01 = -0.2032;
            p20 = 0.1435;
            p11 = -0.1148;
            p02 = 0.007497;
            p30 = 0.01415;
            p21 = -0.01705;
            p12 = 0.002012;
        }
    }

    if (T < -3.4 && T >= -4.2 && rho <= -14.697 * pow(T, 2) - 99.364 * T - 169.27) { // reg 6
        p00 = -0.08332;
        p10 = -0.06331;
        p01 = -0.1547;
        p20 = -0.008806;
        p11 = -0.07599;
        p02 = 0.003598;
        p21 = -0.009395;
        p12 = 0.0006994;
        p03 = 1.589e-05;
    }

    if (T < -3.4 && T >= -4.2 && rho > -14.697 * pow(T, 2) - 99.364 * T - 169.27 && rho <= -14.167 * pow(T, 2) - 96.917 * T - 167.09) { // reg 7
        p00 = 0.0866;
        p10 = 0.4228;
        p01 = -0.4032;
        p20 = 0.6716;
        p11 = -1.414;
        p02 = 0.6675;
        p30 = -0.0273;
        p21 = -0.2214;
        p12 = 0.2476;
        p03 = -0.0004677;
        p40 = 0.3077;
        p31 = -1.131;
        p22 = 1.496;
        p13 = -0.7515;
        p04 = 0.101;
        p50 = 1.268;
        p41 = -7.067;
        p32 = 15.55;
        p23 = -16.92;
        p14 = 9.141;
        p05 = -1.972;

         x = (x + 3.975) / 0.1941;
        y = (y + 6.632) / 2.954;
    }
if (T < -3.4 && T >= -4.2 && rho < -1.3 && rho > -14.167 * pow(T, 2) - 96.917 * T - 167.09) { // reg 8
    p00 = 0.01001;
    p10 = 0.001134;
    p01 = 0.0008743;
    p20 = 0.0002612;
    p11 = -3.232e-05;
    p02 = 0.0001316;
    p30 = 0.004425;
    p21 = -0.0117;
    p12 = 0.009942;
    p03 = -0.002865;
    p40 = 0.00334;
    p31 = -0.01372;
    p22 = 0.01893;
    p13 = -0.01084;
    p04 = 0.002044;
    p50 = 0.0001318;
    p41 = -0.003017;
    p32 = 0.008457;
    p23 = -0.007669;
    p14 = 0.002452;
    p05 = -0.000202;

    x = (x + 4.042) / 0.1581;
    y = (y + 4.243) / 2.211;
}

    s = p00 + p10 * x + p01 * y + p20 * pow(x, 2) + p11 * x * y + p02 * pow(y, 2) + p30 * pow(x, 3) + p21 * pow(x, 2) * y
         + p12 * x * pow(y, 2) + p03 * pow(y, 3) + p40 * pow(x, 4) + p31 * pow(x, 3) * y + p22 * pow(x, 2) * pow(y, 2)
         + p13 * x * pow(y, 3) + p04 * pow(y, 4) + p50 * pow(x, 5) + p41 * pow(x, 4) * y + p32 * pow(x, 3) * pow(y, 2)
         + p23 * pow(x, 2) * pow(y, 3) + p14 * x * pow(y, 4) + p05 * pow(y, 5);


// s in units of (jrk/(g*keV))
    if (s <= 0) {
      s = 0;
    }

 s = s * 8.61698 * pow(10, 8); // to erg/(g*K)

    return s;
}

double S_AQUA(double rho, double T) {
    // AQUA form Jonnas

    // if (rho>=7.935){
    //     cout<<"rho="<<rho <<endl;
    //     cout <<"T="<<T << endl;
    //     cin.ignore();
    // }
    T= log10(T); //to log scale
    rho=log10(1000*rho); //form gr/cm3 to log10(kg/m3)
    
    T = round(T * 10000) / 10000; // Small rounding here
    rho = round(rho * 10000) / 10000;// Small rounding here

    //cout<<"inside S_AQU
    double x = T;
    double y = rho;
    double s = 1e-16;//zero

    if (rho > 3.90) {
        cout<<"Warning! in S_AQUA"<< endl;
        cout << "rho outside the range" << endl;
        cout<<"rho="<<rho << " in log10(kg/m3)"<< endl;
        cout <<"T="<<T <<" in log10(K)"<< endl;
        cout<<"setting rho=3.90 in log10(kg/m3) for s calcuation"<<endl;
        //cout<<"returning s=~0"<<endl;
      //  return s;
    }

    double p00 = 0;
    double p10 = 0;
    double p01 = 0;
    double p20 = 0;
    double p11 = 0;
    double p02 = 0;
    double p30 = 0;
    double p21 = 0;
    double p12 = 0;
    double p03 = 0;
    double p40 = 0;
    double p31 = 0;
    double p22 = 0;
    double p13 = 0;
    double p04 = 0;
    double p50 = 0;
    double p41 = 0;
    double p32 = 0;
    double p23 = 0;
    double p14 = 0;
    double p05 = 0;

    if (T <= 2.43 && rho <= 2.96 && rho > 2.01) { // reg 0a
        p00 = 1.297e+04;
        p10 = -2991;
        p01 = -3573;
        p20 = 502.2;
        p11 = 1029;
        p02 = -1143;
        p21 = -161.4;
        p12 = 359.2;
        p03 = -272;
        p22 = -49.78;
        p13 = 74.98;
        p04 = -42.86;

        x = (x - 2.215) / 0.127;
        y = (y - 2.485) / 0.2771;
    }

    if (T <= 2.43 && rho <= 2.01 && rho >= (41.827 * pow(T, 3) - 341.32 * pow(T, 2) + 937.89 * T - 866.24)) { // reg 0b
        p00 = 1.984e+04;
        p10 = -4535;
        p01 = 255;
        p20 = 667.7;
        p11 = -84.19;
        p02 = 165;
        p30 = -51.74;
        p21 = -36.02;
        p12 = 84.42;
        p03 = -276.1;
        p40 = 5.095;
        p31 = 4.106;
        p22 = -20.62;
        p13 = 74.88;
        p04 = -160.2;
        p50 = -1.384;
        p41 = 0.04435;
        p32 = -1.939;
        p23 = 22.81;
        p14 = -14.36;
        
        x = (x - 2.181) / 0.1175;
        y = (y + 3.14) / 3.265;
    }

    if (rho > 2.15 && T > 2.43 && T <= 2.8 && rho < (-12.45 * pow(T, 4) + 118.3 * pow(T, 3) - 419.24 * pow(T, 2) + 655.96 * T - 378.81)) {
        p00 = 8491;
        p10 = -401.7;
        p01 = -1060;
        p20 = 75.72;
        p11 = 323.8;
        p02 = -278.7;
        p21 = -41.37;
        p12 = 88.17;
        p03 = -56;
        p22 = -14.5;
        p13 = 13.14;
        p04 = -7.816;

        x = (x - 2.612) / 0.1047;
        y = (y - 2.549) / 0.2306;
    }

    if (rho >= (41.827 * pow(T, 3) - 341.32 * pow(T, 2) + 937.89 * T - 866.24) && T > 2.43 && T <= 2.8 && rho <= 2.15) { // 1b
        p00 = 1.105e+04;
        p10 = -970;
        p01 = -59.9;
        p20 = 154.6;
        p11 = 0.3956;
        p02 = -150;
        p21 = -40.05;
        p12 = 71.06;
        p03 = -140.5;
        p22 = -10.68;
        p13 = 32.03;
        p04 = -37.83;

        x = (x - 2.553) / 0.08687;
        y = (y - 0.7782) / 1.018;
    }

    if (T > 2.15 && T <= 2.81 && rho < (41.827 * pow(T, 3) - 341.32 * pow(T, 2) + 937.89 * T - 866.24)) { // 3aI
        p00 = -2.217e+04;
        p10 = 3.134e+04;
        p01 = -1596;
        p20 = -1.087e+04;
        p11 = 308.3;
        p02 = -21.04;
        p30 = 1414;
        p21 = -45.77;
        p12 = 4.962;
        p03 = -0.4591;
    }

    if (T > 2.81 && T <= 3.18 && rho <= (3.7695 * pow(T, 4) - 52.814 * pow(T, 3) + 275.6 * pow(T, 2) - 636.26 * T + 550.74)) { // 3aII
        p00 = -2.217e+04;
        p10 = 3.134e+04;
        p01 = -1596;
        p20 = -1.087e+04;
        p11 = 308.3;
        p02 = -21.04;
        p30 = 1414;
        p21 = -45.77;
        p12 = 4.962;
        p03 = -0.4591;
    }
    if (T >= 3.12 && T < 3.26 && rho <= (3.7695 * pow(T, 4) - 52.814 * pow(T, 3) + 275.6 * pow(T, 2) - 636.26 * T + 550.74)) { // 3b
        p00 = 1.777e+04; // (1.775e+04, 1.778e+04)
        p10 = 16.08; // (-15.23, 47.38)
        p01 = -3406; // (-3437, -3376)
        p20 = 8.772; // (-15.58, 33.12)
        p11 = 590.4; // (571, 609.8)
        p02 = -864; // (-887.3, -840.7)
        p30 = 198.7; // (161.8, 235.5)
        p21 = -322.4; // (-349.9, -294.9)
        p12 = 98.94; // (72.24, 125.6)
        p03 = 3.642; // (-29.83, 37.11)
        p40 = -48.7; // (-57.55, -39.84)
        p31 = -115.3; // (-122.8, -107.9)
        p22 = 704.8; // (697.7, 712)
        p13 = -1165; // (-1172, -1157)
        p04 = 725.6; // (717.3, 733.8)
        p50 = -36.55; // (-47.41, -25.68)
        p41 = 91.85; // (82.99, 100.7)
        p32 = 24.48; // (16.17, 32.78)
        p23 = -393.4; // (-401.5, -385.2)
        p14 = 578.5; // (570.2, 586.7)
        p05 = -262.8; // (-272.3, -253.3)
        
        x = (x - 3.185) / 0.04031;
        y = (y + 4.257) / 3.319;
    }

    if (T < 3.7 && T >= 3.26 && rho > (-11.389 * pow(T, 2) + 96.06 * T - 198.63) && rho < (3.7695 * pow(T, 4) - 52.814 * pow(T, 3) + 275.6 * pow(T, 2) - 636.26 * T + 550.74)) { // 3c
        p00 = 1.595e+04; // (1.595e+04, 1.595e+04)
        p10 = 1224; // (1223, 1225)
        p01 = -2440; // (-2441, -2439)
        p20 = 682.2; // (680.3, 684.1)
        p11 = -1237; // (-1240, -1234)
        p02 = 531.5; // (529.6, 533.5)
        p30 = 430.1; // (429, 431.2)
        p21 = -1533; // (-1535, -1531)
        p12 = 1559; // (1557, 1561)
        p03 = -473.7; // (-474.9, -472.6)
        p40 = 96.76; // (95.68, 97.83)
        p31 = -756; // (-759, -753)
        p22 = 1419; // (1415, 1424)
        p13 = -983; // (-986.1, -979.9)
        p04 = 215.2; // (214.2, 216.3)
        p50 = -1.358; // (-1.645, -1.071)
        p41 = -64.54; // (-65.73, -63.34)
        p32 = 349.4; // (346.9, 351.9)
        p23 = -466; // (-468.5, -463.4)
        p14 = 201.2; // (200, 202.5)
        p05 = -24.29; // (-24.59, -23.98)

        x = (x - 3.396) / 0.1021;
        y = (y + 1.304) / 1.809;
    }    
    if (T < 3.7 && T >= 3.26 && rho < (-11.389 * pow(T, 2) + 96.06 * T - 198.63) && rho > (-0.4373 * pow(T, 2) + 26.273 * T - 90.521)) { // 3d
        p00 = 3.374e+04;
        p10 = 2.022e+04;
        p01 = -2.578e+04;
        p20 = -1.239e+04;
        p11 = 1.469e+04;
        p02 = -4722;
        p30 = -2.869e+04;
        p21 = 8.852e+04;
        p12 = -8.845e+04;
        p03 = 2.927e+04;
        p40 = 1.619e+04;
        p31 = -4.407e+04;
        p22 = 4.27e+04;
        p13 = -1.827e+04;
        p04 = 3371;
        p50 = 3.603e+04;
        p41 = -1.64e+05;
        p32 = 2.932e+05;    
        p23 = -2.589e+05;
        p14 = 1.137e+05;
        p05 = -1.998e+04;
        
        x = (x - 3.423) / 0.1109;
        y = (y + 4.552) / 2.37;
    }

    if (T < 3.7 && T >= 3.26 && rho <= (-0.4373 * pow(T, 2) + 26.273 * T - 90.521)) { // 3e
        p00 = 4.631e+04;
        p10 = 504.7;
        p01 = -7875;
        p20 = 118.2;
        p11 = -361;
        p02 = 362.8;
        p30 = -39.44;
        p21 = -87.72;
        p12 = 139.6;
        p03 = -92.31;
        p40 = -99.7;
        p31 = 510.4;
        p22 = -963;
        p13 = 876.1;
        p04 = -378.5;
        p50 = 79.56;
        p41 = -252.2;
        p32 = 445.7;
        p23 = -560.7;
        p14 = 432.9;
        p05 = -118.8;
        
        x = (x - 3.543) / 0.1071;
        y = (y + 6.465) / 2.497;
    }

    if (T >= 3.7 && rho <= (3.7695 * pow(T, 4) - 52.814 * pow(T, 3) + 275.6 * pow(T, 2) - 636.26 * T + 550.74) && rho > (-7.9704 * pow(T, 2) + 79.442 * T - 194.84) && rho >= (32.543 * T - 138.57)) { // 3f
        p00 = 3.754e+04;
        p10 = 1033;
        p01 = -7900;
        p20 = 977.7;
        p11 = -1722;
        p02 = 804;
        p30 = 788.4;
        p21 = -2795;
        p12 = 2949;
        p03 = -943;
        p40 = 149.4;
        p31 = -1434;
        p22 = 2791;
        p13 = -1942;
        p04 = 389.4;
        p50 = 79.56;
        p41 = -142;
        p32 = 856.9;
        p23 = -1251;
        p14 = 591.8;
        p05 = -90.07;
        
        x = (x - 3.884) / 0.1339;
        y = (y + 3.14) / 2.395;
    }    
    if (T >= 3.7 && rho <= (-7.9704 * pow(T, 2) + 79.442 * T - 194.84) && rho >= (72.232 * pow(T, 3) - 859.81 * pow(T, 2) + 3429.8 * T - 4591.5)) { // 3g
        p00 = 6.692e+04;
        p10 = 2.596e+04;
        p01 = -3.301e+04;
        p20 = -5133;
        p11 = -1530;
        p02 = 3023;
        p30 = -2.375e+04;
        p21 = 6.254e+04;
        p12 = -5.359e+04;
        p03 = 1.549e+04;
        p40 = -3448;
        p31 = 2.604e+04;
        p22 = -4.728e+04;
        p13 = 3.228e+04;
        p04 = -7464;
        p41 = 2.217e+04;
        p32 = -7.909e+04;
        p23 = 1.031e+05;
        p14 = -5.827e+04;
        p05 = 1.207e+04;

        x = (x - 3.938) / 0.117;
        y = (y + 6.581) / 2.086;
    }

    if (T >= 3.7 && T <= 4.26 && ((rho < (72.232 * pow(T, 3) - 859.81 * pow(T, 2) + 3429.8 * T - 4591.5) && T < 4.21) ||
        (T >= 4.21 && rho < (32.543 * T - 138.57)))) { // 3h
        p00 = 8.127e+04;
        p10 = 971.3;
        p01 = -1.381e+04;
        p20 = -17.03;
        p11 = 242.7;
        p02 = 40.66;
        p30 = 201.8;
        p21 = -1284;
        p12 = 2445;
        p03 = -1675;
        p40 = -333.7;
        p31 = 1485;
        p22 = -2380;
        p13 = 1877;
        p04 = -830.1;
        p50 = 155.5;
        p41 = -518;
        p32 = 873.6;
        p23 = -809.1;
        p14 = 158.2;
        p05 = 193.9;

        x = (x - 4.112) / 0.1046;
        y = (y + 6.94) / 2.241;
    }

    if (T > 4.26 && rho <= (-1.7922 * pow(T, 3) + 25.745 * pow(T, 2) - 124.2 * T + 200.36)) {
        p00 = 7.63e+04;
        p10 = 2107;
        p01 = -1.706e+04;
        p20 = 134.9;
        p11 = -431.2;
        p02 = -80.25;
        p30 = -8.168;
        p21 = -6.508;
        p12 = -106.4;
        p03 = -281.2;
        p40 = -45.09;
        p31 = 103.6;
        p22 = -209.5;
        p13 = 462;
        p04 = 139;
        p50 = 19.11;
        p41 = -26.75;
        p32 = 44.75;
        p23 = -94.98;
        p14 = 243.6;
        p05 = 146.2;

        x = (x - 4.624) / 0.2168;
        y = (y + 5.325) / 2.707;
    }
    if (T >= 2.41 && T <= 2.8 && rho >= (-12.45 * pow(T, 4) + 118.3 * pow(T, 3) - 419.24 * pow(T, 2) + 655.96 * T - 378.81) &&
        rho < (5.2542 * pow(T, 3) - 41.011 * pow(T, 2) + 106.95 * T - 90.052)) { // phase 4
        p00 = -1.635e+06;
        p10 = 8.573e+05;
        p01 = 8.15e+05;
        p20 = -9.464e+04;
        p11 = -3.836e+05;
        p02 = -8.947e+04;
        p30 = -947.1;
        p21 = 3.258e+04;
        p12 = 3.383e+04;
        p03 = -1227;
    }

    if (rho >= 3.91 && T <= 2.61) { // 5aI
        p00 = -2.725e+04;
        p10 = 3.072e-10;
        p01 = 3.391e+04;
        p20 = -2.681e-11;
        p11 = -4.453e-11;
        p02 = -6613;
    }

    if (rho >= 3.91 && T > 2.61) { // 5aII
        p00 = -2.551e+06;
        p10 = 2.973e+05;
        p01 = 2.077e+06;
        p20 = 1.303e+05;
        p11 = -4.788e+05;
        p02 = -4.642e+05;
        p30 = -2.39e+04;
        p21 = 1.281e+04;
        p12 = 9.964e+04;
        p03 = 3.239e+04;
        p40 = 595.5;
        p31 = 2350;
        p22 = -4130;
        p13 = -5363;
    }

    if (rho < 3.91 && rho >= 3.39 && T >= 3.35) { // 5b
        p00 = 1.389e+04;
        p10 = 4440;
        p01 = -1155;
        p20 = -566;
        p11 = 376.7;
        p02 = -100.9;
        p30 = -126;
        p21 = 45.79;
        p12 = 72.35;

        x = (x - 4.175) / 0.4792;
        y = (y - 3.65) / 0.153;
    }

    if (((T >= 2.81 && T < 3.12) && (rho <= (8.9265 * pow(T, 4) - 110.55 * pow(T, 3) + 511.88 * pow(T, 2) - 1050 * T + 808.23)) &&
        rho > (3.7695 * pow(T, 4) - 52.814 * pow(T, 3) + 275.6 * pow(T, 2) - 636.26 * T + 550.74)) ||
        (T >= 3.12 && T < 3.35 && rho <= 3.38 && rho > (3.7695 * pow(T, 4) - 52.814 * pow(T, 3) + 275.6 * pow(T, 2) - 636.26 * T + 550.74))) { // 5c
        p00 = 9790;
        p10 = 794.1;
        p01 = -1054;
        p20 = -32.32;
        p11 = 31.67;
        p02 = -256.4;
        p30 = 29;
        p21 = 57.09;
        p12 = 92.99;
        p03 = -8.592;
        p40 = 6.919;
        p31 = 4.675;
        p22 = 20.4;
        p13 = 55.59;
        p04 = -43.54;
        p50 = -1.047;
        p41 = -6.047;
        p32 = -18.59;
        p23 = -33.15;
        p14 = -11.51;
        p05 = -30.98;

        x = (x - 3.099) / 0.1511;
        y = (y - 2.476) / 0.5138;
    }

    if (T >= 3.35 && rho < 3.39 && rho >= 2.65) { // 5d
        p00 = 1.709e+04;
        p10 = 3763;
        p01 = -593.4;
        p20 = -1149;
        p11 = -49.49;
        p02 = -80.06;
        p30 = -177.1;
        p21 = 92.47;
        p12 = -29.55;
        p03 = 26.64;
        p40 = 230.4;
        p31 = 33.71;
        p22 = 1.571;
        p13 = 27.38;
        p04 = -14.12;
        p41 = -20.57;
        p32 = 4.922;
        p23 = -15.8;
        p14 = 20.76;
        p05 = -11.17;

        x = (x - 4.175) / 0.4792;
        y = (y - 3.02) / 0.2107;
    }    
    if (rho < 2.65 && T >= 3.35 &&
        ((T <= 4.26 && rho > (3.7695 * pow(T, 4) - 52.814 * pow(T, 3) + 275.6 * pow(T, 2) - 636.26 * T + 550.74)) ||
         (T > 4.26 && rho > (-1.7922 * pow(T, 3) + 25.745 * pow(T, 2) - 124.2 * T + 200.36)))) {
        // 5e
        if (rho > 1.33) { // 5eI
            p00 = 2.136e+04;
            p10 = 5226;
            p01 = -2423;
            p20 = 102.3;
            p11 = -1599;
            p02 = -23.1;
            p30 = 1416;
            p21 = -1044;
            p12 = 153.5;
            p03 = 313.7;
            p40 = 679.1;
            p31 = -714;
            p22 = -200.9;
            p13 = 251;
            p04 = -1.438;
            p41 = -123.7;
            p32 = -183.2;
            p23 = 115.6;
            p14 = -17.8;
            p05 = -63.94;

            x = (x - 4.175) / 0.4792;
            y = (y - 1.99) / 0.3839;
        } else { // 5eII
            p00 = 3.406e+04;
            p10 = 1.354e+04;
            p01 = -3943;
            p20 = 7918;
            p11 = -1281;
            p02 = -519.7;
            p30 = 495.7;
            p21 = 1270;
            p12 = -1150;
            p03 = -28.08;
            p40 = -2267;
            p31 = 1065;
            p22 = 488.9;
            p13 = 403;
            p04 = 235.7;
            p50 = -769.1;
            p41 = -104.1;
            p32 = 59.78;
            p23 = -378.4;
            p14 = -195.9;
            p05 = -7.794;

            x = (x - 4.443) / 0.383;
            y = (y - 0.4717) / 0.5813;
        }
    }

    if (T >= 2.0 && T <= 2.63 && rho >= 3.19 && rho <= 3.47) { // M7A
        p00 = 698.7;
        p10 = 484.4;
        p01 = -279.3;
        p20 = 88.82;
        p11 = -155.2;
        p02 = 45.49;
        p30 = -4.357;
        p21 = -11.68;
        p12 = 23.48;

        x = (x - 2.31) / 0.1824;
        y = (y - 3.331) / 0.08323;
    }

    if (T < 3.35 && T > 3.12 && rho <= 3.47 && rho > 3.38) { // M7BI
        p00 = 2630;
        p10 = 1068;
        p01 = -520.3;
        p20 = 165.1;
        p11 = -78.71;
        p02 = 71.5;

        x = (x - 2.854) / 0.1732;
        y = (y - 3.377) / 0.07291;
    }
    if (T <= 3.12 && T >= 2.81 && rho <= 3.47 && rho > (8.9265 * pow(T, 4) - 110.55 * pow(T, 3) + 511.88 * pow(T, 2) - 1050 * T + 808.23)) {
        // M7BII
        p00 = 2630;
        p10 = 1068;
        p01 = -520.3;
        p20 = 165.1;
        p11 = -78.71;
        p02 = 71.5;

        x = (x - 2.854) / 0.1732;
        y = (y - 3.377) / 0.07291;
    }

    if (T < 2.81 && T > 2.71 && rho <= 3.47 && rho > 3.18 &&
        (rho >= (5.2542 * pow(T, 3) - 41.011 * pow(T, 2) + 106.95 * T - 90.052))) {
        // M7BIII (between phase 4 and ice 1)
        p00 = 2630;
        p10 = 1068;
        p01 = -520.3;
        p20 = 165.1;
        p11 = -78.71;
        p02 = 71.5;

        x = (x - 2.854) / 0.1732;
        y = (y - 3.377) / 0.07291;
    }

    if (T > 2.63 && T <= 2.71 && rho <= 3.47 && rho >= 3.19) {
        // M7BIV (between ice 7 and ice 1)
        p00 = 2630;
        p10 = 1068;
        p01 = -520.3;
        p20 = 165.1;
        p11 = -78.71;
        p02 = 71.5;

        x = (x - 2.854) / 0.1732;
        y = (y - 3.377) / 0.07291;
    }

    if (rho <= 3.64 && rho > 3.47 && T < 3.35) {
        // Ice 1
        p00 = 949.3;
        p10 = 1321;
        p01 = -152.5;
        p20 = 631.2;
        p11 = -90.43;
        p02 = 12.78;
        p30 = 93.7;
        p21 = -4.125;
        p12 = 4.974;

        x = (x - 2.67) / 0.3898;
        y = (y - 3.555) / 0.05189;
    }

    if (rho > 3.64 && rho <= 3.75 && T < 2.47) {
        // Ice 2
        p00 = -1.13e+06;
        p10 = 1.199e+04;
        p01 = 5.305e+05;
        p11 = -3190;
        p02 = -6.056e+04;
    }

    if (rho < 3.91 && rho > 3.75 && T < 2.47) {
        // Ice 3
        p00 = -2.196e+05;
        p10 = 2.592;
        p01 = 1.335e+05;
        p11 = -0.672;
        p02 = -1.938e+04;
    }

    if (rho < 3.91 && rho > 3.75 && T >= 2.47 && T < 3.35) {
        // Ice 4
        p00 = 1.523e+06;
        p10 = -9.901e+05;
        p01 = -9.609e+05;
        p20 = -1.457e+05;
        p11 = 9.717e+05;
        p02 = 7.121e+04;
        p30 = 1.451e+05;
        p21 = -2.462e+05;
        p12 = -7.504e+04;
        p40 = -1.676e+04;
        p31 = 1.406e+04;
        p22 = 1.52e+04;
    }

    if (rho > 3.64 && rho <= 3.75 && T >= 2.47 && T < 3.35) {
        // Ice 5
        p00 = 2536;
        p10 = 39.14;
        p01 = 428.9;
        p20 = 705.5;
        p11 = -770.5;
        p02 = -42.18;
        p30 = 20.41;
        p21 = 225.1;
        p12 = -0.09814;

        x = (x - 2.905) / 0.2541;
        y = (y - 3.7) / 0.03164;
    }
    if (rho > 2.96 && rho < 3.19 && T <= 2.43) {
        // Ice 6b
        p00 = 1280;
        p10 = 357.5;
        p01 = 95.98;
        p20 = 2.97;
        p11 = -18.69;
        p02 = 165.9;
        p30 = -11.84;
        p21 = -44.07;
        p12 = 44.24;
        p03 = -97.18;
        p40 = 66.96;
        p31 = -52.57;
        p22 = -44.32;
        p13 = 53.23;
        p04 = -64.5;
        p50 = 43.29;
        p41 = -42.72;
        p32 = -48.66;
        p23 = 61.54;
        p14 = 6.362;

        x = (x - 2.215) / 0.127;
        y = (y - 3.075) / 0.06348;
    }

    if (rho > (5.2542 * pow(T, 3) - 41.011 * pow(T, 2) + 106.95 * T - 90.052) && T > 2.43 && T <= 2.8 && rho < 3.19) {
        // Ice 7
        p00 = 2993;
        p10 = 1075;
        p01 = -540.6;
        p20 = 146.1;
        p11 = -369.1;
        p02 = 221;
        p30 = -146.5;
        p21 = 86.27;
        p12 = -67.21;
        p40 = 13.6;
        p31 = 130.3;
        p22 = -137.8;

        x = (x - 2.544) / 0.07807;
        y = (y - 3.149) / 0.02604;
    }

    // Calculate s using the appropriate polynomials and variables
    s = p00 + p10 * x + p01 * y + p20 * pow(x, 2) + p11 * x * y + p02 * pow(y, 2) + p30 * pow(x, 3) + p21 * pow(x, 2) * y
        + p12 * x * pow(y, 2) + p03 * pow(y, 3) + p40 * pow(x, 4) + p31 * pow(x, 3) * y + p22 * pow(x, 2) * pow(y, 2)
        + p13 * x * pow(y, 3) + p04 * pow(y, 4) + p50 * pow(x, 5) + p41 * pow(x, 4) * y + p32 * pow(x, 3) * pow(y, 2)
        + p23 * pow(x, 2) * pow(y, 3) + p14 * x * pow(y, 4) + p05 * pow(y, 5);
    // s is in J/(kg*K)

    if (s <= 0) {
        s = 0;
    }
    s=10000.0 *s; // to erg/(g*K)
    return s;    

}