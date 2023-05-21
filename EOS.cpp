#include "EOS.h"
#define DEBUG_LEVEL 2
EOS::EOS():phasetype(""),eqntype(0), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), Debye_approx(false), thermal_type(0), rhotable(NULL), Ptable(NULL), bn(0), acc(NULL), spline(NULL), nline(0)
{
  density_extern=NULL;
  entropy_extern=NULL;
  dTdP = NULL;
}

EOS::EOS(string phaseinput, double params[][2], int length):phasetype(phaseinput),eqntype(0), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), Debye_approx(false), thermal_type(0), rhotable(NULL), Ptable(NULL), bn(0), acc(NULL), spline(NULL), nline(0)
{
  density_extern=NULL;
  entropy_extern=NULL;
  dTdP = NULL;
  
  for(int i=0;i<length;i++)
  {
    switch(lround(params[i][0]))
    {
    case 0:
      eqntype=round(params[i][1]);
      break;
    case 1:
      V0=params[i][1];
      break;
    case 2:
      K0=params[i][1];
      break;
    case 3:
      K0p=params[i][1];
      break;
    case 4:
      K0pp=params[i][1];
      break;
    case 5:
      mmol=params[i][1];
      break;
    case 6:
      P0=params[i][1];
      break;
    case 7:
      Theta0=params[i][1];
      break;
    case 8:
      gamma0=params[i][1];
      break;
    case 9:
      beta=params[i][1];
      break;
    case 10:
      gammainf=params[i][1];
      break;
    case 11:
      gamma0p=params[i][1];
      break;
    case 12:
      e0=params[i][1];
      break;
    case 13:
      g=params[i][1];
      break;
    case 14:
      n=round(params[i][1]);
      break;
    case 15:
      Z=round(params[i][1]);
      break;
    case 16:
      T0=params[i][1];
      break;
    case 17:
      alpha0=params[i][1];
      break;
    case 18:
      alpha1=params[i][1];
      break;
    case 19:
      xi=params[i][1];
      break;
    case 20:
      cp_a=params[i][1];
      break;
    case 21:
      cp_b=params[i][1];
      break;
    case 22:
      cp_c=params[i][1];
      break;
    case 23:
      Debye_approx = params[i][1]>0 ? true : false;
      break;
    case 24:
      thermal_type = params[i][1];
      break;
    case 25:
      at1 = params[i][1];
      break;
    case 26:
      at2 = params[i][1];
      break;
    case 27:
      at3 = params[i][1];
      break;
    case 28:
      at4 = params[i][1];
      break;
    case 29:
      ap1 = params[i][1];
      break;
    case 30:
      ap2 = params[i][1];
      break;
    case 31:
      ap3 = params[i][1];
      break;
    case 32:
      ap4 = params[i][1];
      break;
    default:
      cout<<"Error: Incorrect index "<<round(params[i][0])<<" for EOS constructor "<<phasetype<<endl;
      exit(1);
    };
  }
  if (eqntype == 6)
  {
    if (!gsl_finite(mmol))	// default mean molecular weight of gas.  Mix of hydrogen and helium
      mmol = 2.3;
    if (n<0)		// number of atom
      n = 2;
  }

  if (eqntype == 6)		// ideal gas
    thermal_type = 3;
  else if (thermal_type!=1 && thermal_type!=2 && thermal_type!=4 && gsl_finite(gamma0)) // thermal type not specified and have enough information to calculate using Grueneisen parameter
  {
    if (gsl_finite(Theta0))
    {
      if (gsl_finite(e0) && gsl_finite(g)) // Also have enough information to get Pe
	thermal_type = 7;
      else
	thermal_type = 6;
    }
    else
      thermal_type = 5;
  }
  else if (thermal_type!=1 && thermal_type!=2 && thermal_type!=4 && gsl_finite(cp(300)) && gsl_finite(alpha(10,300)) && gsl_finite(T0)) // thermal type not specified and have enough information to calculate using thermal expansion
    thermal_type = 9;
  
  if (eqntype >= 8)		// RTpress EOS style
  {
    cout<<"Error: Incorrect equation type. Type "<<eqntype<<" is an index for RTpress style EOS. It cannot be initialized by this constructor."<<endl;
    exit(1);
  }
}

EOS::EOS(string phaseinput, string filename):phasetype(phaseinput),eqntype(7), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), Debye_approx(false), thermal_type(0), bn(0)
{
  ifstream fin;
  string sline;
  fin.open(filename.c_str());
  density_extern=NULL;
  entropy_extern=NULL;
  dTdP = NULL;
  
  if(!fin)
  {
    if (verbose)
      cout<<"Warning: Failed to open input EOS file "<<filename<<" \n This MAY cause segmentation fault in the run time if this phase "<<phaseinput<<" is used in the runtime!!"<<endl;
    return;
  }

  nline = 0;
  getline(fin,sline);
  streampos beginpos=fin.tellg();

  while(getline(fin,sline))
  {
   if(!sline.empty())
     nline++;
  }
  fin.clear();
  fin.seekg(beginpos);

  rhotable=new double[nline];
  Ptable=new double[nline];

  // grl interpolation require strictly increasing x.  Reverse the array order if x is decreasing.
  
  fin>>rhotable[0]>>Ptable[0];
  fin>>rhotable[1]>>Ptable[1];
  fin.seekg(beginpos);

  if(Ptable[0]>Ptable[1])    
    for(int i=nline-1; i>=0; i--)
      fin>>rhotable[i]>>Ptable[i];
  else
    for(int i=0; i<nline; i++)
      fin>>rhotable[i]>>Ptable[i];
  
  fin.close();
  
  acc = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_steffen, nline);
  gsl_spline_init (spline, Ptable, rhotable, nline);
}

EOS::EOS(string phaseinput, double (*f)(double P, double T), double (*g)(double rho, double T)):phasetype(phaseinput),eqntype(0), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), Debye_approx(false), rhotable(NULL), Ptable(NULL),  bn(0), acc(NULL), spline(NULL), nline(0)
{
  density_extern=f;
  entropy_extern=g;
  dTdP = NULL;

  if (entropy_extern)
    thermal_type = 1;
  else
    thermal_type = 0;
}

EOS::EOS(string phaseinput, double *Plist, double *rholist, int len_list):phasetype(phaseinput),eqntype(7), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), Debye_approx(false), thermal_type(0), bn(0), nline(len_list)
{
  density_extern=NULL;
  entropy_extern=NULL;
  dTdP = NULL;

  rhotable=new double[nline];
  Ptable=new double[nline];

  // grl interpolation require strictly increasing x.  Reverse the array order if x is decreasing.
  
  if(Plist[0]>Plist[1])    
    for(int i=0; i<nline; i++)
    {
      rhotable[i] = rholist[nline-1-i];
      Ptable[i]   = Plist[nline-1-i];
    }
  else
    for(int i=0; i<nline; i++)
    {
      rhotable[i] = rholist[i];
      Ptable[i]   = Plist[i];
    }
  
  acc = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_steffen, nline);
  gsl_spline_init (spline, Ptable, rhotable, nline);
}

EOS::EOS(string phaseinput, double params[][2], double bparams[], int length, int blength):phasetype(phaseinput),eqntype(8), V0(numeric_limits<double>::quiet_NaN()), K0(numeric_limits<double>::quiet_NaN()), K0p(numeric_limits<double>::quiet_NaN()), K0pp(numeric_limits<double>::quiet_NaN()), mmol(numeric_limits<double>::quiet_NaN()), P0(0), Theta0(numeric_limits<double>::quiet_NaN()), gamma0(numeric_limits<double>::quiet_NaN()), beta(numeric_limits<double>::quiet_NaN()), gammainf(numeric_limits<double>::quiet_NaN()), gamma0p(numeric_limits<double>::quiet_NaN()), e0(numeric_limits<double>::quiet_NaN()), g(numeric_limits<double>::quiet_NaN()), T0(300), alpha0(numeric_limits<double>::quiet_NaN()), alpha1(0), xi(0), cp_a(numeric_limits<double>::quiet_NaN()), cp_b(0), cp_c(0), at1(numeric_limits<double>::quiet_NaN()), at2(numeric_limits<double>::quiet_NaN()), at3(numeric_limits<double>::quiet_NaN()), at4(numeric_limits<double>::quiet_NaN()), ap1(numeric_limits<double>::quiet_NaN()), ap2(numeric_limits<double>::quiet_NaN()), ap3(numeric_limits<double>::quiet_NaN()), ap4(numeric_limits<double>::quiet_NaN()), n(-1), Z(-1), Debye_approx(false), thermal_type(8), rhotable(NULL), Ptable(NULL), bn(blength), acc(NULL), spline(NULL), nline(0)		
{
  // construction EOS for RTpress
  density_extern=NULL;
  entropy_extern=NULL;
  dTdP = NULL;
  
  for(int i=0;i<length;i++)
  {
    switch(lround(params[i][0]))
    {
    case 0:
      eqntype=round(params[i][1]);
      break;
// eqntype has to be equal or larger than 8
    case 1:
      V0=params[i][1];
      break;
    case 2:
      K0=params[i][1];
      break;
    case 3:
      K0p=params[i][1];
      break;
    case 4:
      K0pp=params[i][1];
      break;
    case 5:
      mmol=params[i][1];
      break;
    case 6:
      P0=params[i][1];
      break;
    case 7:
      Theta0=params[i][1];
      break;
    case 8:
      gamma0=params[i][1];
      break;
    case 9:
      beta=params[i][1];
      break;
    case 10:
      gammainf=params[i][1];
      break;
    case 11:
      gamma0p=params[i][1];
      break;
    case 12:
      e0=params[i][1];
      break;
    case 13:
      g=params[i][1];
      break;
    case 14:
      n=round(params[i][1]);
      break;
    case 15:
      Z=round(params[i][1]);
      break;
    case 16:
      T0=params[i][1];
      break;
    case 17:
      alpha0=params[i][1];
      break;
    case 18:
      alpha1=params[i][1];
      break;
    case 19:
      xi=params[i][1];
      break;
    case 20:
      cp_a=params[i][1];
      break;
    case 21:
      cp_b=params[i][1];
      break;
    case 22:
      cp_c=params[i][1];
      break;
    case 23:
      Debye_approx = params[i][1]>0 ? true : false;
      break;
    case 24:
      break;		// thermal_type has to be 8
    case 25:
      at1 = params[i][1];
      break;
    case 26:
      at2 = params[i][1];
      break;
    case 27:
      at3 = params[i][1];
      break;
    case 28:
      at4 = params[i][1];
      break;
    case 29:
      ap1 = params[i][1];
      break;
    case 30:
      ap2 = params[i][1];
      break;
    case 31:
      ap3 = params[i][1];
      break;
    case 32:
      ap4 = params[i][1];
      break;
    default:
      cout<<"Error: Incorrect index "<<round(params[i][0])<<" for EOS constructor "<<phasetype<<endl;
      exit(1);
    };
  }

  if (eqntype<8)
  {
    cout<<"Error: Incorrect equation type. Type "<<eqntype<<" cannot be initialized by the constructor for RTpress style EOS."<<endl;
    exit(1);
  }

  // fitted polynomial parameters of the thermal coefficients b(V) in erg/mol.  Convert eV/atom to erg/mol need to multiply eV_erg*n*NA. For example, for MgSiO3, 0.9821 eV/atom = 4.824E12 *0.9821 erg/mol = 4.738E12 erg/mol.
  b = new double[bn];
  for (int i=0; i<bn; i++)
    b[i] = bparams[i];
}


EOS::~EOS()
{
  if(rhotable)
    delete[] rhotable;
  if(Ptable)
    delete[] Ptable;
  if(spline)
    gsl_spline_free (spline);
  if(acc)
    gsl_interp_accel_free (acc);
  if(bn>0)
    delete[] b;
}

void EOS::modifyEOS(double params[][2], int length)  // modify the constructed EOS parameters
{
  if (density_extern)		// Warning for EOS modification can't overwrite external EOS functions.
  {
    if (verbose)
      cout<<"Warning: External EOS function is set for state "<<phasetype<<". The modification on the EOS fitting parameters won't change the EOS function.  Consider change the external function if really need to modify this EOS."<<endl;
    return;
  }

  int index;
  
  for(int i=0;i<length;i++)
  {
    index = lround(params[i][0]);
    
    switch(index)
    {
    case 0:
      eqntype=round(params[i][1]);
      break;
    case 1:
      V0=params[i][1];
      break;
    case 2:
      K0=params[i][1];
      break;
    case 3:
      K0p=params[i][1];
      break;
    case 4:
      K0pp=params[i][1];
      break;
    case 5:
      mmol=params[i][1];
      break;
    case 6:
      P0=params[i][1];
      break;
    case 7:
      Theta0=params[i][1];
      break;
    case 8:
      gamma0=params[i][1];
      break;
    case 9:
      beta=params[i][1];
      break;
    case 10:
      gammainf=params[i][1];
      break;
    case 11:
      gamma0p=params[i][1];
      break;
    case 12:
      e0=params[i][1];
      break;
    case 13:
      g=params[i][1];
      break;
    case 14:
      if(eqntype>=8)		// RTpress style
      {
	cout<<"Error: The number of atoms of a RTpress style EOS is not allowed to be modified."<<endl;
	exit(1);
      }
      n=round(params[i][1]);
      break;
    case 15:
      Z=round(params[i][1]);
      break;
    case 16:
      T0=params[i][1];
      break;
    case 17:
      alpha0=params[i][1];
      break;
    case 18:
      alpha1=params[i][1];
      break;
    case 19:
      xi=params[i][1];
      break;
    case 20:
      cp_a=params[i][1];
      break;
    case 21:
      cp_b=params[i][1];
      break;
    case 22:
      cp_c=params[i][1];
      break;
    case 23:
      Debye_approx = params[i][1]>0 ? true : false;
      break;
    case 24:
      thermal_type = params[i][1];
      break;
    case 25:
      at1 = params[i][1];
      break;
    case 26:
      at2 = params[i][1];
      break;
    case 27:
      at3 = params[i][1];
      break;
    case 28:
      at4 = params[i][1];
      break;
    case 29:
      ap1 = params[i][1];
      break;
    case 30:
      ap2 = params[i][1];
      break;
    case 31:
      ap3 = params[i][1];
      break;
    case 32:
      ap4 = params[i][1];
      break;

    default:
      cout<<"Error: Incorrect index "<<round(params[i][0])<<" for EOS constructor "<<phasetype<<endl;
      exit(1);
    };
  }

  if (eqntype == 6)		// ideal gas
    thermal_type = 3;
  else if (thermal_type!=1 && thermal_type!=2 && thermal_type!=4 && gsl_finite(gamma0)) // thermal type not specified and have enough information to calculate using Grueneisen parameter
  {
    if (gsl_finite(Theta0))
    {
      if (gsl_finite(e0) && gsl_finite(g)) // Also have enough information to get Pe
	thermal_type = 7;
      else
	thermal_type = 6;
    }
    else
      thermal_type = 5;
  }
  else if (thermal_type!=1 && thermal_type!=2 && thermal_type!=4 && gsl_finite(cp(300)) && gsl_finite(alpha(10,300)) && gsl_finite(T0)) // thermal type not specified and have enough information to calculate using thermal expansion
    thermal_type = 9;

  
  if (eqntype >= 8)		// RTpress EOS style
    thermal_type = 8;
}

void EOS::modifyEOS(int index, double value)	     // modify one value of the EOS
{
  if (density_extern)		// Warning for EOS modification can't overwrite external EOS functions.
  {
    if (verbose)
      cout<<"Warning: External EOS function is set for state "<<phasetype<<". The modification on the EOS fitting parameters won't change the EOS function.  Consider change the external function if really need to modify this EOS."<<endl;
    return;
  }

  switch(index)
  {
  case 0:
    eqntype=round(value);
    break;
  case 1:
    V0=value;
    break;
  case 2:
    K0=value;
    break;
  case 3:
    K0p=value;
    break;
  case 4:
    K0pp=value;
    break;
  case 5:
    mmol=value;
    break;
  case 6:
    P0=value;
    break;
  case 7:
    Theta0=value;
    break;
  case 8:
    gamma0=value;
    break;
  case 9:
    beta=value;
    break;
  case 10:
    gammainf=value;
    break;
  case 11:
    gamma0p=value;
    break;
  case 12:
    e0=value;
    break;
  case 13:
    g=value;
    break;
  case 14:
    if(eqntype>=8)		// RTpress style
    {
      cout<<"Error: The number of atoms of a RTpress style EOS is not allowed to be modified."<<endl;
      exit(1);
    }
    n=round(value);
    break;
  case 15:
    Z=round(value);
    break;
  case 16:
    T0=value;
    break;
  case 17:
    alpha0=value;
    break;
  case 18:
    alpha1=value;
    break;
  case 19:
    xi=value;
    break;
  case 20:
    cp_a=value;
    break;
  case 21:
    cp_b=value;
    break;
  case 22:
    cp_c=value;
    break;
  case 23:
    Debye_approx = value>0 ? true : false;
    break;
  case 24:
    thermal_type = value;
    break;
  case 25:
    at1 = value;
    break;
  case 26:
    at2 = value;
    break;
  case 27:
    at3 = value;
    break;
  case 28:
    at4 = value;
    break;
  case 29:
    ap1 = value;
    break;
  case 30:
    ap2 = value;
    break;
  case 31:
    ap3 = value;
    break;
  case 32:
    ap4 = value;
    break;
  
  default:
    cout<<"Error: Incorrect index "<<index<<" for EOS constructor "<<phasetype<<endl;
    exit(1);
  };

  if (eqntype == 6)		// ideal gas
    thermal_type = 3;
  else if (thermal_type!=1 && thermal_type!=2 && thermal_type!=4 && gsl_finite(gamma0)) // thermal type not specified and have enough information to using Grueneisen parameter
  {
    if (gsl_finite(Theta0))
    {
      if (gsl_finite(e0) && gsl_finite(g)) // Also have enough information to get Pe
	thermal_type = 7;
      else
	thermal_type = 6;
    }
    else
      thermal_type = 5;
  }
  else if (thermal_type!=1 && thermal_type!=2 && thermal_type!=4 && gsl_finite(cp(300)) && gsl_finite(alpha(10,300)) && gsl_finite(T0)) // thermal type not specified and have enough information to calculate using thermal expansionmake
    thermal_type = 9;
  
  if (eqntype >= 8)		// RTpress EOS style
  {
    if (!gsl_finite(beta))
      beta = 0.6;
    thermal_type = 8;
  }
}

double EOS::BM3(double rho)
// input rho in g/cm^3, return pressure in GPa
{
  if (!gsl_finite(V0) || !gsl_finite(K0))
  {
    cout<<"Error: "<<this->phasetype<<" missing input parameter(s) for BM3 EOS."<<endl;
    exit(1);
  }
  if (!gsl_finite(K0p))
    K0p=4;			// reduced to BM2
  
  double V = mmol/rho;

  double eta = V0/V;
  
  return P0 + 1.5*K0 * (pow(eta,7./3.)-pow(eta,5./3.)) * (1+0.75*(K0p-4)*(pow(eta,2./3.)-1));
}


double EOS::BM4(double rho)
{
  if (!gsl_finite(V0) || !gsl_finite(K0) || !gsl_finite(K0p) || !gsl_finite(K0pp))
  {
    cout<<"Error: "<<this->phasetype<<" missing input parameter(s) for BM4 EOS."<<endl;
    exit(1);
  }
  double V = mmol/rho;

  double eta = V0/V;
  return P0 + 1.5*K0 * (pow(eta,7./3.)-pow(eta,5./3.)) * (1 + 0.75*(K0p-4)*(pow(eta,2./3.)-1) + 0.375*sq(pow(eta,2./3.)-1)*(K0*K0pp+K0p*(K0p-7)+143./9.));
}


double EOS::Vinet(double rho)
{
  if (!gsl_finite(V0) || !gsl_finite(K0) || !gsl_finite(K0p))
  {
    cout<<"Error: "<<this->phasetype<<" missing input parameter(s) for Vinet EOS."<<endl;
    exit(1);
  }

  double V = mmol/rho;

  double eta = V0/V;
  return P0 + 3*K0 * pow(eta,2./3.) * (1-pow(eta,-1./3.)) * exp(1.5*(K0p-1)*(1-pow(eta,-1./3.)));
}

double EOS::Holzapfel(double rho)
{
  if (!gsl_finite(V0) || !gsl_finite(K0) || !gsl_finite(K0p) || Z<0)
  {
    cout<<"Error: "<<this->phasetype<<" missing input parameter(s) for Holzapfel EOS."<<endl;
    exit(1);
  }
  
  double V = mmol/rho;

  double x = pow(V/V0,1./3.);
  double c0= -log(3*K0/(1003.6*pow(Z/V0,5./3.)));
  double c2= 1.5*(K0p-3)-c0;

  return 3*K0 * pow(x,-5) * (1-x) * exp(c0*(1-x)) * (1+c2*x*(1-x));
}

double EOS::Keane(double rho)
// Keane 1954, Aust. J. Phys. 7, 322-333
{
  if (!gsl_finite(V0) || !gsl_finite(K0) || !gsl_finite(K0p) || !gsl_finite(gammainf))
  {
    cout<<"Error: "<<this->phasetype<<" missing input parameter(s) for Keane EOS."<<endl;
    exit(1);
  }

  double V = mmol/rho;

  double y = V0/V;
  double Kinfp = 2*(gammainf+1./6.);

  return K0p*K0/sq(Kinfp)*(pow(y,Kinfp)-1) - (K0p-Kinfp)*K0/Kinfp*log(y);
}


void EOS::DebyeT(double x, double &gamma, double &Theta)  // return the Grueneisen parameter, Debye temperature or Einstein temperature according to Altshuler form.
// If Theta0 is not available, a Debye temperature scaling factor is returned
{
  if ((!gsl_finite(V0) || thermal_type < 4) && !(thermal_type==2 && gsl_finite(gamma0) && gsl_finite(Theta0))) // don't have thermal pressure data
  {
    cout<<"Error: Cannot calculate the Debye temperature for phase "<<phasetype<<".  Lack of necessary information."<<endl;
    gamma = numeric_limits<double>::quiet_NaN();
    Theta = numeric_limits<double>::quiet_NaN();
    return;
  }

  // set up default parameters for gammainf, beta, and n if not provided by the input data.
  
  if (!gsl_finite(gammainf))	// According to Al'tshuler et al. 1987, gammainf = 2/3 for all elements except alkali elments, for which gammainf = 0.5
    gammainf = 2./3.;
  if (!gsl_finite(beta))	// Altshuler form.
    beta = gamma0 / (gamma0-gammainf);

  gamma = gammainf + (gamma0-gammainf)*pow(x,beta);

  if (thermal_type == 5 || thermal_type == 4)
    Theta = pow(x,-gammainf) * exp((gamma0-gammainf)/beta*(1-pow(x,beta))); 
  else
    Theta = Theta0 * pow(x,-gammainf) * exp((gamma0-gammainf)/beta*(1-pow(x,beta))); // gammma = - d ln(Theta) / d ln(V)
}


double EOS::Pth(double V, double T)
// calculate the thermal pressure in GPa, ref. Bouchet et al. 2013 PRB 87, 094102, Shim & Duffy, 2000, American Mineralogist
// Or RTpress style Pth.
{
  if (!gsl_finite(V0) || thermal_type < 2 || thermal_type == 3)  // don't have thermal pressure data
    return 0;

  if ((thermal_type == 2 && !gsl_finite(gamma0)) || thermal_type == 9)
    return 0;

  if (gsl_finite(T) && getthermal() == 8) // RTpress style
  {
    double cf = 1E-10;
// conversion factor from erg/cm^3 to GPa
    double T_OS = TOS(V);
    double fTp_OS = fTp(T_OS);
    double bVpV = bVp(V);	// in  erg/cm^3 (microbar)

    return - cf*bVpV*fT(T) + cf*gamma0S(V)*(T-T0)/V*Cv(V,T_OS) + cf*bVpV/(beta-1)*(T*(fTp(T)-fTp_OS) - T0*(fTp(T0)-fTp_OS));
  }

  // set up default parameters for n if not provided by the input data.
  if (n<0)
    n = 1;

  double gamma, Theta, x = V/V0;
  DebyeT(x, gamma, Theta);

  if (thermal_type == 4 || thermal_type == 5 || (thermal_type == 2 && !gsl_finite(Theta0)))	// Only gamma available
    // ref. Ono & Oganov 2005, Eq. 3, Earth Planet. Sci. Lett. 236, 914
    return 3E-10*gamma*n*R*(T-T0)/V;
  
  double Eth = 3*n*R*T*gsl_sf_debye_3(Theta/T);
  double Eth0 = 3*n*R*T0*gsl_sf_debye_3(Theta/T0);

  if (thermal_type == 7  || (thermal_type==2 && gsl_finite(e0) && gsl_finite(g)))
    return 1E-10*gamma*(Eth-Eth0)/V + 1.5E-16*n*R*e0*pow(x,g)*g/V*(sq(T)-sq(T0));
  // 1E-16 = 1E-10 * 1E-6.  GPa -> microbar and e0 in 10^-6 K^-1
  else // don't have enough information to get Pe
    return 1E-10 * gamma*(Eth-Eth0)/V;	  // convert to GPa
}

double EOS::adiabatic_index()	    // get the adiabatic index for ideal gas.  Vibrational freedom is always ignored.
{
  if (eqntype != 6 || n<0)
  {
    cout<<"Error: "<<phasetype<<" is not ideal gas or number of atom per molecule unknown.  No adiabatic index applied."<<endl;
    return numeric_limits<double>::quiet_NaN();
  }

  if (n == 1)			// monatomic ideal gas
    return 5./3.;
  else if (n == 2)		//  diatomic gas and collinear molecules e.g. carbon dioxide
    return 1.4;
  else if (n == 0)		// isothermal atmosphere
    return 1.;
  else			// polyatomic gas
    return 4./3.;
}
  
double EOS::density(double P, double T, double rho_guess)
// input P in cgs (microbar), return density in g/cm^3
{
  if(!gsl_finite(P) || !gsl_finite(T)) // Check if P, T, or rho_guess is infinite or nan due to some error.  Stop code to avoid further error.
  {
    if (verbose)
      cout<<"Warning: Request density for "<<phasetype<<" at infinite/nan value.  P="<<P/1E10<<" T="<<T<<endl;
    return numeric_limits<double>::quiet_NaN();
  }

  int status;
  
  if(P < 0 || P > 1E16)		// unrealistic pressure
    return numeric_limits<double>::quiet_NaN();

  else if(density_extern)
    return density_extern(P, T);
  
  else if(eqntype == 7)		// interpolate an input file
  {

    P /= 1E10;
    double rho;
    
    status = gsl_spline_eval_e(spline, P, acc, &rho);

    if(status == GSL_EDOM)
    {
      if (verbose)
	cout<<"Warning: Pressure "<<P<<"GPa is outside the tabulated range for "<<this->phasetype<<". The density at the end point is returned"<<endl;
      if(P < Ptable[0])
	return rhotable[0];
      else
	return rhotable[nline-1];
    }
    else	
      return rho;
  }

  else if(eqntype == 6)		// ideal gas
  {
    if (!gsl_finite(mmol))
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;
    return P*mmol*mp/(kb*T);
  }

  else
  {
    if (!gsl_finite(mmol))
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;

    if (eqntype >= 8 && (!gsl_finite(n)||!gsl_finite(gamma0)||!gsl_finite(gamma0p)||!gsl_finite(V0)||!gsl_finite(beta)||!gsl_finite(T0)||!(bn>0)))
      cout<<"Error: Don't have enough input parameters to calculate the density of "<<phasetype<<" using RTpress style EOS."<<endl;

    P /= 1E10;			// convert pressure from microbar to GPa
    // if no temperature information, T should be numeric_limits<double>::quiet_NaN()
    struct EOS_params params = {{P, T}, this};

    if(rho_guess < 0.5 || !gsl_finite(rho_guess) || dP_EOS(rho_guess, &params) < 0)	// rho_guess will be set to negative if it is unknown. Ideal gas doesn't need a rho_guess.
      // if rho_guess is too small, dP/drho can be negative, and the solver may be tricked to the unphysical branch of the solution.
    {
      rho_guess = density(V0) + P/1E3;
    }

    int iter = 0, max_iter = 100;
    const gsl_root_fdfsolver_type *TPL = gsl_root_fdfsolver_newton;
    gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc (TPL);
    gsl_function_fdf FDF;

    double rho = rho_guess, rho0;

    FDF.f = &P_EOS;
    FDF.df = &dP_EOS;
    FDF.fdf = &PdP_EOS;
    FDF.params = &params;

    gsl_root_fdfsolver_set (s, &FDF, rho);

    do
    {
      iter++;

      status = gsl_root_fdfsolver_iterate (s);
      rho0 = rho;
      rho = gsl_root_fdfsolver_root (s);
      if (rho<0.95*rho0)// limit the step size of each iteration to increase stability.
      {
	rho = 0.95*rho0;
	gsl_root_fdfsolver_set (s, &FDF, rho);
      }
      else if (rho>1.05*rho0)
      {
	rho = 1.05*rho0;
	gsl_root_fdfsolver_set (s, &FDF, rho);
      }

      status = gsl_root_test_delta (rho, rho0, 1E-16, rho_eps_rel);
    }
    while (status == GSL_CONTINUE && gsl_finite(rho) && iter < max_iter);

    if (!gsl_finite(rho))
    {
      if (verbose)
	cout<<"Warning: Can't find the density for "<<phasetype<<" at pressure "<<P<<" GPa and temperature "<<T<<" K, initial guessed rho:"<<rho_guess<<". V0, K0, K0p: "<<V0<<' '<<K0<<' '<<K0p<<". Likely no solution exist for this physical condition under the EOS used."<<endl;
      
      gsl_root_fdfsolver_free (s);
      return numeric_limits<double>::quiet_NaN();
    }
    else if (status == GSL_CONTINUE)
    {
      if (verbose)
	cout<<"Warning: Can't find the density for "<<phasetype<<" at pressure "<<P<<" GPa and temperature "<<T<<" K within maximum interation "<<max_iter<<", initial guessed rho:"<<rho_guess<<". V0, K0, K0p: "<<V0<<' '<<K0<<' '<<K0p<<endl;
      
      gsl_root_fdfsolver_free (s);
      return numeric_limits<double>::quiet_NaN();
    }

    gsl_root_fdfsolver_free (s);

    if(thermal_type == 9 && T>T0)	// thermal expansion
      // convert alpha to K^-1
      return rho*exp(-1E-6*pow(1+K0p*P/K0, -xi)*(alpha0*(T-T0)+0.5*alpha1*(sq(T)-sq(T0))));
    else
      return rho;
  }
}

void EOS::printEOS()
// Create a tabulated EOS at  ./tabulated/phasename.txt
// The table from a pressure of 0.1GPa (or P0 if it is larger) to 2000 GPa at the temperature T0 (default 300 K) of the EOS.
{
  string filename = "./tabulated/" + phasetype + ".txt";
  ofstream fout(filename.c_str(),ofstream::trunc);
  if(!fout)
  {
    cout<<"Error: Failed to open file "<<filename<<" to save the EOS data"<<endl;
    return;
  }
  fout<<"Pressure (GPa)\t Density (g/cm^3)"<<endl;
  double rho_guess = 1, P = max(0.1,P0);
  for(int i = 0; i <= 100; i++ )
  {
    if(P > 2000)
      break;
    rho_guess = density(P,T0,rho_guess);
    fout << P << "\t " << rho_guess << endl;
    P *= 1.1;
  }
}

double EOS::entropy(double rho, double T)
// Given the density and temperature, calculate the entropy over n*R, or P V^{7/5} / R = T V^0.4 for ideal gas.
// At either constant P or constant V, entropy always increases with temperature (y decreases with temperature).
{
  double gamma;
  if (eqntype == 6)		// ideal gas,  S ~ nR log(T^(1/(gamma-1))*V) + const.  For better performance and more concise code, here returns T rho^{1-gamma}.
  {
    gamma = this->adiabatic_index();
    return T*pow(rho,1-gamma);
  }

  if (entropy_extern)
    return entropy_extern(rho, T);

  double Theta, V = mmol/rho, x = V/V0;
   
  if (!gsl_finite(V0) || thermal_type < 5 || thermal_type == 8 || thermal_type == 9) // don't have thermal pressure data, or outside the range of the EOS, return negative number
    return -1;
  
  DebyeT(x, gamma, Theta);
  double y = Theta/T;

  double Sth;
  if (Debye_approx == true)		// Debye model
  {
    Sth = 4*gsl_sf_debye_3(y) - 3*log(1-exp(-y));
  }
  else				// Einstein model
  {
    Sth = 3 * (y/(exp(y)-1) - log(1-exp(-y)));
  }

  if (!gsl_finite(e0) || !gsl_finite(g)) // don't have enough information to get Pe
    return Sth;
  else
  {
    return Sth + 3E-6*e0*pow(x,g)*T;
  }
}

double P_EOS(double rho, void *params)
// function of EOS, given a rho (in cgs), return the difference between pressure from EOS and target P (in GPa).  Let this function equals 0 to solve for the correct rho.
{
  struct EOS_params *p = (struct EOS_params *) params;

  double P = p->x[0];
  double T = p->x[1];
  EOS* Phase = p->Phase;
  if (P < Phase->getP0())
  {
    if (verbose)
      cout<<"Warning: Incorrect phase diagram. "<<Phase->getEOS()<<" doesn't exist at pressure "<<P<<" GPa, which is smaller than its transition pressure at "<<Phase->getP0()<<" GPa. The transition pressure applied."<<endl;
    return Phase->getP0() - P;
  }

  return Phase->Press(rho, T) - P;
}

double dP_EOS(double rho, void *params)
{
  gsl_function F;
  double result, abserr;

  F.function = &P_EOS;
  F.params = params;
  gsl_deriv_central(&F, rho, 1E-4, &result, &abserr);
  return result;
}

void PdP_EOS(double rho, void *params, double *P, double *dP)
{
  *P=P_EOS(rho,params);
  *dP=dP_EOS(rho,params);
}

double EOS::gamma0S(double V)
// Grueneisen parameter along the reference adiabat (eq A.3), take volume in cm^3 / mol
{
  double a1 = 6*gamma0;
  double a2 = -12*gamma0 + 36*sq(gamma0) -18*gamma0p;
  double f = fV(V);
  return (2*f + 1) * (a1 + a2*f) / 6 / (1 + a1*f + 0.5*a2*sq(f));
}

double EOS::bV(double V)
// thermal coefficients b(V) in erg/mol (Eq.10), take volume in cm^3 / mol
{
  double sum=0;
  double Vdev = (V/V0-1);
  for (int i=0; i<bn; i++)
    sum += b[i]*pow(Vdev,i);
  return sum;
}

double EOS::bVp(double V)
// derivative of b(V) in erg/cm^3 (microbar) (Eq. B.2), take volume in cm^3 / mol
{
  double sum=0;
  double Vdev = (V/V0-1);
  for (int i=1; i<bn; i++)
    sum += b[i]*i/V0*pow(Vdev,i-1);
  return sum;
}

double EOS::TOS(double V)
// reference adiabat temperature profile (Eq. 7), take volume in cm^3 / mol
{
  double a1 = 6*gamma0;
  double a2 = -12*gamma0 + 36*sq(gamma0) -18*gamma0p;
  double f = fV(V);
  return T0 * sqrt(1+a1*f+0.5*a2*sq(f));
}

double EOS::Cv (double V, double T)
//total heat capacity in erg/mol/K (Eq. B.4), take volume in cm^3 / mol
{
  return bV(V)*fTp(T) + 1.5*n*R;
}

double EOS::Spot(double V, double T)
// potential contribution of entropy in erg/mol/K, take volume in cm^3 / mol
{
  return bV(V)/(beta-1)*fTp(T);
}

double EOS::gamma (double V, double T)
// Grueneisen parameter, take volume in cm^3 / mol
{
  if (!gsl_finite(gamma0))
  {
    cout<<"Error: gamma0 of phase "<<phasetype<<" is not available.  Thermal type is "<<thermal_type<<". Can't calculate Grueneisen parameter at V="<<V<<" T="<<T<<endl;
    return numeric_limits<double>::quiet_NaN();
  }
  
  if (eqntype >= 8)
  {
    // RTpress Eq. 17
    double T_OS = TOS(V);
    return gamma0S(V) * Cv(V,T_OS) / Cv(V,T) + V * bVp(V) / bV(V) * (Spot(V,T)-Spot(V,T_OS)) / Cv(V,T);
  }
  
  if (!gsl_finite(gammainf))	// According to Al'tshuler et al. 1987, gammainf = 2/3 for all elements except alkali elments, for which gammainf = 0.5
    gammainf = 2./3.;
  if (!gsl_finite(beta))	// Altshuler form.
    beta = gamma0 / (gamma0-gammainf);

  return gammainf + (gamma0-gammainf)*pow(V/V0,beta);
}

double EOS::cp(double T)
// specific heat capacity in J/g/K at constant pressure
{
  if (!gsl_finite(cp_a) && cp_b==0 && cp_c==0)
    return numeric_limits<double>::quiet_NaN();

  if (!gsl_finite(cp_a))
    return cp_b*T - cp_c/sq(T);
  else
    return cp_a + cp_b*T - cp_c/sq(T);
}

double EOS::alpha (double P, double T)
// coefficient of thermal expansion in K^-1. Input P in GPa, T in K
{
  if (!gsl_finite(alpha0) && alpha1==0)
    return numeric_limits<double>::quiet_NaN();

  else if (K0p<0 || K0<0 || P<0)
  {
    if (verbose)
      cout<<"Warning: thermal expansion of phase "<<phasetype<<" is not available because some physical parameters are negative, which is nonphysical.  Thermal type is "<<thermal_type<<"."<<endl;
    return numeric_limits<double>::quiet_NaN();
  }
  
  double alphaP0;		// alpha at P=0
  if (!gsl_finite(alpha0))
  {
    alphaP0 = alpha1*T;
  }
  else
  {
    if (T>T0)
      alphaP0 = alpha0 + alpha1*T;
    else				// avoid alpha becomes negative at low temperature
      alphaP0 = alpha0 + alpha1*T0;
  }
  return 1E-6 * alphaP0 * pow(1+K0p*P/K0, -xi);
}
//	#if DEBUG_LEVEL == 2
		
//		std::ofstream output_file;
		//output_file.open("outputTPrho.txt");
		
//	#endif
double oldrho=0, oldP=0; 


		
double EOS::Press(double rho, double T)
// pressure in GPa (Eq. 6, 13, 14) in Wolf&Bower 2018, take density in g/cm^3. For thermal expansion representation, this return the pressure at T0.
{
		

  double P;
  double V = volume(rho);	// volume in cm^3/mol
//  cout << "geteqntype()= " << geteqntype() << endl;
 //cin.ignore(); // Ignore any input from the user
  switch(geteqntype())
  {
  case 0:			// BM3
  case 8:
    P = BM3(rho);
    break;
  case 1:			// BM4
  case 9:
    P = BM4(rho);
    break;
  case 2:			// Vinet
  case 10:
    P = Vinet(rho);
    break;
  case 3:			// Holzapfel
  case 11:
    P = Holzapfel(rho);
    break;
  case 4:			// Keane
  case 12:
    P = Keane(rho);
    break;
  case 6:			// ideal gas
    P = rho*kb*T/(mmol*mp);
    return P;
  case 5:
		//ofstream output_file3;
		//output_file3.open("outputTPrho.txt");
	 //@@@@
		// if (rho<0){
			// cout << "out with P= " << P <<endl;
			// cout << "rho= " << rho <<endl;
			// cout << "T= " << T <<endl;
			// cout << "oldP= " << oldP <<endl;
			// cout << "oldrho= " << oldrho <<endl;
			
		// }
	 P= P_mixEOS(T,rho, 0.5,       1.5, 0.001, 1.0e-3 );
	// #if DEBUG_LEVEL == 1
		// cout << "out with P= " << P <<endl;
		// cout << "Rho_MixEOS(500, P, 0.5); ="<< Rho_MixEOS(500,  P, 0.5)<<endl;
		// cin.ignore(); 
	//	output_file3 << "T= "<<T<<" | P= "<<P<< " | rho= "<<rho <<endl;
		oldrho=rho;
		oldP=P;
		
	// #endif

	return P;
  default:
    cout<<"Error: No such EOS type "<<geteqntype()<<" used in "<<phasetype<<endl;
    P = -1;
    exit(1);
  };

  if (gsl_finite(T) && (getthermal() > 4 || getthermal() == 2)) // Pth = 0 if getthermal == 9
// thermal pressure
    P+= Pth(V,T);

  return P;
}

void printline( double T, double P, double rho) {
	
	 std::fstream outputFile666("output.txt", std::ios::app); 
       if (outputFile666.is_open()) {
         
          outputFile666 <<  "output: T= | " << T << " | P = " << P <<" rho= " << rho <<std::endl; // Add a line to the file
          outputFile666.close(); // Close the file
      } else {
        std::cout << "Failed to open the file." << std::endl;
      }
	  
	  //cin.ignore();  
}


double  P_mixEOS(double T, double rho0, double Xr, double Pguess , double epsilon, double DeltaRho) {
	std::fstream outputFile666("output.txt", std::ios::app);
    outputFile666 <<  "input: T= | " << T << " | P = " << Pguess <<" rho= " << rho0 <<std::endl; // Add a line to the file
    outputFile666.close(); // Close the file	
	if ( (Xr>1)||(T<0)||(rho0<0)){
		cout <<"Error in input of P_mixEOS ! returning default" <<endl;
		cout <<"T= "<<T <<" rho= "<<rho0<< " Xr= "<<Xr <<endl;
		cin.ignore();
		return Pguess;
	}
	//int printDebug=0; // 1 - turn on, 0 -turn off.
    double p =Pguess;
	//double POld= Pguess;
    double rho = Rho_MixEOS(T, p, Xr);
	double toll = rho0*DeltaRho; //tollarance 
    double alpha = 0.1; // step size	
	int errorflag=0;
	double resmin=1e100;//, res_pre_min=	1e100; //huge numbers
	double Pmin=0, P_pre_min=0;
	double resOld= 1e100;
	


	#if DEBUG_LEVEL == 1
	ofstream output_file2;
	ofstream output_file;
		
	//if (printDebug==1){ // create output file
		
		output_file.open("outputDelta.txt");
		output_file2.open("resedue.txt");
		output_file << " P_guess= "<< p  <<"  | tolarance= " <<  toll <<" | rho0= " << rho0 <<"\n"; // write line to file		
	//}
	#endif


	for (int i = 0; i < 1e4; i++){
    
        double drho_dp = (Rho_MixEOS(T, p+epsilon, Xr) - rho) / epsilon; // gradient approximation
		//cout << "T= "<<T<< " p= "<< p << " Xr=" << Xr<< " Rho_MixEOS(T, p, Xr)="<< Rho_MixEOS(T, p, Xr) << " stepsize= "<<alpha * (rho - rho0) / drho_dp<<"\n";
		//cin.ignore();
		
		double stepsize = alpha * (rho - rho0) / drho_dp;	
			while (stepsize>= p) {
				stepsize=stepsize*0.5;
			}
        p = p -  stepsize;// update p
        rho =  Rho_MixEOS(T, p, Xr);
		
		//if (printDebug==1){
		#if DEBUG_LEVEL == 1	
			output_file << i << " | P= "<< p<< " | rhoMix= " << rho << "  | drho_dp= "<< drho_dp << " | abs(rho - rho0) = "<< abs(rho - rho0)  <<" | step = " << alpha * (rho - rho0) / drho_dp<< " tolerance = " <<toll <<"\n"; // write line to file	
			output_file2 << i <<" "<<abs(rho - rho0) << " "<<alpha * (rho - rho0) / drho_dp<<"\n"; 	
		//}
		#endif

		
		if (abs(rho - rho0) < toll){
			#if DEBUG_LEVEL == 1
				cout <<"converged after "<< i <<" steps \n";
			#endif
			break; 
		}
		if ( (abs(rho - rho0) >resOld) ){
			alpha=alpha/5.0; 
			errorflag=errorflag+1;
			#if DEBUG_LEVEL == 1
				cout<< "i= "<<i;
				cout<< " deresing alpha" <<"\n";
			#endif
			if (errorflag>3){
				if (rho >=min(Rho_MixEOS(T, P_pre_min, Xr), Rho_MixEOS(T, Pmin, Xr) ) && rho <= max( Rho_MixEOS(T, P_pre_min, Xr), Rho_MixEOS(T, Pmin, Xr) )){
					printline (T,P_pre_min +(rho0-Rho_MixEOS(T, P_pre_min, Xr) ) * (P_pre_min-Pmin)/(Rho_MixEOS(T, P_pre_min, Xr)- Rho_MixEOS(T, Pmin, Xr)  ), rho );
					return P_pre_min +(rho0-Rho_MixEOS(T, P_pre_min, Xr) ) * (P_pre_min-Pmin)/(Rho_MixEOS(T, P_pre_min, Xr)- Rho_MixEOS(T, Pmin, Xr) );
				}else{
					printline (T,Pmin,  rho );
					return Pmin;
				}
			}
		}
		if ( abs(stepsize )< 0.00005*epsilon ){ // the step becomes zero
			#if DEBUG_LEVEL == 1
				cout<< "i= "<<i;
				cout<< " step too small "<<"\n";
			#endif

			if (rho >=min(Rho_MixEOS(T, P_pre_min, Xr), Rho_MixEOS(T, Pmin, Xr) ) && rho <= max( Rho_MixEOS(T, P_pre_min, Xr), Rho_MixEOS(T, Pmin, Xr) )){
				printline (T,P_pre_min +(rho0-Rho_MixEOS(T, P_pre_min, Xr) ) * (P_pre_min-Pmin)/(Rho_MixEOS(T, P_pre_min, Xr)- Rho_MixEOS(T, Pmin, Xr)  ), rho );
				return P_pre_min +(rho0-Rho_MixEOS(T, P_pre_min, Xr) ) * (P_pre_min-Pmin)/(Rho_MixEOS(T, P_pre_min, Xr)- Rho_MixEOS(T, Pmin, Xr) );
			}else{
				printline (T,Pmin,  rho );
				return Pmin;
			}
		}
		if (p<0){
				cout <<"emergency stop! P<0"<<"\n";
				cin.ignore(); 
		}		
		//POld=p;
		resOld=abs(rho - rho0);
		if ( abs(rho - rho0) <= resmin ){
			//res_pre_min=resmin;
			P_pre_min=Pmin;
			resmin=abs(rho - rho0);
			Pmin=p;

		}
    }
	printline (T,p,  rho );
    return p;
}

//ofstream output_file("outputDelta.txt"); // create output file
	//	output_file << i << " | P= "<< P<< " | Rho_AQUAEOS(T, P)= "<<Rho_AQUAEOS(T, P) <<" | Rho_QEOS(T, P)= "<< Rho_QEOS(T, P)<<" | rhoMix= " << rhoMix  << "  | abs(rhoMix-rho_input)= "<< abs(rhoMix-rho)<< " | epsilon= " << epsilon <<"\n"; // write line to file
		
		
		

//double P_mixEOS(double T, double rho, double Xr, double Pguess, double epsilon, double DeltaRho  ){
	// fuction that caculates P (in Gpa) of mix of rock and water(Xr, Xw=1-Xr) given T and rho in cgs. Pguess is the initial guess. 
	//epsilon is initial guess of timestep (?) and DeltaRho (in percent)
	// @@
//	double Xw, rhoMix,RhoW,RhoR, P,Pfinal, rhoMixFinal;
	// int tstflg=0; //test flag; 1 - prints output
	// int monflag; // 0 -monotonic up, 1 -monotonic down
	// int i;
	
	// Xr=1-Xw;
	// DeltaRho=DeltaRho*rho;

	
	// ofstream output_file("outputDelta.txt"); // create output file
	//ofstream output_file2("tmp.txt"); // create output file
	// P=Pguess;

	// rhoMix= Rho_MixEOS( T, P,  Xr);
	
	// if ((rhoMix-rho)>0 ){
		// monflag=0; //going up
	// }else{
		// monflag=1; //going down
	// }
	
	// if (tstflg==1) {

		// cout <<"input: T=" <<T<<" rho= "<< rho<< " Pguess= "<< Pguess << endl ;
		// cout << "initial: RhoW= "<< RhoW << " RhoR = " <<RhoR << " rhoMix= "<< rhoMix<<endl ;
		// cin.ignore(); // Ignore any input from the user		
		// output_file <<"input: T=" <<T<<" rho= "<< rho<< " Pguess= "<< Pguess << "\n";
		
		
	// }
	
	// for (i = 0; i < 1e4; i++) {

		// rhoMix= Rho_MixEOS( T, P,  Xr);
		//cout << " RhoW= "<< RhoW << " RhoR = " <<RhoR << " rhoMix= "<< rhoMix<<endl ;
	    //cin.ignore(); // Ignore any input from the user	
		
		// output_file << i << " | P= "<< P<< " | Rho_AQUAEOS(T, P)= "<<Rho_AQUAEOS(T, P) <<" | Rho_QEOS(T, P)= "<< Rho_QEOS(T, P)<<" | rhoMix= " << rhoMix  << "  | abs(rhoMix-rho_input)= "<< abs(rhoMix-rho)<< " | epsilon= " << epsilon <<"\n"; // write line to file
		
				
		// output_file2 << i << "  |  Rho_pre_old= "<< Rho_pre_old << " | Rho_old = " << Rho_old << " | rhoMix= " << rhoMix<<endl;
		// if (rhoMix==Rho_pre_old){
				// cout << "Error in  P_mixEOS!  inifinite loop" << endl;
				// cout << " RhoW= "<< RhoW << " RhoR = " <<RhoR << " rhoMix= "<< rhoMix<<endl ;
				// cout << "abs(rhoMix-rho)= "<< abs(rhoMix-rho) << " DeltaRho =" << DeltaRho << endl;
				// cout << "stuck in a loop with P= " << P << endl;
				// cout << "test:"<< endl;
				// cout <<"Rho_AQUAEOS(T, P) = " <<Rho_AQUAEOS( T, P)<< endl;
				// cout <<"Rho_QEOS( T, P) = " <<Rho_QEOS( T, P)<< endl;
				// cout << "1/(Xr/Rho_QEOS( T, P) +Xw/Rho_AQUAEOS( T, P) = " << 1.0/(Xr/Rho_QEOS( T, P) +Xw/Rho_AQUAEOS( T, P) )<<endl;	
				//cin.ignore(); // Ignore any input from the user
				// cout << "epsilon= "<< epsilon <<  endl;
				// if (epsilon > 1e-10){
					// epsilon=0.5*epsilon;
				// }else{
					// DeltaRho= 1.1* DeltaRho;
				// }
				
				
		// }
		
		// if (abs(rhoMix-rho)< DeltaRho) { //found!
			// if (tstflg==1) {
				// cout << "Solution found."<<endl ;
				// cout <<"number of initrations= " <<i<< endl ;
				// cout << " RhoW= "<< RhoW << " RhoR = " <<RhoR << " rhoMix= "<< rhoMix<<endl ;
				// cout << "abs(rhoMix-rho)= "<< abs(rhoMix-rho) << " DeltaRho =" << DeltaRho << endl;
				// cout << "im out with P= " << P << endl;
				// cout << "test:"<< endl;
				// cout <<"Rho_AQUAEOS(T, P) = " <<Rho_AQUAEOS( T, P)<< endl;
				// cout <<"Rho_QEOS( T, P) = " <<Rho_QEOS( T, P)<< endl; 
				// cout << "1/(Xr/Rho_QEOS( T, P) +Xw/Rho_AQUAEOS( T, P) = " << 1.0/(Xr/Rho_QEOS( T, P) +Xw/Rho_AQUAEOS( T, P) )<<endl;
			// }

			// return P;
		// } else {
			// epsilon= 0.0004*(abs(rhoMix-rho)/rho);
			// if (rhoMix> rho) {
				// P -=epsilon;
				// if (monflag==1){
					// break ;
				//}
				
			// } else {
				// P +=epsilon;
				// if (monflag==0){
					// break ;
			//	}
				
			// }
	
			
			

		// }

	// } //end main loop
	
	// if (rhoMix> rho) {
				// Pfinal =P-epsilon;
				
				
	// } else {
				// Pfinal =P+epsilon;
	// }	
		// RhoW=Rho_AQUAEOS(T, Pfinal); 
		// RhoR=Rho_QEOS(T, Pfinal); 
		// rhoMix= Rho_MixEOS( T, Pfinal,  Xr);;
		
		// cout <<" P= "<< P << " Pfinal= " <<Pfinal << " rhoMix = "<< rhoMix << " rhoMixFinal= " << rhoMixFinal <<	endl;
	////	cout <<	" final P= " << (P + (rho -rhoMix)* (Pfinal-P)/(rhoMixFinal -rhoMix) )<<	endl;
		// return (P + (rho -rhoMix)* (Pfinal-P)/(rhoMixFinal -rhoMix));
	
	//cout << "Error in  P_mixEOS! no solution found" << endl;
				// cout << "out of loop with i= "<<i <<endl ;
				// cout << " RhoW= "<< RhoW << " RhoR = " <<RhoR << " rhoMix= "<< rhoMix<<endl ;
				// cout << "abs(rhoMix-rho)= "<< abs(rhoMix-rho) << " DeltaRho =" << DeltaRho << endl;
				// cout << "finised a loop with P= " << P << endl;
				// cout << "test:"<< endl;
				// cout <<"Rho_AQUAEOS(T, P) = " <<Rho_AQUAEOS( T, P)<< endl;
				// cout <<"Rho_QEOS( T, P) = " <<Rho_QEOS( T, P)<< endl;
				// cout << "1/(Xr/Rho_QEOS( T, P) +Xw/Rho_AQUAEOS( T, P) = " << 1.0/(Xr/Rho_QEOS( T, P) +Xw/Rho_AQUAEOS( T, P) )<<endl;	
	
//}


double Rho_MixEOS(double T, double P, double Xr) {
	// funtion that calles two rhows and finds mixed rho
	// formula: 1/rho= Xr/rhor +Xw/rhow 
	    double	RhoW, RhoR,Xw;
		Xw=1-Xr;
		RhoW=Rho_AQUAEOS(T, P); //using guess
		RhoR=Rho_QEOS(T, P); //using guess
		return 1.0 / (Xr/RhoR + Xw/RhoW);	
}

double Rho_AQUAEOS(double Tlin, double Plin) {
	//input in K, Gpa
	//cout << " input: T= " <<Tlin << " P = "<<Plin << endl;
    Plin *= 1e9; // from GPa to Pa
    double P = log10(Plin);
    double T = log10(Tlin);
	
    double rho = NAN;
    //int linFLG = 0;
 //   cout << " input (Aqua units): T= " <<T << " P = "<<P << endl;
    double p00 = 0, p10 = 0, p01 = 0, p20 = 0, p11 = 0, p02 = 0, p30 = 0, p21 = 0, p12 = 0, p03 = 0, p40 = 0, p31 = 0, p22 = 0, p13 = 0, p04 = 0, p50 = 0, p41 = 0, p32 = 0, p23 = 0, p14 = 0, p05 = 0;
	double x, y;
	x=T;
	y=P;
	
    bool aboveLine1 = (P > -26.471 * pow(T, 2) + 150.19 * T - 206.1);
    bool aboveLine2 = (P > -1.0682 * pow(T, 4) + 17.137 * pow(T, 3) - 102.01 * pow(T, 2) + 267.83 * T - 252.75);
    bool aboveLineAB = (P > -13.285 * pow(T, 3) + 102.28 * pow(T, 2) - 252.5 * T + 204.01);
    int reg4flg = 0;

    if (aboveLine1 && (T <= 2.43) && (P <= 8.17)) // A
    {
        p00 = 2.794;
        p10 = 0.1732;
        p01 = 0.0001124;
        p20 = -0.04267;
        p11 = -4.571e-05;
    }
    if ((aboveLine1 || aboveLine2) && (P <= 8.17) && (T > 2.43) && (T <= 2.71)) // A1
    {
        p00 = -4.261;
        p10 = 5.952;
        p01 = -0.01192;
        p20 = -1.216;
        p11 = 0.00115;
        p02 = 0.0009103;
    }
    if ((aboveLineAB > 0) && (P <= 8.17) && (T <= 2.82) && (T > 2.71)) // A2
    {
        p00 = -41.37;
        p10 = 43.67;
        p01 = -3.82;
        p20 = -9.93;
        p11 = 1.291;
        p02 = 0.02265;
    }
    if ((T <= 2.71) && (!aboveLine1)) // B1
    {
        p00 = -3.188;
        p10 = -0.5929;
        p01 = 1.001;
        p20 = -0.07902;
        p11 = -0.0009546;
        p02 = 0.0005481;
    }
if (T <= 2.81 && T > 2.71 && !aboveLineAB) // B2
{
    p00 = 0.868;
    p10 = -3.499;
    p01 = 0.9022;
    p20 = 0.4418;
    p11 = 0.03184;
    p02 = 0.00241;
}

if (T > 2.81 && T <= 3.25 && P < -1.2381 * pow(T, 2) + 9.7711 * T - 10.187) // B3
{
    p00 = -2.464;
    p10 = -1.148;
    p01 = 1.048;
    p20 = 0.02678;
    p11 = -0.01683;
    p02 = 0.0007109;
}

if (T > 3.25 && T <= 3.35 && P < -1.2381 * pow(T, 2) + 9.7711 * T - 10.187) // B4
{
    p00 = 8.679; // (8.553, 8.805)
    p10 = -4.485; // (-4.524, -4.447)
    p01 = -3.669; // (-3.744, -3.594)
    p11 = 1.447; // (1.424, 1.469)
    p02 = 0.4196; // (0.4102, 0.429)
    p12 = -0.1335; // (-0.1363, -0.1306)
    p03 = 0.001257; // (0.001227, 0.001287)
}	

if ( (T <= 3.75) && (T > 3.35) && (P < -1.2381*pow(T, 2) + 9.7711*T - 10.187) ) { //B5

if (P > 3.75) {
                p00 = 0.01413; //(0.01316, 0.01511)
                p10 = -0.296; //(-0.2979, -0.2941)
                p01 = 1.623; //(1.621, 1.625)
                p20 = -0.05786; //(-0.05932, -0.05641)
                p11 = 0.1294; //(0.1282, 0.1306)
                p02 = -0.02996; //(-0.03142, -0.02849)
                p30 = 0.03236; //(0.03027, 0.03445)
                p21 = -0.03754; //(-0.03921, -0.03588)
                p12 = 0.0195; //(0.01782, 0.02119)
                p03 = -0.0009692; //(-0.00304, 0.001102)
                p40 = 0.001306; //(0.000788, 0.001824)
                p31 = -0.026; //(-0.02645, -0.02555)
                p22 = 0.03599; //(0.03554, 0.03644)
                p13 = -0.01107; //(-0.01154, -0.0106)
                p04 = -0.008997; //(-0.009516, -0.008477)
                p50 = -0.005768; //(-0.006362, -0.005174)
                p41 = 0.00668; //(0.006163, 0.007197)
                p32 = 0.003481; //(0.002973, 0.003989)
                p23 = 0.0001127; //(-0.0004005, 0.0006258)
                p14 = -0.001201; //(-0.001736, -0.0006663)
                p05 = -0.005141; //(-0.005727, -0.004556)
       x = (x - 3.553) / 0.1181;
       y = (y - 6.317) / 1.483;
   }
   else
   {
	   
                p00 = -5.341;
                p10 = -0.08695;
                p01 = 1.359;
                p20 = 0.06176;
                p11 = -0.06896;
                p02 = 0.02323;
                p30 = -0.06136;
                p21 = 0.0959;
                p12 = -0.05158;
                p03 = 0.007803;
                p40 = -0.001382;
                p31 = 0.0005141;
                p22 = -0.00809;
                p13 = 0.006997;
                p04 = -0.0005287;
                p50 = 0.005795;
                p41 = -0.01622;
                p32 = 0.02178;
                p23 = -0.0122;
                p14 = 0.000881;
                p05 = 0.0007325;

		x = (x - 3.55) / 0.1183;
		y = (y - 1.372) / 1.373;	   
   }
}
if (T > 3.75 && T <= 3.9 && P < -0.3858 * pow(T, 2) + 4.512 * T - 2.5347) { //B6
	p00 = -2.934;
	p10 = -0.03621;
	p01 = 2.878;
	p20 = 0.001311;
	p11 = -0.01117;
	p02 = 0.0006643;
	p30 = -0.0004874;
	p21 = 0.00704;
	p12 = -0.02717;
	p03 = 0.06056;
	p31 = 0.0004619;
	p22 = -0.003914;
	p13 = 0.004785;
	p04 = 0.01941;
	x = (x - 3.825) / 0.04031;
	y = (y - 4.036) / 2.912; //11.18 and std 0.399
}
if (T > 3.9 && P < -0.3858 * pow(T, 2) + 4.512 * T - 2.5347) { //B7

	p00 = -1.288;
	p10 = -1.962;
	p01 = 5.449;
	p20 = 0.1075;
	p11 = -1.997;
	p02 = -0.4678;
	p21 = 0.2234;
	p12 = 0.2416;
	p03 = -0.007573;
	p22 = -0.03111;
	p13 = -0.0005944;
	p04 = 0.0005344;
	p23 = 0.0008553;
	p14 = -0.000379;
	p05 = 6.461e-05;
}   


if (P > 11.83) //C1
{
    p00 = 4.301; p10 = -0.009416; p01 = 0.3692; p20 = -0.006416; p11 = 0.0113;
    p02 = 0.02318; p30 = -0.002028; p21 = 0.008771; p12 = -0.005377; p03 = -0.001204;
    p40 = -0.0005432; p31 = 0.00269; p22 = -0.004188; p13 = 0.0006328; p04 = -0.0007244;
    p50 = -0.000137; p41 = -4.564e-05; p32 = -0.001235; p23 = 0.0007373; p14 = 0.0002877;
    p05 = 4.814e-05;
    x = (x - 3.5) / 0.8689;
    y = (y - 13.22) / 0.8043;
}
else if (P <= 11.83 && P > (9.003 * pow(T, 4) - 81.988 * pow(T, 3) + 279.25 * pow(T, 2) - 420.9 * T + 245.64) && T < 2.81) //C2
{
    p00 = 3.394; p10 = -0.00158; p01 = 0.1796; p20 = -0.0007509; p11 = 0.002232;
    p02 = 0.0005657; p30 = -0.0002373; p21 = 0.001701; p12 = -0.0005835; p03 = -0.02224;
    p40 = -0.0001095; p31 = 0.0005192; p22 = -0.0007121; p13 = 5.797e-05; p04 = 0.009211;
    p50 = 4.605e-06; p41 = 0.0001306; p32 = -0.0003046; p23 = -0.0002056; p14 = -0.0001364;
    p05 = 0.005902;
    x = (x - 2.381) / 0.2325;
    y = (y - 10.56) / 0.7428;
}

if ((P<=11.83)&&(T<=3.35)&&(P> 5.5834*pow(T,4) - 68.396*pow(T,3) + 310.69*pow(T,2) - 618.37*T + 463.86)&&(T>=2.81)) //C3
{
	p00 = 3.535;
	p10 = -0.002158;
	p01 = 0.09431;
	p20 = -0.0009023;
	p11 = 0.002154;
	p02 = 0.00352;
	p21 = 2.74e-05;
	p12 = -0.001015;
	p03 = 0.008594;
	p22 = 0.0003741;
	p13 = -0.0002965;
	p04 = 0.002108;
	x = (x - 3.051) / 0.1538; // 3.051 and std 0.1538
	y = (y - 11.18) / 0.399; //11.18 and std 0.399
}
if ((P<=11.83)&&(P>10.45)&&(T>3.75)) //C4
{
	p00 = 3.326;
	p10 = -0.1155;
	p01 = 0.1809;
	p20 = -0.02035;
	p11 = 0.03136;
	p02 = -0.00219;
	p30 = 0.001435;
	p21 = 0.001197;
	p12 = 0.0001679;
	p03 = 0.0003041;
	p40 = 0.0005927;
	p31 = -0.002242;
	p22 = 0.0015;
	p13 = -0.0005015;
	x = (x - 4.38) / 0.3608; // mean 4.38 and std 0.3608
	y = (y - 11.14) / 0.3959; // mean 11.14 and std 0.3959
}
if (P<=10.45 && (T>3.9) && (P>= -0.3858*T*T + 4.512*T - 2.5347)) //C5
{
    p00 = 2.864;
    p10 = -0.1374;
    p01 = 0.1495;
    p20 = -0.009803;
    p11 = 0.01471;
    p02 = -0.007024;
    p30 = 0.003783;
    p21 = -0.003452;
    p12 = 0.0008165;
    p03 = 0.0005644;
    p31 = -0.0004699;
    p22 = 0.0003421;
    p13 = 0.0003195;
    p04 = 6.51e-05;
    p32 = -0.0001348;
    p23 = 0.0002483;
    p14 = -0.0003249;
    p05 = 0.000253;

    x = (x - 4.261) / 0.2694; // 4.261 and std 0.2694
    y = (y - 10.05) / 0.2919; // 10.05 and std 0.2919
}

if (P<=10.45 && T<=3.9 && T>=3.75 && P>-0.3858*T*T + 4.512*T - 2.5347) //C6
{
    p00 = 2.921;
    p10 = -0.01531;
    p01 = 0.1826;
    p20 = -0.0004291;
    p11 = 0.003723;
    p02 = -0.01612;
    p30 = -3.244e-06;
    p21 = 8.14e-05;
    p12 = -0.001205;
    p03 = 0.004603;
    p31 = -5.641e-06;
    p22 = 2.432e-05;
    p13 = 0.000668;
    p04 = -0.001128;

    x = (x - 3.823) / 0.04027; // 3.823 and std 0.04027
    y = (y - 9.764) / 0.3989; // 9.764 and std 0.3989
}
if (P <= 11.83 && T < 3.75 && T > 3.35 && P >= -1.2381 * pow(T, 2) + 9.7711 * T - 10.187) { //C7
      p00 = 3.24;
      p10 = -0.01657;
      p01 = 0.265;
      p20 = -0.005132;
      p11 = 0.02098;
      p02 = -0.02605;
      p30 = -0.002369;
      p21 = 0.003013;
      p12 = -0.01378;
      p03 = 0.03458;
      p31 = -0.001151;
      p22 = 0.003333;
      p13 = 0.0002948;
      p04 = -0.003777;
      p32 = 0.00236;
      p23 = -0.002254;
      p14 = 0.00112;
      p05 = -0.002434;
      x = (x - 3.556) / 0.1127; // 3.556 and std 0.1127
      y = (y - 10.36) / 0.8498; // 10.36 and std 0.8498
}
if ((T>3.0)&&(T<=3.35)&&(P>=-1.2381*pow(T,2) + 9.7711*T - 10.187)&&(P< 5.5834*pow(T,4) - 68.396*pow(T,3) + 310.69*pow(T,2) - 618.37*T + 463.86)) { //C8
         p00 = 3.057; //(3.057, 3.058)
         p10 = -0.02485; //(-0.02504, -0.02467)
         p01 = 0.2136; //(0.2134, 0.2139)
         p20 = -0.002074; //(-0.002154, -0.001994)
         p11 = 0.02951; //(0.02933, 0.02969)
         p02 = -0.02923; //(-0.02944, -0.02902)
         p30 = 1.454e-05; //(-7.38e-05, 0.0001029)
         p21 = 0.001657; //(0.001512, 0.001802)
         p12 = -0.01872; //(-0.01901, -0.01844)
         p03 = 0.01692; //(0.01664, 0.01721)
         p31 = 1.108e-05; //(-5.469e-05, 7.685e-05)
         p22 = 0.00111; //(0.001033, 0.001187)
         p13 = 0.002147; //(0.002056, 0.002238)
         p04 = -0.005896; //(-0.005971, -0.00582)
         p32 = 0.0001517; //(7.338e-05, 0.00023)
         p23 = -0.001184; //(-0.001284, -0.001083)
         p14 = 0.001556; //(0.001441, 0.00167)
         p05 = 0.00133; //(0.001247, 0.001413)
       x=(x-3.171)/0.09773; //3.171 and std 0.09773
       y=(y-9.536)/0.7075; //9.536 and std 0.7075     
}
if ((T>2.81)&&(T<=3.0)&&(P>=-1.2381*pow(T,2) + 9.7711*T - 10.187)&&(P< 5.5834*pow(T,4) - 68.396*pow(T,3) + 310.69*pow(T,2) - 618.37*T + 463.86)) { //C9
         p00 = 3; //(3, 3.001)
         p10 = -0.01957; //(-0.02048, -0.01865)
         p01 = 0.1627; //(0.1618, 0.1636)
         p20 = -0.002044; //(-0.002755, -0.001332)
         p11 = 0.01685; //(0.01621, 0.01749)
         p02 = -0.001894; //(-0.002597, -0.00119)
         p30 = -0.001173; //(-0.002223, -0.0001223)
         p21 = 0.004752; //(0.003868, 0.005637)
         p12 = -0.01165; //(-0.01265, -0.01065)
         p03 = 0.006814; //(0.005787, 0.00784)
         p40 = -0.0004157; //(-0.0006718, -0.0001597)
         p31 = -0.001004; //(-0.00124, -0.0007685)
         p22 = 0.0005406; //(0.0002834, 0.0007978)
         p13 = 0.01852; //(0.01822, 0.01881)
         p04 = -0.01631; //(-0.01657, -0.01605)
         p50 = 0.0002879; //(-1.704e-05, 0.0005929)
         p41 = 0.0005869; //(0.000312, 0.0008617)
         p32 = 0.001409; //(0.001115, 0.001703)
         p23 = -0.002036; //(-0.002372, -0.0017)
         p14 = -0.009735; //(-0.01012, -0.009352)
         p05 = 0.008759; //(0.008457, 0.009061)
       x=(x-2.91)/0.05485; //2.91 and std 0.05485
       y=(y-9.017)/0.7393; //9.017 and std 0.7393  
}
if ((P>8.17)&&(P<=9.003*pow(T,4) - 81.988*pow(T,3) + 279.25*pow(T,2) - 420.9*T+ 245.64)&&(T<2.81)) { //triple point region
    int otherflg=1;
    if ((P<=8.18)&&(T<=2.22)) { //add to region a1
               p00 =      -4.261  ;
               p10 =       5.952  ;
               p01 =    -0.01192  ;
               p20 =      -1.216  ;
               p11 =     0.00115  ;
               p02 =   0.0009103  ;
               otherflg=0;
    }
    if ((P< -21.532*pow(T,3) + 148.41*pow(T,2) - 340.08*T + 267.34)&&(T>2.22)&&(T<2.4)) { //add to region a1
               p00 =      -4.261  ;
               p10 =       5.952  ;
               p01 =    -0.01192  ;
               p20 =      -1.216  ;
               p11 =     0.00115  ;
               p02 =   0.0009103  ;
               otherflg=0;
    }
    if ((T>=2.4)&& (P< -607.22*pow(T,2) + 2912.5*T - 3484)) { //add to region a1
               p00 =      -4.261  ;
               p10 =       5.952  ;
               p01 =    -0.01192  ;
               p20 =      -1.216  ;
               p11 =     0.00115  ;
               p02 =   0.0009103  ;
               otherflg=0;
    }
    if ((T>2.4)&&(P< -12.514*pow(T,2) + 67.481*T - 81.378)) { //prob4 - 1
       p00 =       3.046  ;//(3.046, 3.046)
       p10 =    -0.02306  ;//(-0.02308, -0.02305)
       p01 =      0.0523  ;//(0.05229, 0.05232)
       p20 =   -0.003779  ;//(-0.003798, -0.003761)
       p11 =    0.007491  ;//(0.00747, 0.007512)
       p02 =    0.008425  ;//(0.008409, 0.008441)
       p30 =  -0.0006698  ;//(-0.0006796, -0.0006601)
       p21 =    0.003463  ;//(0.003449, 0.003477)
       p12 =   -0.001183  ;//(-0.001198, -0.001167)
       p03 =  -0.0003367  ;//(-0.0003471, -0.0003262)
       p40 =   -0.000185  ;//(-0.000192, -0.000178)
       p31 =    0.000936  ;//(0.0009243, 0.0009478)
       p22 =   -0.001465  ;//(-0.00148, -0.001451)
       p13 =   0.0003629  ;//(0.0003495, 0.0003762)
       p04 =  -4.152e-05  ;//(-4.77e-05, -3.533e-05)
           x=(x- 2.646 )/0.1105; //  2.646 and std 0.1105
           y=(y- 8.838)/0.433 ;  //  8.838 and std 0.433 
           reg4flg=1;
           otherflg=0;

    }
    if ((T>=2.68)&&(reg4flg==0)) { //prob4 - 2
       
       p00 =       3.046  ;//(3.046, 3.046)
       p10 =    -0.02306  ;//(-0.02308, -0.02305)
       p01 =      0.0523  ;//(0.05229, 0.05232)
       p20 =   -0.003779  ;//(-0.003798, -0.003761)
       p11 =    0.007491  ;//(0.00747, 0.007512)
       p02 =    0.008425  ;//(0.008409, 0.008441)
       p30 =  -0.0006698  ;//(-0.0006796, -0.0006601)
       p21 =    0.003463  ;//(0.003449, 0.003477)
       p12 =   -0.001183  ;//(-0.001198, -0.001167)
       p03 =  -0.0003367  ;//(-0.0003471, -0.0003262)
       p40 =   -0.000185  ;//(-0.000192, -0.000178)
       p31 =    0.000936  ;//(0.0009243, 0.0009478)
       p22 =   -0.001465  ;//(-0.00148, -0.001451)
       p13 =   0.0003629  ;//(0.0003495, 0.0003762)
       p04 =  -4.152e-05  ;//(-4.77e-05, -3.533e-05)
           x=(x- 2.646 )/0.1105; //  2.646 and std 0.1105
           y=(y- 8.838)/0.433 ;  //  8.838 and std 0.433
           otherflg=0;
    }
	if (((P>3.9155*pow(T,3) - 26.102*pow(T,2) + 57.565*T - 33.101)&&(T<2.41) ) || ((P>= -12.514*pow(T,2) + 67.481*T - 81.378)&&(T>=2.41)&&(T<2.68))) { //prob5
         p00 =       3.141 ;//(3.141, 3.141)
         p10 =   -0.003927 ;//(-0.003984, -0.003871)
         p01 =     0.00657 ;//(0.006506, 0.006635)
         p20 =   -0.001098 ;//(-0.001136, -0.00106)
         p11 =   0.0005272 ;//(0.0004633, 0.0005911)
         p02 =   0.0005593 ;//(0.0005096, 0.000609)
         p30 =  -0.0001789 ;//(-0.0002144, -0.0001434)
         p21 =    7.68e-05 ;//(6.062e-06, 0.0001475)
         p12 =   2.701e-05 ;//(-5.154e-05, 0.0001056)
         p03 =   2.437e-05 ;//(-1.754e-05, 6.629e-05)
         p31 =  -4.894e-05 ;//(-9.034e-05, -7.544e-06)
         p22 =   -1.85e-05 ;//(-8.825e-05, 5.125e-05)
         p13 =  -8.868e-06 ;//(-6.103e-05, 4.33e-05)
         p04 =   8.046e-06 ;//(-1.125e-05, 2.734e-05)
           x=(x- 2.326 )/0.1203; //  2.326 and std 0.1203
           y=(y- 9.049)/0.1139  ;  // 9.049 and std 0.1139     
           otherflg=0;
    }
    if (otherflg==1) {
         p00 =      -4.112 ;//(-4.633, -3.591)
         p10 =      0.2761 ;//(0.1127, 0.4396)
         p01 =       2.608 ;//(2.447, 2.77)
         p20 =      0.1004 ;//(0.06525, 0.1355)
         p11 =       -0.11 ;//(-0.1373, -0.08273)
         p02 =     -0.3121 ;//(-0.3297, -0.2945)
         p30 =    -0.03736 ;//(-0.04133, -0.03338)
         p21 =     0.01356 ;//(0.01113, 0.01598)
         p12 =    0.003287 ;//(0.001895, 0.00468)
         p03 =     0.01279 ;//(0.01213, 0.01345)        
    }
}

     double rhoLog = p00 + p10*x + p01*y + p20*pow(x,2) + p11*x*y + p02*pow(y,2) + p30*pow(x,3) + p21*pow(x,2)*y  
         + p12*x*pow(y,2) + p03*pow(y,3) + p40*pow(x,4) + p31*pow(x,3)*y + p22*pow(x,2)*pow(y,2)
         + p13*x*pow(y,3) + p04*pow(y,4) + p50*pow(x,5) + p41*pow(x,4)*y + p32*pow(x,3)*pow(y,2)
         + p23*pow(x,2)*pow(y,3) + p14*x*pow(y,4) + p05*pow(y,5);


// cout << "p00: " << p00 << endl;
// cout << "p10: " << p10 << endl;
// cout << "p01: " << p01 << endl;
// cout << "p20: " << p20 << endl;
// cout << "p11: " << p11 << endl;
// cout << "p02: " << p02 << endl;
// cout << "p30: " << p30 << endl;
// cout << "p21: " << p21 << endl;
// cout << "p12: " << p12 << endl;
// cout << "p03: " << p03 << endl;
// cout << "p40: " << p40 << endl;
// cout << "p31: " << p31 << endl;
// cout << "p22: " << p22 << endl;
// cout << "p13: " << p13 << endl;
// cout << "p04: " << p04 << endl;
// cout << "p50: " << p50 << endl;
// cout << "p41: " << p41 << endl;
// cout << "p32: " << p32 << endl;
// cout << "p23: " << p23 << endl;
// cout << "p14: " << p14 << endl;
// cout << "p05: " << p05 << endl;
					

if (p00==0) {
    cout << "error in AQUA" << endl;
    cout <<"P= "<< P << endl;
    cout <<"T= "<< T << endl;
}
  //  cout << " output (Aqua units): rho= " <<rhoLog << endl;

    rho= pow(10,rhoLog); //to lin
	rho=rho*0.001 ;//in gr/cm^3
//    cout << " output (cgs): rho= " <<rho << endl;
//	cin.ignore(); // Ignore any input from the user	
	 return rho;
	 
}

double Rho_QEOS(double T, double P) {
 // input in : K, Gpa
    double p00 = 0, p10 = 0, p01 = 0, p20 = 0, p11 = 0, p02 = 0, p30 = 0, p21 = 0, p12 = 0, p03 = 0, p40 = 0, p31 = 0, p22 = 0, p13 = 0, p04 = 0, rho = NAN;
	double x, y;
 //cout << "geteqntype()= " << geteqntype() << endl;
 // cout << " input: P= " <<P << " T = "<<T << endl;

 
  	T= T*8.61697e-8; // form K to Kev
    T=log10(T); //to Log
	P=P*1e-6; //to PPa
	P=log10(P); //to Log
 // cout << " Attay Units: P= " <<P << " T = "<<T << endl;	
 // cin.ignore(); // Ignore any input from the user	
  if (P >= -4.8 && (P > -0.6836 * pow(T, 4) - 5.7785 * pow(T, 3) - 17.732 * pow(T, 2) - 21.736 * T - 11.52)) {
        if (T <= -2.6) { //A1
            p00 = 2.233;
            p10 = -0.0359;
            p01 = 0.5274;
            p20 = -0.009718;
            p11 = 0.01016;
            p02 = 0.007094;
            p30 = -0.0008563;
            p21 = 0.00142;
            p12 = -0.0008899;
            p03 = -0.003685;
        }
        else { //A2
            p00 = 1.609;
            p10 = -1.037;
            p01 = 1.389;
            p20 = -0.4969;
            p11 = 1.238;
            p02 = -0.2777;
            p30 = -0.07421;
            p21 = 0.6532;
            p12 = -0.394;
            p03 = 0.08212;
            p31 = 0.1228;
            p22 = -0.121;
            p13 = 0.04475;
            p04 = -0.004055;
        }
    }
    else {
        if (P <= -0.6836 * pow(T, 4) - 5.7785 * pow(T, 3) - 17.732 * pow(T, 2) - 21.736 * T - 11.52) { //reg B
            p00 = 1.131;
            p10 = -1.237;
            p01 = 1.005;
            p20 = 0.07896;
            p11 = -0.0125;
        }
        else { //reg C
            p00 = 0.5057;
            p10 = -0.01474;
            p01 = 0.009411;
            p20 = -0.001885;
            p11 = 0.0005793;

            if (P - (-0.6836 * pow(T, 4) - 5.7785 * pow(T, 3) - 17.732 * pow(T, 2) - 21.736 * T - 11.52) < 0.0666) { //small patch - jump to region B
                p00 = 1.253;
                p10 = -1.314;
                p01 = 1.051;
                p20 = 0.09291;
                p11 = -0.02485;
                p02 = 0.003414;
            }
        }
    }

    x = T;
    y = P;
    rho =  p00 + p10*x + p01*y + p20*pow(x,2) + p11*x*y + p02*pow(y,2) + p30*pow(x,3)
            + p21*pow(x,2)*y + p12*x*pow(y,2) + p03*pow(y,3) + p40*pow(x,4) + p31*pow(x,3)*y
            + p22*pow(x,2)*pow(y,2) + p13*x*pow(y,3) + p04*pow(y,4); //in log
//	cout <<  "result:"<< endl;
//	cout << " Attay Units, pre log: rho= " <<rho << endl;
	rho=pow(10,rho); //in gr/cm^3
//	cout << " cgs: rho= " <<rho << endl;
//	  cin.ignore(); // Ignore any input from the user	

if (P <= 0.0937*pow(T,2) + 2.1263*T - 11.9835) {
    rho = NAN;
    // cout << "out" << endl;
}

	return rho;
}  



double S_T(double V, void *params)
// entropy at constant T, volume in cm^3/mol
{
  struct EOS_params *p = (struct EOS_params *) params;
  
  EOS* Phase = p->Phase;
  double rho = Phase->density(V);
  double T = p->x[0];

  return Phase->entropy(rho, T);
}

double S_V(double T, void *params)
// entropy at constant V, volume in cm^3/mol
{
  struct EOS_params *p = (struct EOS_params *) params;

  EOS* Phase = p->Phase;
  double rho = Phase->density(p->x[0]);
  
  return Phase->entropy(rho, T);
}

double EOS::pSpV_T(double V, double T)
// partial S (entropy) partial V at constant T
{
  gsl_function F;
  double result, abserr;
  struct EOS_params params = {{T}, this};
  
  F.function = &S_T;
  F.params = &params;
  gsl_deriv_central(&F, V, 1E-4, &result, &abserr);
  return result;
}

double EOS::pSpT_V(double V, double T)
// partial S (entropy) partial T at constant V
{
  gsl_function F;
  double result, abserr;
  struct EOS_params params = {{V}, this};
  
  F.function = &S_V;
  F.params = &params;
  gsl_deriv_central(&F, T, 1E-2, &result, &abserr);
  return result;
}


double EOS::dTdV_S(double V, double P, double T)
  // adiabatic temperature gradient in K mol/cm^3, take volume in cm^3 / mol, P in GPa
{
  if (thermal_type == 1)	// has external entropy
    // dT/dV_S = - (dS/dV_T) / (dS/dT_V)
  {
    return - pSpV_T(V,T) / pSpT_V(V,T);
  }
  
  if (eqntype == 6)		// ideal gas
  {
    if (!gsl_finite(mmol))
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;

    double gamma = adiabatic_index();
    return (1-gamma)*T/V;
  }

  if (thermal_type == 9)	// thermal expansion
  {
    if (!gsl_finite(cp(300)) || !gsl_finite(alpha(10,300)) || !gsl_finite(mmol))
    {
      cout<<"Error: Information of phase "<<phasetype<<" is not enough to calculate temperature gradient using the thermal expension method."<<endl;
      return 0;
    }
    double a=alpha(P,T);
    // cp to GPa cm^3 g^-1 K^-1
    return (a*T*K0)/(mmol*(sq(a)*T*K0*V/mmol-1E-3*cp(T)));
  }
  
  return -gamma(V,T)*T/V;
}

double P_T(double rho, void *params)
// pressure at constant T in GPa, density in g/cm^3
{
  struct EOS_params *p = (struct EOS_params *) params;
  
  EOS* Phase = p->Phase;
  double T = p->x[0];

  return Phase->Press(rho, T);
}

double P_rho(double T, void *params)
// pressure at constant rho in GPa, density in g/cm^3
{
  struct EOS_params *p = (struct EOS_params *) params;

  EOS* Phase = p->Phase;
  double rho = p->x[0];
  
  return Phase->Press(rho, T);
}

double EOS::pPprho_T(double rho, double T)
  // partial P partial rho at constant T in GPa / g/cm^3
{
  gsl_function F;
  double result, abserr;
  struct EOS_params params = {{T}, this};
  
  F.function = &P_T;
  F.params = &params;
  gsl_deriv_central(&F, rho, 1E-4, &result, &abserr);
  return result;
}

double EOS::pPpT_rho(double rho, double T)
// partial P partial T at constant rho in GPa / K
{
  gsl_function F;
  double result, abserr;
  struct EOS_params params = {{rho}, this};
  
  F.function = &P_rho;
  F.params = &params;
  gsl_deriv_central(&F, T, 1E-2, &result, &abserr);
  // 1E-2 is the step size in conducting derivative. A step size too small will increase the error due to machine error. A value too large may also increase the error and may even crash the code.
  return result;
}


double EOS::dTdm(double m, double r, double rho, double P, double T)
  // adiabatic temperature gradient in K/g, P in cgs
{
  if (eqntype == 6)		// ideal gas
  {
    if (!gsl_finite(mmol))
    {
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;
      return 0;
    }

    double gamma = adiabatic_index();
    return - (gamma-1)*mp*mmol*G*m/(gamma*kb*rho*4*pi*pow(r,4));
  }

  else if (thermal_type == 9)	// thermal expension
  {
    if (!gsl_finite(cp(300)) || !gsl_finite(alpha(10,300)))
    {
      cout<<"Error: Information of phase "<<phasetype<<" is not enough to calculate temperature gradient using the thermal expension method."<<endl;
      return 0;
    }

    return -1E-7*(alpha(P/1E10,T)*T*G*m)/(4*pi*pow(r,4)*rho*cp(T));
  }
  
  double V = volume(rho);
  double dTdV = dTdV_S(V, P/1E10, T);
  if (r<1)		// At the center of the planet where dTdm has a 0/0 limit
  {
    if (m>400 && verbose)
      cout<<"Warning: At the center of of the planet when conducting the first step ODE integration, the material density is "<<m*3./4./pi<<"g/cm^3, which seems to be too high."<<endl;
    return 0;
  }
  return 1E-10*dTdV*G*m/(4*pi*pow(r,4)) / (rho/V*pPprho_T(rho,T) - dTdV*pPpT_rho(rho,T));
}

double EOS::dTdP_S(double P, double T, double &rho_guess)
// partial T partial P along isentrope in K / GPa, given pressure in GPa
{
  if (eqntype == 6)		// ideal gas
  {
    if (!gsl_finite(mmol))
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;

    double gamma = adiabatic_index();
    return - (gamma-1)*T/(gamma*P);
  }

  rho_guess = density(P*1E10, T, rho_guess);
  double V = volume(rho_guess);
  double dTdV = dTdV_S(V, P, T);

  return -V*dTdV/(rho_guess*pPprho_T(rho_guess,T));
}

double P_EOS_S(double rho, void *params)
// function used to solve volume and temperature along adiabatic temperature profile with known temperature gradient.  Read in rho, return the difference between pressure from EOS and target P (in GPa).  Let this function equals 0 to solve for the correct rho.
{
  struct EOS_params *p = (struct EOS_params *) params;

  double P2 = p->x[0];
  double T1 = p->x[1];
  double rho1 = p->x[2];
  double dTdV = p->x[3];
  EOS* Phase = p->Phase;
  double V = Phase->volume(rho);
  double V1 = Phase->volume(rho1);
  
  return Phase->Press(rho, T1+(V-V1)*dTdV) - P2;
}

double dP_EOS_S(double rho, void *params)
{
  gsl_function F;
  double result, abserr;

  F.function = &P_EOS_S;
  F.params = params;
  gsl_deriv_central(&F, rho, 1E-4, &result, &abserr);
  return result;
}

void PdP_EOS_S(double rho, void *params, double *P, double *dP)
{
  *P=P_EOS_S(rho,params);
  *dP=dP_EOS_S(rho,params);
}

double EOS::density(double P1, double T1, double rho, double P2, double &T2)
// Given the pressure (cgs), temperature, density of the previous step, the pressure of the next step, return the temperature and density at the new pressure.  This solver doesn't conserve the entropy well enough. Only used as an approximation in the first integration step from the core of the planet where dTdm has 0/0 limit.
{
  if( !gsl_finite(P1) || !gsl_finite(P2) || !gsl_finite(T1) || !gsl_finite(rho)) // Check if P, the guess of T and rho is infinite or nan due to some error.  Stop code to avoid further error.
  {
    if (verbose)
      cout<<"Warning: Request density for "<<phasetype<<" at infinite/nan value.  P="<<P2/1E10<<" T="<<T1<<" rho_guess="<<rho<<endl;
    T2 = numeric_limits<double>::quiet_NaN();
    return numeric_limits<double>::quiet_NaN();
  }

  if(P2 < 0)
  {
    T2 = P2;
    return P2;
  }

  if (!entropy_extern && thermal_type < 3)
    // Using external density function, but no external entropy function. assuming isothermal.
// don't have thermal pressure data, isothermal applied.
  {
    T2 = T1;
    return density(P2, T2, rho);
  }

  else if (eqntype == 6)		// ideal gas
  {
    if (!gsl_finite(mmol))
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;

    double gamma = adiabatic_index();
    T2 = T1 * pow(P1/P2, (1-gamma)/gamma);
    return P2*mmol*mp/(kb*T2);
  }

  else
  {
    if (!gsl_finite(mmol))
      cout<<"Error: The mean molecular weight of "<<phasetype<<" unknown."<<endl;
	
    if(rho < 0.5 || !gsl_finite(rho))		// rho will be set to negative if it is unknown.
      rho = density(V0) + P2/1E13;
  
    P1 /= 1E10;			// convert pressure from microbar to GPa
    P2 /= 1E10;

    if (eqntype >= 8 && (!gsl_finite(n)||!gsl_finite(gamma0)||!gsl_finite(gamma0p)||!gsl_finite(V0)||!gsl_finite(beta)||!gsl_finite(T0)||!(bn>0)))
    {
      cout<<"Error: Don't have enough input parameters to calculate the density of "<<phasetype<<" using RTpress style EOS."<<endl;
      T2 = numeric_limits<double>::quiet_NaN();
      return numeric_limits<double>::quiet_NaN();
    }
  
    else if (thermal_type >=4 && thermal_type <8 && (!gsl_finite(V0) || !gsl_finite(mmol)||!gsl_finite(gamma0)))
    {
      cout<<"Error: Don't have enough input parameters to calculate the density of "<<phasetype<<" using Debye temperature method."<<endl;
      T2 = numeric_limits<double>::quiet_NaN();
      return numeric_limits<double>::quiet_NaN();
    }
  
    else if (thermal_type ==1 && !entropy_extern)
    {
      cout<<"Error: Don't have user defined external entropy function to calculate the density of "<<phasetype<<"."<<endl;
      T2 = numeric_limits<double>::quiet_NaN();
      return numeric_limits<double>::quiet_NaN();
    }

    else if (thermal_type == 9 && (!gsl_finite(cp_a) || !gsl_finite(alpha0) || !gsl_finite(T0) || !gsl_finite(mmol)))
    {
      cout<<"Error: Information of phase "<<phasetype<<" is not enough to calculate temperature gradient using the thermal expension method."<<endl;
      T2 = numeric_limits<double>::quiet_NaN();
      return numeric_limits<double>::quiet_NaN();
    }
    
    double dTdV = dTdV_S(volume(rho), P1, T1);

    int iter = 0, max_iter = 100;
    const gsl_root_fdfsolver_type *TPL = gsl_root_fdfsolver_newton;
    gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc (TPL);
    gsl_function_fdf FDF;

    double rho2 = rho, rho1=rho;

    struct EOS_params params = {{P2, T1, rho, dTdV}, this};

    FDF.f = &P_EOS_S;
    FDF.df = &dP_EOS_S;
    FDF.fdf = &PdP_EOS_S;
    FDF.params = &params;
  
    int status;
    gsl_root_fdfsolver_set (s, &FDF, rho2);
    
    do
    {
      iter++;

      status = gsl_root_fdfsolver_iterate (s);
      rho1 = rho2;
      dTdV = dTdV_S(volume(rho1), P1, T1);
      params.x[3] = dTdV;
      rho2 = gsl_root_fdfsolver_root (s);
      if (rho2<0.95*rho1)// limit the step size of each iteration to increase stability.
      {
	rho2 = 0.95*rho1;
	gsl_root_fdfsolver_set (s, &FDF, rho2);
      }
      else if (rho2>1.05*rho1)
      {
	rho2 = 1.05*rho1;
	gsl_root_fdfsolver_set (s, &FDF, rho2);
      }

      T2 = T1 + mmol/sq(rho2)*(rho-rho2)*dTdV; // 
      status = gsl_root_test_delta (rho1, rho2, 1E-16, rho_eps_rel);
    }
    while (status == GSL_CONTINUE && gsl_finite(rho) && iter < max_iter);

    if (!gsl_finite(rho2))
    {
      if (verbose)
	cout<<"Warning: Can't find the density for "<<phasetype<<" at pressure "<<P2<<" GPa and temperature "<<T1<<" K, initial guessed density:"<<rho<<". V0, K0, K0p: "<<V0<<' '<<K0<<' '<<K0p<<". Likely no solution exist for this physical condition under the EOS used."<<endl;
      
      gsl_root_fdfsolver_free (s);
      return numeric_limits<double>::quiet_NaN();
    }
    else if (status == GSL_CONTINUE)
    {
      if (verbose)
	cout<<"Warning: Can't find the density for "<<phasetype<<" at pressure "<<P2<<" GPa and temperature "<<T1<<" K within maximum interation "<<max_iter<<", initial guessed density:"<<rho<<". V0, K0, K0p: "<<V0<<' '<<K0<<' '<<K0p<<endl;
      
      gsl_root_fdfsolver_free (s);
      return numeric_limits<double>::quiet_NaN();
    }
    gsl_root_fdfsolver_free (s);
    return rho2;
  }
}
