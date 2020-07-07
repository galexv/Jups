#include <stdio.h>

double Planck(double T, double lambda);

double two_stream(int NLAYER, int kmin, double w0_val, double g0_val, \
                 double *temperature_array, double *tau_array, \
                 double NU, double NU_BIN, double* TMI)
{
  double mu_1 = 1.0;
  double mu_0 = 1.0;  // Needs to be adjusted for the solid angle stuff

  // These are indexing values
  int J, L, KINDEX, Z;
  //int NEW_NLAYER;
  int NEW_NLAYER = 51;

  // These are just constants
  double bolz_constant = 1.380649e-23;
  double h_constant    = 6.62607015e-34;

  // These are boundary conditions are values for the flux stuff
  double EMIS = 1.;
  double RSFX = 1.;
  double redistribution_param = 1.0;

  double SFCS;
  double DIRECT;
  double STELLAR_BB;
  double BB_TOP_OF_ATM;
  double BB_BOTTOM_OF_ATM;

  // Frequency Stuff
  double B_SPECTRAL_DENSITY_VAL;

  // Scattering and atmosphere parameters
  //double W0[NLAYER - kmin];
  //double G0[NLAYER - kmin];

  // Scattering and atmosphere parameters
  //double TEMPS[NLAYER - kmin];
  //double TAUCS[NLAYER - kmin];
  //double TAULS[NLAYER - kmin];

  // How to generate the matrix from the toon paper
  double y1[NLAYER - kmin];
  double y2[NLAYER - kmin];
  double LAMBDAS[NLAYER - kmin];
  double GAMMA[NLAYER - kmin];
  double temp_e_val[NLAYER - kmin];
  double e1[NLAYER - kmin];
  double e2[NLAYER - kmin];
  double e3[NLAYER - kmin];
  double e4[NLAYER - kmin];

  // Matrix Coefficients
  double A[2 * (NLAYER - kmin)];
  double B[2 * (NLAYER - kmin)];
  double D[2 * (NLAYER - kmin)];
  double E[2 * (NLAYER - kmin)];
  
  // Solution Coefficients
  double AS[2 * (NLAYER - kmin)];
  double DS[2 * (NLAYER - kmin)];
  double X[2 * (NLAYER - kmin)];
  double Y[2 * (NLAYER - kmin)];

  // More temperatoray matrix soluton stuff
  double temp_gamma_val[NLAYER - kmin];
  double CP[NLAYER - kmin];
  double CPB[NLAYER - kmin];
  double CM[NLAYER - kmin];
  double CMB[NLAYER - kmin];

  // Planck Function Stuff, and the slope
  // These are strane and hard
  double B0[NLAYER - kmin];
  double B1[NLAYER - kmin];
  double FNET[NLAYER - kmin];

  double Bnu, twohnu3_c2, hc_Tkla;
  double temp_val_1, temp_val_2;

  double TWO_STREAM_INTENSITY;
  double SOURCE_INTENSITY;

  
  // These variables solve the source function technique
  // I really should figure out a way to check these
  // All these vary for each layer, but as long as I do
  // Everything in a single loop it's fine
  double SOURCE_G[NLAYER - kmin];
  double SOURCE_H[NLAYER - kmin];
  double SOURCE_J[NLAYER - kmin];
  double SOURCE_K[NLAYER - kmin];

  double ALPHA_1[NLAYER - kmin];
  double ALPHA_2[NLAYER - kmin];
  double SIGMA_1[NLAYER - kmin];
  double SIGMA_2[NLAYER - kmin];
  double SOURCE_Y1[NLAYER - kmin];
  double SOURCE_Y2[NLAYER - kmin];
  double source_temp[NLAYER - kmin];
  double mu = 1.0;

  double INTENSITY_DOWN[NLAYER - kmin];
  double INTENSITY_UP[NLAYER - kmin];

  /*
  NEW_NLAYER = NLAYER - kmin;

  for (J=0; J<NEW_NLAYER; J++)
  {
    W0[J] = w0_val;
    G0[J] = g0_val;
  
    TEMPS[J] = temperature_array[J+kmin];
    TAUCS[J] = tau_array[J+kmin];

    TAULS[J] =  -1.0 * (tau_array[J+kmin] - tau_array[J+kmin + 1]);
  }

  // Data Sanitation
  if (TAULS[1] < 1e-8)
  {
    TAULS[1] = TAULS[2];
  }

  if (TAULS[0] < 1e-8)
  {
    TAULS[0] = TAULS[1];
  }

  if (TAULS[NEW_NLAYER-1] < 1e-8)
  {
    TAULS[NEW_NLAYER-1] = TAULS[NEW_NLAYER-2];
  }

  if (TEMPS[0] < 1)
  {
    TEMPS[0] = TEMPS[1];
  }
  */


  double W0[] = {7.817670187228165E-002,5.979540194933853E-002,\
  3.390544549138463E-002,0.000000000000000E+000,\
  0.000000000000000E+000,0.000000000000000E+000,\
  0.000000000000000E+000,7.541926930821227E-002,8.610388742978484E-002,\
  8.114495479445624E-002,8.114503856306433E-002,8.114510210816728E-002,\
  8.114515031212381E-002,8.114518687860246E-002,8.114521461713511E-002,\
  8.114523565897663E-002,8.114525162085297E-002,8.114526372917864E-002,\
  8.114527291428536E-002,8.114527988190275E-002,8.114528516738183E-002,\
  8.114528917682824E-002,8.114529221830458E-002,8.114529452550048E-002,\
  8.598749266923190E-002,9.442147489112089E-002,0.105814531639559,\
  0.121776416165926,0.144302039555417,0.178533838444239,0.229830533607487,\
  0.300075311627056,0.394466324831345,0.511311766327278,0.629838681235809,\
  0.728208632988234,0.797517784277552,0.843341307465869,0.877649130407956,\
  0.899337921747103,0.906316319343882,0.903345976275421,0.892383990523738,\
  0.873074095623763,0.842198182958891,0.796579853091006,0.739415708498399,\
  0.684823856196835,0.633282310127787,0.350731998560003,0.156325746000884};

  double G0[] = {0.110326882158668 ,0.110326956918972 ,0.110326772109971,\
  0.110326873438366 ,0.110326950303938 ,0.110327008612508 ,\
  0.110327052844123 ,0.110327086397258 ,0.110327111849928 ,\
  0.110327131157760 ,0.110327145804253 ,0.110327156914757 ,\
  0.110327165342937 ,0.110327171736366 ,0.110327176586279 ,\
  0.110327180265315 ,0.110327183056148 ,0.112254085593167 ,\
  0.114959453958986 ,0.118358383332087 ,0.122183827464402 ,\
  0.126155103807160 ,0.130313245599292 ,0.133748779520810 ,\
  0.135375962548969 ,0.134183699412485 ,0.129238005738520 ,\
  0.121600030011815 ,0.115855879644724 ,0.122830902173550 ,\
  0.154889553165107 ,0.225947922855962 ,0.315123504542173 ,\
  0.358428805263750 ,0.388659181988818 ,0.408227725579857 ,\
  0.418568543084830 ,0.419002832244987 ,0.410708458415535 ,\
  0.402779488738229 ,0.408189397139313 ,0.421549423462755 ,\
  0.421108712660989 ,0.375638623840812 ,0.379455427694435 ,\
  0.381757089883631 ,0.345678203649185 ,0.000000000000000E+000,\
  0.000000000000000E+000,0.000000000000000E+000,0.000000000000000E+000};

  double TAULS[] = {3.785206799298454E-003,4.369965972797778E-003,\
  2.595438703931929E-003,3.421454547729317E-003,4.510355499969544E-003,\
  5.945806517880663E-003,7.838099495109633E-003,1.033262746366466E-002,\
  1.362105576772974E-002,1.795604853882313E-002,2.367068197091029E-002,\
  3.120403600113891E-002,4.113493071541886E-002,5.422639959241263E-002,\
  7.148431665711115E-002,9.423468212109742E-002,0.124225504698618,\
  0.164346850475370,0.217802839656924,0.288312940957087,\
  0.381698208868834,0.504313298985020,0.664601927142552,\
  0.870448149355705, 1.13363591367261, 1.46837177270562,\
  1.89879903080063, 2.46800061035889, 3.26133011659987, \
  4.44504276493956, 6.24956799366401, 8.99077758138895, \
  11.5088199165136, 13.7880334308601, 14.9228878340194, \
  15.4293229311154, 16.0014607231307, 17.2370839545222, \
  19.2839577291430, 21.6324646844230, 23.9839555160984, \
  26.5834619303397, 23.1690453010524, 24.7407793725123, \
  31.5058061322877, 40.5853307621816, 51.4183088548791, \
  65.3288657406461, 86.1202174848536, 113.528557025352, \
  139.338065787836};


  double TEMPS[] = {928.7, 929.1, 929.2, 929.4, 929.7, 930.0, 930.4,\
  931.0, 931.7, 932.7, 933.9, 935.6, 937.7, 940.5, 944.2,\
  948.9, 955.1, 962.9, 972.9, 985.4, 1001.0, 1020.2, 1043.3,\
  1070.8, 1102.7, 1139.0, 1179.1, 1222.4, 1267.7, 1313.2,\
  1357.3, 1398.0, 1433.7, 1463.6, 1488.3, 1509.8, 1531.0,\
  1555.1, 1584.4, 1620.4, 1664.5, 1717.7, 1781.0, 1855.3,\
  1941.3, 2039.5, 2150.6, 2274.8, 2412.6, 2564.3, 2722.3};


  double TAUCS[] = {3.785206799298454E-003,4.369965972797778E-003,\
  2.595438703931929E-003,3.421454547729317E-003,4.510355499969544E-003,\
  5.945806517880663E-003,7.838099495109633E-003,1.033262746366466E-002,\
  1.362105576772974E-002,1.795604853882313E-002,2.367068197091029E-002,\
  3.120403600113891E-002,4.113493071541886E-002,5.422639959241263E-002,\
  7.148431665711115E-002,9.423468212109742E-002,0.124225504698618,\
  0.164346850475370,0.217802839656924,0.288312940957087,\
  0.381698208868834,0.504313298985020,0.664601927142552,\
  0.870448149355705, 1.13363591367261, 1.46837177270562,\
  1.89879903080063, 2.46800061035889, 3.26133011659987, \
  4.44504276493956, 6.24956799366401, 8.99077758138895, \
  11.5088199165136, 13.7880334308601, 14.9228878340194, \
  15.4293229311154, 16.0014607231307, 17.2370839545222, \
  19.2839577291430, 21.6324646844230, 23.9839555160984, \
  26.5834619303397, 23.1690453010524, 24.7407793725123, \
  31.5058061322877, 40.5853307621816, 51.4183088548791, \
  65.3288657406461, 86.1202174848536, 113.528557025352, \
  139.338065787836};









  // Calculate the intensity at the top of the atmosphere
  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[NEW_NLAYER-1])) - 1.0;
  BB_TOP_OF_ATM = temp_val_1 * (1.0 / temp_val_2);

  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[NEW_NLAYER-1])) - 1.0;
  BB_BOTTOM_OF_ATM = temp_val_1 * (1.0 / temp_val_2);

  // Calculate the flux at the top of the atmosphere from the star
  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * STELLAR_TEMP)) - 1.0;
  STELLAR_BB = temp_val_1 * (1.0 / temp_val_2);

  // Define the energy from other sources (eq 37 and 38, Toon)
  if (NU > 430.0e12)
  {
    DIRECT = redistribution_param * PI * STELLAR_BB * pow(R_STAR / ORB_SEP, 2.0);
    SFCS = RSFX * DIRECT * mu_0 * exp(-(TAUCS[NEW_NLAYER-1] + TAULS[NEW_NLAYER-1]) / mu_0);
  }
  else
  {
    SFCS = EMIS * PI * BB_BOTTOM_OF_ATM;
    DIRECT = 0.0;
  }

  // HERE WE FIND LAYER PROPERTIES FOLLOWING GENERAL SCHEME
  // OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
  // NEEDED FOR MATRIX.

  for(J=0; J<NEW_NLAYER; J++)
  {
    y1[J]    =  2.0 - (W0[J] * (1.0 + G0[J]));
    y2[J]    =  W0[J] * (1.0 - G0[J]);

    LAMBDAS[J]    =  sqrt(fabs(pow(y1[J], 2.0) - pow(y2[J], 2.0)));
    GAMMA[J]  =  y2[J] / (y1[J] + LAMBDAS[J]);
    temp_e_val[J]   =  exp(-LAMBDAS[J] * TAULS[J]);

    e1[J]   =  1.0 + GAMMA[J] * temp_e_val[J];  //e1                          
    e2[J]   =  1.0 - GAMMA[J] * temp_e_val[J]; //e2                      
    e3[J]   =  GAMMA[J] + temp_e_val[J];       //e3                        
    e4[J]   =  GAMMA[J] - temp_e_val[J];       //e4
  }


  J = 0;
  for(L=1; L<2*NEW_NLAYER -1; L+=2)
  {
    // HERE ARE THE EVEN MATRIX ELEMENTS
    A[L]   =  e2[J+1] * e1[J]   - e4[J+1] * e3[J];
    B[L]   =  e2[J+1] * e2[J]   - e4[J+1] * e4[J];
    D[L]   =  e1[J+1] * e4[J+1] - e3[J+1] * e2[J+1];
  
    // HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP. 
    A[L+1] =  e2[J]   * e3[J]   - e1[J]   * e4[J];
    B[L+1] =  e1[J+1] * e1[J]   - e3[J+1] * e3[J]; 
    D[L+1] =  e3[J]   * e4[J+1] - e1[J]   * e2[J+1];
    J = J + 1;
   }

  // HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
  // BEGINNING OF THE TRIDIAGONAL SOLUTION DDINITIONS. I ASSUME
  // NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
  A[0] = 0.0;
  B[0] = e1[0];
  D[0] = -e2[0];

  A[2*NEW_NLAYER-1] = e1[NEW_NLAYER-1] - RSFX * e3[NEW_NLAYER-1];
  B[2*NEW_NLAYER-1] = e2[NEW_NLAYER-1] - RSFX * e4[NEW_NLAYER-1];
  D[2*NEW_NLAYER-1] = 0.0;


  // This is the part of the code that solves for the blackbody stuff
  for(J=0; J<NEW_NLAYER; J++)
  {
    if(0 >= J-1)
      KINDEX = 0;
	  else
      KINDEX = J-1;
    
    temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
    temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[J])) - 1.0;

    B_SPECTRAL_DENSITY_VAL = temp_val_1 * (1.0 / temp_val_2);

    //B0[J] = B_SPECTRAL_DENSITY_VAL * NU_BIN;
    B0[J] = B_SPECTRAL_DENSITY_VAL;
    B1[J] = (B0[J] - B0[KINDEX]) / TAULS[J];
  }

  // You need this because sometimes TAU=0 and you don't want a NAN
  // THIS IS ME GUESSING!!!!
  //if (B1[0] < 1e-15)
  //  B1[0] = B1[1];
  B1[0] = 0;
  //B1[1] = -1.0 * B1[1];

  // This solves for the C values in the toon code
  // I don't like this part but I think it works
  for(J=0; J<NEW_NLAYER; J++)
  {
    if(0 > J-1)
      KINDEX = 0;
    else
      KINDEX = J-1;

    temp_gamma_val[J]   = 1.0 / (y1[J] + y2[J]);

    CP[J]  = (B0[KINDEX] + B1[J] * temp_gamma_val[J]) * 2.0 * PI * mu_1;
    CPB[J] = CP[J] + B1[J] * TAULS[J] * 2.0 * PI * mu_1;

    CM[J]  = (B0[KINDEX] - B1[J] * temp_gamma_val[J]) * 2.0 * PI * mu_1;
    CMB[J] = CM[J] + B1[J] * TAULS[J] * 2.0 * PI * mu_1;
  }

  J = 0;
  for(L=1; L<2*NEW_NLAYER; L+=2)
  {
    // HERE ARE THE EVEN MATRIX ELEMENTS
    E[L]   = (CP[J+1] - CPB[J]) * e2[J+1] - (CM[J+1] - CMB[J]) * e4[J+1];

    // HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
    E[L+1] = e3[J] * (CP[J+1] - CPB[J]) + e1[J] * (CMB[J] - CM[J+1]);
    J = J + 1;
  }

  // HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
  // BEGINNING OF THE TRIDIAGONAL SOLUTION DDINITIONS. I ASSUME NO
  // DIFFUSE RADIATION IS INCIDENT AT THE TOP.

  E[0] = -CM[0];
  E[2*NEW_NLAYER-1]  = SFCS + RSFX * CMB[NEW_NLAYER-1] - CPB[NEW_NLAYER-1];
  DS[2*NEW_NLAYER-1] = E[2*NEW_NLAYER-1] / B[2*NEW_NLAYER-1];
  AS[2*NEW_NLAYER-1] = A[2*NEW_NLAYER-1] / B[2*NEW_NLAYER-1];



  //********************************************
  //*     WE SOLVE THE TRIDIAGONAL EQUATIONS   *
  //********************************************

  for(L=2; L<2*NEW_NLAYER+1; L++)
  {
    X[2*NEW_NLAYER-L]  = 1.0 / (B[2*NEW_NLAYER-L] - D[2*NEW_NLAYER-L] * AS[2*NEW_NLAYER-L+1]);
    AS[2*NEW_NLAYER-L] = A[2*NEW_NLAYER-L] * X[2*NEW_NLAYER-L];
    DS[2*NEW_NLAYER-L] = (E[2*NEW_NLAYER-L]-D[2*NEW_NLAYER-L] * DS[1+2*NEW_NLAYER-L]) * X[2*NEW_NLAYER-L];


    //printf("%d %.8e %.8e %.8e %.8e\n", 2*NEW_NLAYER-L, E[2*NEW_NLAYER-L], D[2*NEW_NLAYER-L], DS[1+2*NEW_NLAYER-L], X[2*NEW_NLAYER-L]);

  }
  
  Y[0] = DS[0];
  for(L=1; L<2*NEW_NLAYER; L++)
  {
    Y[L] = DS[L] - AS[L] * Y[L-1];
  }

  //***************************************************************
  //  CALCULATE LAYER CODFICIENTS, NET FLUX AND MEAN INTENSITY
  //***************************************************************

  // I reverse the list here because that's how it wants it in Eliza's code
  for(J=1; J<NEW_NLAYER+1; J++)
  {
    FNET[NEW_NLAYER-J] = Y[2*J-2] * (e1[J-1]-e3[J-1]) + Y[2*J-2] \
                 * (e2[J-1]-e4[J-1]) + CPB[J-1] - CMB[J-1] - DIRECT;

    TMI[NEW_NLAYER-J] = (1.0 / mu_1) * (Y[2*J-2]*(e1[J-1] + e3[J-1]) + Y[2*J-1] * (e2[J-1]+e4[J-1]) \
             + CPB[J-1] + CMB[J-1]) +  (DIRECT / mu_0);

    TMI[NEW_NLAYER-J] = TMI[NEW_NLAYER-J] / (4.0 * PI);

  }

  //Pretty sure this is right now
  for(J=1; J<NEW_NLAYER+1; J++)
  {
    SOURCE_Y1[J-1] = Y[2*J-2];
    SOURCE_Y2[J-1] = Y[2*J-1];
  }
  
  for(J=0; J<NEW_NLAYER; J++)
  {
    SOURCE_G[J] = (SOURCE_Y1[J] + SOURCE_Y2[J]) * (1.0/mu_1 - LAMBDAS[J]);
    SOURCE_H[J] = (SOURCE_Y1[J] - SOURCE_Y2[J]) * GAMMA[J] * (1.0/mu_1 + LAMBDAS[J]);
    SOURCE_J[J] = (SOURCE_Y1[J] + SOURCE_Y2[J]) * GAMMA[J] * (1.0/mu_1 + LAMBDAS[J]);
    SOURCE_K[J] = (SOURCE_Y1[J] - SOURCE_Y2[J]) * GAMMA[J] * (1.0/mu_1 - LAMBDAS[J]);

    source_temp[J] = (1.0 / (y1[J] + y2[J])) - mu_1;
    ALPHA_1[J]     = 2.0 * PI * (B0[J] + (B1[J] * source_temp[J]));
    ALPHA_2[J]     = 2.0 * PI * B1[J];

    SIGMA_1[J] = 2.0 * PI * (B0[J] - (B1[J] * source_temp[J]));
    SIGMA_2[J] = 2.0 * PI * B1[J];
  }


  INTENSITY_DOWN[0] = BB_TOP_OF_ATM * exp(-TAULS[0]) + \
                      SOURCE_J[0]/(LAMBDAS[0] + 1.0) * (1.0 - exp(-TAULS[0]*LAMBDAS[0]+1.0)) + \
                      SOURCE_K[0]/(LAMBDAS[0] - 1.0) * (exp(-TAULS[0]) - exp(-TAULS[0]*LAMBDAS[0])) + \
                      SIGMA_1[0] * (1.0 - exp(-TAULS[0])) + \
                      SIGMA_2[0] * (exp(-TAULS[0]) + TAULS[0] + 1.0);


  INTENSITY_UP[NEW_NLAYER-1] = 2.0 * BB_BOTTOM_OF_ATM * EMIS * PI;

  // Do the downward intensity first
  for(J=1; J<NEW_NLAYER; J++)
  {
    INTENSITY_DOWN[J] = INTENSITY_DOWN[J-1] * exp(-TAULS[J]) + \
                    SOURCE_J[J]/(LAMBDAS[J] + 1.0) * (1.0 - exp(-TAULS[J]*LAMBDAS[J]+1.0)) + \
                    SOURCE_K[J]/(LAMBDAS[J] - 1.0) * (exp(-TAULS[J]) - exp(-TAULS[J]*LAMBDAS[J])) + \
                    SIGMA_1[J] * (1.0 - exp(-TAULS[J])) + \
                    SIGMA_2[J] * (exp(-TAULS[J]) + TAULS[J] + 1.0);
    
  }

  // Calculate the upward intensity next
  printf("\n\n\n");
  for(Z=1; Z<NEW_NLAYER; Z++)
  {
    J = NEW_NLAYER - Z - 1;
    
    INTENSITY_UP[J] = INTENSITY_UP[J+1] * exp(-TAULS[J+1]) + \
                      SOURCE_G[J+1]/(LAMBDAS[J+1]-1.0)*(exp(-TAULS[J+1])-exp(-TAULS[J+1]*LAMBDAS[J+1])) + \
                      SOURCE_H[J+1]/(LAMBDAS[J+1]+1.0) * (1.0 - exp(-TAULS[J+1] * (LAMBDAS[J+1] + 1.0))) + \
                      ALPHA_1[J+1] * (1.0 - exp(-TAULS[J+1])) + \
                      ALPHA_2[J+1] * (1.0 - ((TAULS[J+1] + 1.0) * (exp(-TAULS[J+1]))));

    //printf("%d %.8e %.8e \n", J, INTENSITY_UP[NEW_NLAYER-1], INTENSITY_UP[J]);
  }

  TWO_STREAM_INTENSITY = TMI[0];
  SOURCE_INTENSITY     = INTENSITY_UP[0];

  printf("\n\n\n\n\n");
  for(J=0; J<NEW_NLAYER; J++)
  {
    printf("%.8e, \n", INTENSITY_UP[J]);
  //  printf("%.8e, \n", TMI[J]);  
  }

  // Define the energy from other sources (eq 37 and 38, Toon)
  if (NU > 430.0e12)
  {
    return TWO_STREAM_INTENSITY;
  }
  else
  {
    return SOURCE_INTENSITY;
  }

  
}


