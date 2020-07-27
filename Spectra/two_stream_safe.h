#include <stdio.h>

double Planck(double T, double lambda);

double two_stream(int NLAYER, int kmin, double w0_val, double g0_val, \
                 double *temperature_array, double *tau_array, \
                 double NU, double NU_BIN, double* TMI, double incident_frac)
{
  double mu_1 = 1.0;
  double mu_0 = 1.0;  // Needs to be adjusted for the solid angle stuff
  double mu = 1.0;

  // These are indexing values
  int J, L, KINDEX, Z, M;
  int NEW_NLAYER;

  // These are just constants
  double bolz_constant = 1.380649e-23;
  double h_constant    = 6.62607015e-34;

  // These are boundary conditions are values for the flux stuff
  double EMIS = 0.;
  double RSFX = 0.;
  double redistribution_param = 1.0;

  double SFCS;
  double DIRECT;
  double STELLAR_BB;
  double BB_TOP_OF_ATM;
  double BB_BOTTOM_OF_ATM;

  // Frequency Stuff
  double B_SPECTRAL_DENSITY_VAL;

  // Scattering and atmosphere parameters
  double W0[NLAYER - kmin];
  double G0[NLAYER - kmin];

  // Scattering and atmosphere parameters
  double TEMPS[NLAYER - kmin];
  double TAUCS[NLAYER - kmin];
  double TAULS[NLAYER - kmin];

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
  double B0[NLAYER - kmin];
  double B1[NLAYER - kmin];
  double FNET[NLAYER - kmin];
  double Bnu, twohnu3_c2, hc_Tkla;
  double temp_val_1, temp_val_2;

  // The Source Function Variables
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

  // Upward and downwards intensities
  double INTENSITY_DOWN[NLAYER - kmin];
  double INTENSITY_UP[NLAYER - kmin];

  // Top layer values
  double TWO_STREAM_INTENSITY;
  double SOURCE_INTENSITY;
  double RUNNING_SOURCE;

  // The number of layers
  // Sometimes the inputs are bad and the top
  // Of the layers need to be cut off
  NEW_NLAYER = NLAYER - kmin;

  for (J=0; J<NEW_NLAYER; J++)
  {
    W0[J] = w0_val;
    G0[J] = g0_val;

    TEMPS[J] = temperature_array[J+kmin];
    TAUCS[J] = tau_array[J+kmin];

    TAULS[J] =  -1.0 * (tau_array[J+kmin] - tau_array[J+kmin + 1]);
  }

  //**********************************************
  //* Data Sanitation (Sorry, kinda gross code)  *
  //**********************************************

  for (J=0; J<10; J++)
  {
    if (TAULS[9-J] < 1e-10)
    {
      TAULS[9-J] = TAULS[10-J] - fabs(TAULS[10-J] - TAULS[11-J]);
    }
  }

  for (J=0; J<10; J++)
  {
    if (TEMPS[9-J] < 10)
    {
      TEMPS[9-J] = TEMPS[10-J] - fabs(TEMPS[10-J] - TEMPS[11-J]);
    }
  }

  if (TAULS[NEW_NLAYER-1] < 1e-10)
  {
    TAULS[NEW_NLAYER-1] = TAULS[NEW_NLAYER-2] + abs(TAULS[NEW_NLAYER-2] - TAULS[NEW_NLAYER-3]);
  }


  
  // Calculate the intensity at the top of the atmosphere
  temp_val_1 = (2.0 * h_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
  temp_val_2 = exp(h_constant * NU / (bolz_constant * TEMPS[NEW_NLAYER-1])) - 1.0;
  BB_TOP_OF_ATM = temp_val_1 * (1.0 / temp_val_2);

  // Calculate the intensity at the bottom of the atmosphere
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
  	//redistribution param is unused
    DIRECT = incident_frac * PI * STELLAR_BB * pow(R_STAR / ORB_SEP, 2.0);
    SFCS = RSFX * DIRECT * mu_0 * exp(-(TAUCS[NEW_NLAYER-1] + TAULS[NEW_NLAYER-1]) / mu_0);

    //DIRECT=0;
    //SFCS=0;
  }
  else
  {
    SFCS = EMIS * PI * BB_BOTTOM_OF_ATM;
    DIRECT = 0.0;

    //SFCS = 0.0;
    //DIRECT = 0.0;

  }

  // HERE WE FIND LAYER PROPERTIES FOLLOWING GENERAL SCHEME
  // OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
  // NEEDED FOR MATRIX.

  for(J=0; J<NEW_NLAYER; J++)
  {
    y1[J]    =  (2.0 - (W0[J] * (1.0 + G0[J])));
    y2[J]    =  (W0[J] * (1.0 - G0[J]));

    LAMBDAS[J]    =  sqrt(fabs(pow(y1[J], 2.0) - pow(y2[J], 2.0)));
    GAMMA[J]  =  y2[J] / (y1[J] + LAMBDAS[J]);
    temp_e_val[J]   =  exp(-LAMBDAS[J] * TAULS[J]);

    e1[J]   =  1.0 + GAMMA[J] * temp_e_val[J];  //e1                          
    e2[J]   =  1.0 - GAMMA[J] * temp_e_val[J];  //e2                      
    e3[J]   =  GAMMA[J] + temp_e_val[J];        //e3                        
    e4[J]   =  GAMMA[J] - temp_e_val[J];        //e4
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

    //B_SPECTRAL_DENSITY_VAL = temp_val_1 * (1.0 / temp_val_2);
    B_SPECTRAL_DENSITY_VAL = temp_val_1 * (1.0 / temp_val_2) * 1e15; //Magic 1e15 (Dont delete)


    B0[J] = B_SPECTRAL_DENSITY_VAL;
    B1[J] = (B0[J] - B0[KINDEX]) / TAULS[J];
  }

  // The very top of the atmosphere is isothermal
  B1[0] = 0;

  // This solves for the C values in the toon code
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
  }

  //********************************************
  //*    Source Function Technique Solution    *
  //********************************************
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
                    SIGMA_2[J] * (exp(-TAULS[J]) + TAULS[J] - 1.0);
  }

  // Calculate the upward intensity next
  for(Z=1; Z<NEW_NLAYER; Z++)
  {
    J = NEW_NLAYER - Z - 1;
    INTENSITY_UP[J] = INTENSITY_UP[J+1] * exp(-TAULS[J+1]) + \
                      SOURCE_G[J+1]/(LAMBDAS[J+1]-1.0)*(exp(-TAULS[J+1])-exp(-TAULS[J+1]*LAMBDAS[J+1])) + \
                      SOURCE_H[J+1]/(LAMBDAS[J+1]+1.0) * (1.0 - exp(-TAULS[J+1] * (LAMBDAS[J+1] + 1.0))) + \
                      ALPHA_1[J+1] * (1.0 - exp(-TAULS[J+1])) + \
                      ALPHA_2[J+1] * (1.0 - ((TAULS[J+1] + 1.0) * (exp(-TAULS[J+1]))));
  }
  
  
  //printf("\n");
  //for(J=1; J<NEW_NLAYER; J++)
  //{
  //  printf("%.8e, %.8e, \n", (TMI[NEW_NLAYER - J - 1] / (1e15 * 4.0 * PI)), (INTENSITY_UP[J] / (1e15 * 2.0 * PI)));
  //  //printf("%.8e, %.8e, %.8e, %.8e \n", (TMI[NEW_NLAYER - J - 1] / (1e15 * 4.0 * PI)), (INTENSITY_UP[J] / (1e15 * 2.0 * PI)), TEMPS[J], TAULS[J]);
  //}
  //printf("\n\n\n\n\n\n\n");


  TWO_STREAM_INTENSITY = TMI[NEW_NLAYER-1] / (4.0 * PI);
  SOURCE_INTENSITY     = INTENSITY_UP[0] / (2.0 * PI);
  
  // Define the energy from other sources (eq 37 and 38, Toon)
  if (NU > 430.0e12)
  {
    return TWO_STREAM_INTENSITY / 1e15;
  }
  else
  {
  	return TWO_STREAM_INTENSITY / 1e15;
  } 
}