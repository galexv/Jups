#include <stdio.h>

double Planck(double T, double lambda);

double two_stream()
{
  // Here are a chunk of unknown variabels that I really should figure out
  double mu_1 = 0.5773502;

  // This is 1 in Michaels code
  double probably_2pi_mu_1 = 1.0000;
  double mu_0 = 1.0;

  // This is the number of layers
  int NLAYER = 51;

  // These are indexing values
  int J, L, KINDEX;

  // These are just constants
  double bolz_constant   = 1.380649e-23;
  double planck_constant = 6.62607015e-34;

  // These are boundary conditions are values for the flux stuff
  double RSFX;
  double SFCS;
  double EMIS;
  double DIRECT = 0;

  // Frequency Stuff
  double NU;
  double NU_BIN;
  double B_SPECTRAL_DENSITY_VAL;
  double temp_lambda = 0;

  // Scattering and atmosphere parameters
  double W0[51] = {0.7920983991343349, 0.792098418768899, 0.7920983702318101,\
  0.7920983968440899, 0.7920984170315699, 0.7920984323453399,\
  0.792098443962023, 0.7920984527741791, 0.7920984594588829, \
  0.792098464529749, 0.792098468376395, 0.792098471294374, \
  0.7920984735078871, 0.7920984751870089, 0.7920984764607539, \
  0.792098477426987, 0.7920984781599499, 0.7983987264382109, \
  0.8077552840030541, 0.8185675933926501, 0.8306357210395872, \
  0.843746283350503, 0.8583234925153, 0.8734973874694141, \
  0.8874722235231751, 0.900072811490511, 0.9107252787765602, \
  0.9183367659449001, 0.922447681623506, 0.9221620882644591,  \
  0.917492990050582, 0.9100835545057, 0.9015348600388069, \
  0.893449547838393, 0.8866741105572, 0.8804522611815749,\
  0.8734128371968649, 0.8653588930105992, 0.856159631266447,\
  0.845685158729716, 0.8319403025725451, 0.8151284683310831,\
  0.6711648202379829, 0.540089377841582, 0.517763571233066, \
  0.49961811015029794, 0.458273960283873, 0.34007124938820305, \
  0.34007124938820305, 0.34007124938820305, 0.34007124938820305};

  double G0[51] = {0.22168076166565198, 0.221680765809936, 0.221680755565172, \
  0.22168076118224803, 0.221680765443237, 0.22168076867552697, \
  0.22168077112746998, 0.22168077298745897, 0.221680774398406,\
  0.221680775468717, 0.22168077628063199, 0.22168077689653198,\
  0.22168077736374103, 0.221680777718154, 0.221680777987004,\
  0.221680778190948, 0.22168077834565503, 0.23123642759530702,\
  0.24541204997413602, 0.261757159034884, 0.28007854855211606,\
  0.300192214985798, 0.323303335775315, 0.346726866178203,\
  0.36572939220207795, 0.38206457601665994, 0.3936680688270421,\
  0.399353599730723, 0.39727110317405107, 0.38753628471264295,\
  0.378925054069599, 0.38325348969980205, 0.397492604479474,\
  0.409288459759702, 0.415426227651008, 0.419554228800276,\
  0.422859134918534, 0.42569301740957105, 0.427651127381705,\
  0.429119497063642, 0.429882783328733, 0.43001725966865706,\
  0.40874866676711796, 0.356251844291116, 0.341740620713746,\
  0.32524924005012606, 0.257658981581607, 0.0, 0.0, 0.0, 0.0};

  // Info aboutt the profile of the atmosphere
  double TAUS[51]={0.09412313301298277, 0.12407845372340899, \
  0.06944412605637637, 0.0915451865019821, 0.120680058363644, \
  0.15908729952225797, 0.20971790398896198, 0.27646203950437703,\
  0.364447945906697, 0.4804359599195471, 0.6333379409837621,\
  0.834902007860287, 1.1006151972711002, 1.45089339973522,\
  1.91265000050298, 2.52136375103796, 3.3238047545326603,\
  4.39955408098399, 5.82918952737921, 7.715002390531371,\
  10.184131200357902, 13.391965419272198, 17.5210077296356,\
  22.7199868415078, 29.0726583607662, 36.505110340223204,\
  44.487874881804295, 52.0422106252831, 57.690958842437496,\
  59.799270749042996, 58.0501889100439, 54.3170317302161,\
  50.450795142558604, 47.25159852399221, 45.211286170163405,\
  43.6575583750587, 42.4078270414836, 41.3991429925596,\
  40.81856678494329, 40.373963151118005, 40.567464177028896,\
  41.270829823854896, 26.2627218235042, 20.219076205705303,\
  22.614229707154003, 26.3749944502078, 27.341591963452203,\
  25.307281524531103, 3.3614944042572, 3.979014806724, 55.7019444577607};

  double TEMPS[51] = {928.7, 929.1, 929.2, 929.4, 929.7, 930.0, 930.4,\
  931.0, 931.7, 932.7, 933.9, 935.6, 937.7, 940.5, 944.2,\
  948.9, 955.1, 962.9, 972.9, 985.4, 1001.0, 1020.2, 1043.3,\
  1070.8, 1102.7, 1139.0, 1179.1, 1222.4, 1267.7, 1313.2,\
  1357.3, 1398.0, 1433.7, 1463.6, 1488.3, 1509.8, 1531.0,\
  1555.1, 1584.4, 1620.4, 1664.5, 1717.7, 1781.0, 1855.3,\
  1941.3, 2039.5, 2150.6, 2274.8, 2412.6, 2564.3, 2722.3};

  // How to generate the matrix from the toon paper
  double y1[51];
  double y2[51];
  double LAMBDA[51];
  double GAMMA[51];
  double temp_e_val[51];
  double e1[51];
  double e2[51];
  double e3[51];
  double e4[51];

  // Matrix Coefficients
  double A[102];
  double B[102];
  double D[102];
  double E[102];
  
  // Solution Coefficients
  double AS[102];
  double DS[102];
  double X[102];
  double Y[102];

  // More temperatoray matrix soluton stuff
  double temp_gamma_val[51];
  double CP[51];
  double CPB[51];
  double CM[51];
  double CMB[51];

  // Planck Function Stuff, and the slope
  // These are strane and hard
  double B0[51];
  double B1[51];
  
  // Why we came all this way
  double FNET[51];
  double TMI[51];

  double Bnu, twohnu3_c2, hc_Tkla;
  double temp_val_1, temp_val_2;

  NU = 1.0e9;
  NU_BIN = 1.0;

  RSFX = 0.0;
  EMIS = 1.0;

  
  // HERE WE DDINE LAYER PROPERTIES FOLLOWING GENERAL SCHEME
  // OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
  // NEEDED FOR MATRIX.

  for(J=0; J<NLAYER; J++)
  {
    y1[J]    =  2.0 - (W0[J] * (1.0 + G0[J]));
    y2[J]    =  W0[J] * (1.0 - G0[J]);

    LAMBDA[J]    =  sqrt(fabs(pow(y1[J], 2.0) - pow(y2[J], 2.0)));
    GAMMA[J]  =  y2[J] / (y1[J] + LAMBDA[J]);
    temp_e_val[J]   =  exp(-LAMBDA[J] * TAUS[J]);

    e1[J]   =  1.0 + GAMMA[J] * temp_e_val[J];  //e1                          
    e2[J]   =  1.0 - GAMMA[J] * temp_e_val[J]; //e2                      
    e3[J]   =  GAMMA[J] + temp_e_val[J];       //e3                        
    e4[J]   =  GAMMA[J] - temp_e_val[J];       //e4
  }


  J = 0;
  for(L=1; L<2*NLAYER -1; L+=2)
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

  A[2*NLAYER-1] = e1[NLAYER-1] - RSFX * e3[NLAYER-1];
  B[2*NLAYER-1] = e2[NLAYER-1] - RSFX * e4[NLAYER-1];
  D[2*NLAYER-1] = 0.0;


  // This is the part of the code that solves for the blackbody stuff
  for(J=0; J<NLAYER; J++)
  {
    if(0 >= J-1)
      KINDEX = 0;
	  else
      KINDEX = J-1;

    temp_val_1 = (2.0 * planck_constant * (NU * NU * NU)) / (CLIGHT * CLIGHT);
    temp_val_2 = exp(planck_constant * NU / (bolz_constant * TEMPS[J])) - 1.0;
    B_SPECTRAL_DENSITY_VAL = temp_val_1 * (1.0 /  temp_val_2);

    B0[J] = B_SPECTRAL_DENSITY_VAL * NU_BIN;
    B1[J] = (B0[J] - B0[KINDEX]) / TAUS[J];

    if (TAUS[J] <= 1.0e-6)
      B1[J] = 0.0;
  }

  SFCS = EMIS * B0[0] * PI;

  // This solves for the C values in the toon code
  // I don't like this part but I think it works
  for(J=0; J<NLAYER; J++)
  {
    if(0 > J-1)
	    KINDEX = 0;
	  else
      KINDEX = J-1;

    temp_gamma_val[J]   = 1.0 / (y1[J] + y2[J]);

    CP[J]  = (B0[KINDEX] + B1[J] * temp_gamma_val[J]) * probably_2pi_mu_1;
    CPB[J] = CP[J] + B1[J] * TAUS[J] * probably_2pi_mu_1;

    CM[J]  = (B0[KINDEX] - B1[J] * temp_gamma_val[J]) * probably_2pi_mu_1;
    CMB[J] = CM[J] + B1[J] * TAUS[J] * probably_2pi_mu_1;
  }


  J = 0;
  for(L=1; L<2*NLAYER; L+=2)
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

  E[0]           = -CM[0];
  E[2*NLAYER-1]  = SFCS + RSFX * CMB[NLAYER-1] - CPB[NLAYER-1];
  DS[2*NLAYER-1] = E[2*NLAYER-1] / B[2*NLAYER-1];
  AS[2*NLAYER-1] = A[2*NLAYER-1] / B[2*NLAYER-1];


  //********************************************
  //*     WE SOLVE THE TRIDIAGONAL EQUATIONS   *
  //********************************************

  for(L=2; L<2*NLAYER+1; L++)
  {
    X[2*NLAYER-L]  = 1.0 / (B[2*NLAYER-L] - D[2*NLAYER-L] * AS[2*NLAYER-L+1]);
    AS[2*NLAYER-L] = A[2*NLAYER-L] * X[2*NLAYER-L];
    DS[2*NLAYER-L] = (E[2*NLAYER-L]-D[2*NLAYER-L] * DS[1+2*NLAYER-L]) * X[2*NLAYER-L];
  }
  
  Y[0] = DS[0];

  for(L=1; L<2*NLAYER; L++)
  {
  	Y[L] = DS[L] - AS[L] * Y[L-1];
  }

  //***************************************************************
  //  CALCULATE LAYER CODFICIENTS, NET FLUX AND MEAN INTENSITY
  //***************************************************************

  for(J=1; J<NLAYER+1; J++)
  {
  	FNET[J]  = Y[2*J-2] * (e1[J-1]-e3[J-1]) + Y[2*J-2] * (e2[J-1]-e4[J-1]) + CPB[J-1] - CMB[J-1] - DIRECT;

    TMI[J] = (1.0 / mu_1) * (Y[2*J-2]*(e1[J-1] + e3[J-1]) + Y[2*J-1] * (e2[J-1]+e4[J-1]) \
             + CPB[J-1] + CMB[J-1]) +  (DIRECT / mu_0);
  }

  return 0;
}


