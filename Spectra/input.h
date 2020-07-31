/*----------------------- input.h ------------------------------------

Defines input values and files for 3-D emission spectra

---------------------------------------------------------------------- */

#ifndef __INPUT_H__
#define __INPUT_H__

/* I/O SETTINGS. */

/* File names */
#define OUTPUT_PREFIX "OUT/Spec_0_low_grav_clear_phase_90.0_inc_0.3"      /* output name */
#define T_P_3D_FILE "DATA/Final_low_grav_clear_phase_90.0_inc_0.3.txt"         /* input file */

/* Output settings */
#define N_PHASE 1                          /* Number of phases [96 max; lon grid in increments of 3.75] */
#define DOPPLER 0                /* 0:Off; 1:On */
#define CLOUDS 0                           /* 0:Off; 1:On */

/* Grid settings */
#define NTAU 60                            /* Number of altitude points in grid      */
#define NLAT  50                           /* Number of latitude points in 3-D  grid */
#define NLON  99                           /* Number of longitude points in 3-D grid */

#define NTEMP 30                           /* Number of temperature points in grid   */
#define NPRESSURE 13                       /* Number of pressure points in grid   [13/17]   */
#define NLAMBDA 6000                       /* Number of wavelength points in grid [4616/2598]   */

/* Planet parameters */
#define INPUT_INCLINATION 0.3  /* Planet inclination in radians            */
#define INPUT_PHASE 90.0              /* Planet inclination in degrees           */
#define G 12.9324                          /* Planet surface gravity                 */
#define R_PLANET 1.287e+08                 /* Planet radius at base of atmosphere      */

#define ORB_SEP 1.0e+10                    // This is some distance
#define STELLAR_TEMP 6000.0                // Stellar Blackbody temperature
#define R_STAR 1.296e+09                   /* Stellar radius                         */
#define P_ROT 4.6171                       /* Rotation period in days (= P_ORB for tidally locked planet)    */
#define R_VEL 0.0                          /* Radial Velocity                        */
#define MU 2.36                            /* Mean molecular weight                  */
#define FORMAT 2                           /* FORMAT=1 -> small opacity table        */
                                           /* FORMAT=2 -> large opacity table        */

/* Aerosol properties (calculated by the Mischenko Mie code) */

#define PI0_MgSiO3 0.99
#define G0_MgSiO3 0.14
#define QE_MgSiO3 0.07

#define PI0_Fe 0.71
#define G0_Fe -0.13
#define QE_Fe 1.25

#define PI0_Al2O3 0.74
#define G0_Al2O3 0.15
#define QE_Al2O3 0.12

#define PI0_MnS 1.00
#define G0_MnS 0.35
#define QE_MnS 0.56


/* Opacities for hi-res spectra */

#define CHEM_FILE "DATA/eos_solar_doppler_2016_cond.dat"
#define CH4_FILE "DATA/opacCH4_hires_Brogi.dat"
#define CO2_FILE "DATA/opacCO2_hires_Brogi.dat"
#define CO_FILE "DATA/opacCO_hires_Hayley.dat"
#define H2O_FILE "DATA/opacH2O_hires_Brogi_HITEMP.dat"
#define NH3_FILE "DATA/opacNH3_hires_Brogi.dat"
#define O2_FILE "DATA/opacO2_hires_Brogi.dat"
#define O3_FILE "DATA/opacO3_hires_Brogi.dat"

/* Opacities for low-res spectra */
/*
#define CHEM_FILE "DATA/eos_solar_doppler_2016_cond.dat"
#define CH4_FILE "DATA/opacCH4.dat"
#define CO2_FILE "DATA/opacCO2.dat"
#define CO_FILE "DATA/opacCO.dat"
#define H2O_FILE "DATA/opacH2O.dat"
#define NH3_FILE "DATA/opacNH3.dat"
#define O2_FILE "DATA/opacO2.dat"
#define O3_FILE "DATA/opacO3.dat"
*/
 
#endif /* !__INPUT_H__ */

/* ------- end ---------------------------- input.h  ----------------- */


