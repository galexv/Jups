/*Rotation*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "input.h"
#include "opac.h"
#include "atmos.h"
#include "constant.h"
#include "include.h"
#include "nrutil.h"

#include "two_stream.h"

/* --- Global variables ------------------------------------------ */

extern struct Atmos atmos;
extern struct Opac opac;

/* --- Function prototypes --------------------------------------- */

void Locate(int n, double *array, double value, int *ilow);
double Planck(double T, double lambda);

// This is the code that does the two stream calculation
// In progress, by Isaac Malsky
double two_stream(int num_layers, int zero_T_layers, double w0, double g0, double *temperature_array, \
	double *tau_array, double NU, double NU_BIN, double* TMI, double incident_frac);

double lint2D(double x1, double x2, double y1, double y2, double z1,
              double z2, double z3, double z4, double x, double y);
double lint3D(double x1, double x2, double y1, double y2, double z1,
              double z2, double f1, double f2, double f3, double f4,
              double f5, double f6, double f7, double f8, double x,
              double y, double z);
void Angles3d(double ds[], double theta[], double dtheta[], double lat);
double Radius(double R_pl, double ds[]);
double lint(double xa, double ya, double xb, double yb, double x);

/* ------- begin ---------------- RT_Emit -------------------- */

/* Rest frame or Doppler shifted emission spectra for the disk (not including the limb) with clouds turned on or off */

int RT_Emit_3D(double PHASE)
{
    double ***tau_tr_east, ***tau_tr_west, **theta, **dtheta, ***tau_em, ***dtau_em, ***phi_lon_solid, ***theta_lat_solid, ***temperature_3d,
    ***aero_lw_kappa_1, ***aero_lw_tau_1,
    ***aero_lw_kappa_2, ***aero_lw_tau_2,
    ***aero_lw_kappa_3, ***aero_lw_tau_3,
    ***aero_lw_kappa_4, ***aero_lw_tau_4;
    double **intensity, *flux_st, *flux_pl, *flux_tr, *ds, ***dl, **phi,
    *phi_plus_e, *phi_plus_w, *phi_minus_e, *phi_minus_w, **dphi, *theta_lon, *theta_lat, *phi_lon, *dlat_rad, *dlon_rad;
    double R, test, a, b, *lat_rad, kappa_nu_plus_e, kappa_nu_plus_w,
    kappa_nu_minus_e, kappa_nu_minus_w, t_lon_plus_e, t_lon_plus_w,
    t_lon_minus_e, t_lon_minus_w, p_lon_plus_e, p_lon_plus_w,
    p_lon_minus_e, p_lon_minus_w, temperature, pressure, kappa_nu, *lon_rad,
    aero_lw_kappa_interp_1,
    aero_lw_kappa_interp_2,
    aero_lw_kappa_interp_3,
    aero_lw_kappa_interp_4;
    double **I_top, *I_bot, **dkappa_nu;
    int i, j, k, l, m, n, o, c, g, h, ii, kmin, Z;
    double dphid, thetad, dthetad;
    FILE *file;
    FILE *fp;
    double solid;
    double u_vel, v_vel, w_vel, v_los, delta_lam, omega;

    //MALSKY VARIABLES
    double TMI[NTAU];
    double mean_intensity_top;
    double **malsky_intensity;

    double incident_frac;
    double incident_lat;
    double incident_lon;

    char OUTPUT_FILE[200];
    sprintf(OUTPUT_FILE, "%s%06.2f.dat", OUTPUT_PREFIX, PHASE);
    file = fopen(OUTPUT_FILE, "w");
    fp=fopen("inten.txt","w");

    tau_em = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        tau_em[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            tau_em[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    dtau_em = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        dtau_em[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            dtau_em[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    temperature_3d = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        temperature_3d[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            temperature_3d[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    phi_lon_solid = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        phi_lon_solid[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            phi_lon_solid[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    theta_lat_solid = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        theta_lat_solid[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            theta_lat_solid[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    dl = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        dl[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            dl[l][m] = malloc(NTAU*sizeof(double));
        }
    }
    
    /* allocate memory for aero taus and kappas if clouds on */
    
    if(CLOUDS==1){

        /* MgSiO3 */
        aero_lw_kappa_1 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_kappa_1[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_kappa_1[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_lw_tau_1 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_tau_1[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_tau_1[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        
        /* Fe */
        aero_lw_kappa_2 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_kappa_2[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_kappa_2[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_lw_tau_2 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_tau_2[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_tau_2[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        
        /* Al2O3 */
        aero_lw_kappa_3 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_kappa_3[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_kappa_3[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_lw_tau_3 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_tau_3[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_tau_3[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        
        /* MnS */
        aero_lw_kappa_4 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_kappa_4[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_kappa_4[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_lw_tau_4 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_lw_tau_4[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_lw_tau_4[l][m] = malloc(NTAU*sizeof(double));
            }
        }
            printf("clouds: ON\n");
    }
    
    else{
        printf("clouds: OFF\n");
    }
    
    if(DOPPLER==1){
        printf("doppler: ON\n");
    }
    else{
        printf("doppler: OFF\n");
    }
    
    theta = dmatrix(0, NLAT-1, 0, NLON-1);
    dtheta = dmatrix(0, NLAT-1, 0, NLON-1);
    dphi = dmatrix(0, NLAT-1, 0, NLON-1);
    phi = dmatrix(0, NLAT-1, 0, NLON-1);
    
    intensity = dmatrix(0, NLAT-1, 0, NLON-1);
    malsky_intensity = dmatrix(0, NLAT-1, 0, NLON-1);
    I_top = dmatrix(0, NLAT-1, 0, NLON-1);
    
    flux_pl = dvector(0, NLAMBDA-1);
    ds = dvector(0, NTAU-1);
    
    lat_rad = dvector(0, NLAT-1);
    lon_rad = dvector(0, NLON-1);

    dlat_rad = dvector(0, NLAT-1);
    dlon_rad = dvector(0, NLON-1);
    
    /*   Calculate the angular rotational speed omega */
    
    omega = 2.0*PI / (P_ROT*24.0*60.0*60.0);
    
    /*Calculate ds*/
    
    for(j=NTAU-1; j>=0; j--)
    {
        ds[j] = atmos.alt[j-1] - atmos.alt[j];
    }
    ds[0] = ds[1];
    
    /* calculate new aerosol taus and kappas, corrected for wavelength (if clouds on) */
    
    if(CLOUDS==1){
        for(l=0; l<NLAT; l++){
            for(m=0; m<NLON; m++){
                for(j=0; j<NTAU; j++){
                    
                    /* scattering efficiency correction from 5um to 2.3um
                     (PI0, G0, QE calculated from Mie Scattering code of Mischenko, used in Roman 2018) */
                    aero_lw_tau_1[l][m][j] = atmos.aero_sw_tau_1[l][m][j]; //* (1 - PI0_MgSiO3 * G0_MgSiO3 * G0_MgSiO3) * (QE_MgSiO3 / 0.01);
                    aero_lw_tau_2[l][m][j] = atmos.aero_sw_tau_2[l][m][j]; //* (1 - PI0_Fe * G0_Fe * G0_Fe) * (QE_Fe / 0.16);
                    aero_lw_tau_3[l][m][j] = atmos.aero_sw_tau_3[l][m][j]; //* (1 - PI0_Al2O3 * G0_Al2O3 * G0_Al2O3) * (QE_Al2O3 / 0.02);
                    aero_lw_tau_4[l][m][j] = atmos.aero_sw_tau_4[l][m][j]; //* (1 - PI0_MnS * G0_MnS * G0_MnS) * (QE_MnS / 0.02);
                }
            }
        }
        
        for(l=0; l<NLAT; l++){
            for(m=0; m<NLON; m++){
                for(j=0; j<NTAU; j++){
                    
                    aero_lw_kappa_1[l][m][j] = aero_lw_tau_1[l][m][j] / ds[j];
                    aero_lw_kappa_2[l][m][j] = aero_lw_tau_2[l][m][j] / ds[j];
                    aero_lw_kappa_3[l][m][j] = aero_lw_tau_3[l][m][j] / ds[j];
                    aero_lw_kappa_4[l][m][j] = aero_lw_tau_4[l][m][j] / ds[j];
                }
            }
        }
    }
    
    /*Geometry*/
    /*Calculating dl longitude and latitude along line-of-sight*/
    for(l=0;l<NLAT;l++)
    {
        lat_rad[l] = atmos.lat[l] * PI/180.0; /*theta, latitude*/
        for(m=0;m<NLON;m++)
        {
            if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
            {
                lon_rad[m] = atmos.lon[m] * PI/180.0; /*phi, longitude*/
                b = R_PLANET;
                for(j=NTAU-1; j>=0; j--)
                {
                    b += ds[j];
                    
                    dl[l][m][j] = pow(SQ(b) - SQ(R_PLANET * sin(lat_rad[l])) - SQ(R_PLANET * cos(lat_rad[l]) * cos(-PI/2. + lon_rad[m] + PHASE*PI/180.0)), 0.5);
                    
                    phi_lon_solid[l][m][j] = 450. + acos((R_PLANET * cos(-PI/2. + lon_rad[m] + PHASE*PI/180.0) * cos(lat_rad[l]))/(pow(SQ(b)-SQ(R_PLANET * sin(lat_rad[l])), 0.5))) * 180.0/PI;
                    
                    theta_lat_solid[l][m][j] = asin((R_PLANET * sin(lat_rad[l]))/b) * 180.0/PI;
                }
                for(j=0; j<NTAU; j++)
                {
                    if(j!=NTAU-1)
                    {
                        dl[l][m][j] -= dl[l][m][j+1];
                    }
                    else
                    {
                        dl[l][m][j] = dl[l][m][j] - R_PLANET * cos(lat_rad[l]) * sin(lon_rad[m] + PHASE*PI/180.0 - PI/2.);
                    }
                }
            }
        }
    }
    
    /*Calculating solid angle along NLAT*NLON */
    solid = 0.0;
    for(l=0; l<NLAT; l++)
    {
        for(m=0; m<NLON; m++)
        {
            if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
            {
                if(l<NLAT-1 && l>0)
                {
                    theta[l][m] = theta_lat_solid[l][m][0] * PI/180.0;
                    dtheta[l][m] = 0.5*(theta_lat_solid[l+1][m][0] - theta_lat_solid[l-1][m][0])* PI/180.0;
                    dlat_rad[l] = 0.5*(atmos.lat[l+1] - atmos.lat[l-1])* PI/180.0;
                }

                theta[0][m] = theta_lat_solid[0][m][0]* PI/180.0;
                theta[NLAT-1][m] = theta_lat_solid[NLAT-1][m][0]* PI/180.0;
                
                dtheta[0][m] = ( 0.5*(theta_lat_solid[1][m][0]+theta_lat_solid[0][m][0]) + 90)* PI/180.0;
                dtheta[NLAT-1][m] = (90 - 0.5*(theta_lat_solid[NLAT-1][m][0]+theta_lat_solid[NLAT-2][m][0]))* PI/180.0;
                dlat_rad[0] = ( 0.5*(atmos.lat[1]+atmos.lat[0]) + 90.)* PI/180.0;
                dlat_rad[NLAT-1] = (90. - 0.5*(atmos.lat[NLAT-1]+atmos.lat[NLAT-2]))* PI/180.0;
                                
                if(m>0 && m<NLON-1)
                {
                    if(atmos.lon[m]>450.0-PHASE && atmos.lon[m]<630.0-PHASE)
                    {
                        dphi[l][m] = 0.5*(phi_lon_solid[l][m+1][0] - phi_lon_solid[l][m-1][0])* PI/180.0;
                        dlon_rad[m] = 0.5*(atmos.lon[m+1] - atmos.lon[m-1])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 450.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m+1][0] - phi_lon_solid[l][m][0])* PI/180.0;
                        dlon_rad[m] = (atmos.lon[m+1] - atmos.lon[m])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 630.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m][0] - phi_lon_solid[l][m-1][0])* PI/180.0;
                        dlon_rad[m] = (atmos.lon[m] - atmos.lon[m-1])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                }
                
                else if(m==0)
                {
                    if(atmos.lon[m]>450.0-PHASE && atmos.lon[m]<630.0-PHASE)
                    {
                        dphi[l][m] = 0.5*(phi_lon_solid[l][m+1][0] - phi_lon_solid[l][NLON-1][0])* PI/180.0;
                        dlon_rad[m] = 0.5*(atmos.lon[m+1] - atmos.lon[NLON-1])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 450.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m+1][0] - phi_lon_solid[l][m][0])* PI/180.0;
                        dlon_rad[m] = (atmos.lon[m+1] - atmos.lon[m])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 630.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m][0] - phi_lon_solid[l][NLON-1][0])* PI/180.0;
                        dlon_rad[m] = (atmos.lon[m] - atmos.lon[NLON-1])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                }
                
                else if(m==NLON-1)
                {
                    if(atmos.lon[m]>450.0-PHASE && atmos.lon[m]<630.0-PHASE)
                    {
                        dphi[l][m] = 0.5*(phi_lon_solid[l][0][0] - phi_lon_solid[l][m-1][0])* PI/180.0;
                        dlon_rad[m] = 0.5*(atmos.lon[0] - atmos.lon[m-1])* PI/180.0;            
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 450.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][0][0] - phi_lon_solid[l][m][0])* PI/180.0;
                        dlon_rad[m] = (atmos.lon[0] - atmos.lon[m])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 630.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m][0] - phi_lon_solid[l][m-1][0])* PI/180.0;
                        dlon_rad[m] = (atmos.lon[m] - atmos.lon[m-1])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                }
                
        solid += SQ(cos(lat_rad[l]))*cos(lon_rad[m]-PI+PHASE*PI/180.0)*dlat_rad[l]*dlon_rad[m];
            }
        }
    }
    printf("solid %f\n", solid);
        

    // THIS IS TEMPORARY
    for(i=30; i<31; i++)
    //for(i=0; i<NLAMBDA; i++)
    {
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                for(j=0; j<NTAU; j++)
                {
                    tau_em[l][m][j]=0.0;
                }
            }
        }

        /*Optical depth*/
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
                {
                    for(j=0; j<NTAU; j++)
                    {
                        
                        Locate(NLAT, atmos.lat, theta_lat_solid[l][m][j], &o);
                        Locate(NLON, atmos.lon, phi_lon_solid[l][m][j]-PHASE, &c);

                        pressure = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.P_3d[o][c][j], atmos.P_3d[o][c+1][j], atmos.P_3d[o+1][c][j], atmos.P_3d[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                        
                        if(atmos.T_3d[o][c][j] < 100.0 || atmos.T_3d[o][c+1][j] < 100.0)
                        {
                            temperature = 0.0;
                            temperature_3d[l][m][j] = temperature;
                        }
                        else
                        {
                            temperature = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.T_3d[o][c][j], atmos.T_3d[o][c+1][j], atmos.T_3d[o+1][c][j], atmos.T_3d[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            temperature_3d[l][m][j] = temperature;
                        }
                        
                        /* interpolate aero kappas (if clouds on) */
                        if(CLOUDS==1){
                            {
                                aero_lw_kappa_interp_1 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_lw_kappa_1[o][c][j], aero_lw_kappa_1[o][c+1][j], aero_lw_kappa_1[o+1][c][j], aero_lw_kappa_1[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            }
                            {
                                aero_lw_kappa_interp_2 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_lw_kappa_2[o][c][j], aero_lw_kappa_2[o][c+1][j], aero_lw_kappa_2[o+1][c][j], aero_lw_kappa_2[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            }
                            {
                                aero_lw_kappa_interp_3 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_lw_kappa_3[o][c][j], aero_lw_kappa_3[o][c+1][j], aero_lw_kappa_3[o+1][c][j], aero_lw_kappa_3[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            }
                            {
                                aero_lw_kappa_interp_4 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_lw_kappa_4[o][c][j], aero_lw_kappa_4[o][c+1][j], aero_lw_kappa_4[o+1][c][j], aero_lw_kappa_4[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            }
                        }
                    
                        
                        Locate(NTEMP, opac.T, temperature, &g);
                        Locate(NPRESSURE, opac.P, pressure, &h);
                        
                        /* Rotation and wind */
                        if(DOPPLER==1){
                            u_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ew[o][c][j], atmos.vel_ew[o][c+1][j], atmos.vel_ew[o+1][c][j], atmos.vel_ew[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ns[o][c][j], atmos.vel_ns[o][c+1][j], atmos.vel_ns[o+1][c][j], atmos.vel_ns[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            w_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ve[o][c][j], atmos.vel_ve[o][c+1][j], atmos.vel_ve[o+1][c][j], atmos.vel_ve[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                              
                            v_los = u_vel*sin(phi_lon_solid[l][m][j]*PI/180.0) + 
                            v_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*sin(theta_lat_solid[l][m][j]*PI/180.0) - 
                            w_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + 
                            (cos(INPUT_INCLINATION)*(omega*(R_PLANET + atmos.alt[j])*sin(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + 
                            R_VEL*cos((90.0-PHASE)*PI/180.0)));

                            delta_lam = atmos.lambda[i]*v_los/CLIGHT;
                            Locate(NLAMBDA, atmos.lambda, atmos.lambda[i]+delta_lam, &ii);
                            

                            if(temperature < 100.0 || atmos.lambda[i]+delta_lam >= atmos.lambda[NLAMBDA-1] || atmos.lambda[i]+delta_lam + 0.0001e-6 <= atmos.lambda[0])
                            {
                                kappa_nu = 0.0;
                            }
                            else
                            {
                                kappa_nu = lint3D(opac.T[g], opac.T[g+1], opac.P[h],
                                                  opac.P[h+1], atmos.lambda[ii], atmos.lambda[ii+1],
                                                  opac.kappa[ii][h][g],
                                                  opac.kappa[ii][h][g+1],
                                                  opac.kappa[ii][h+1][g],
                                                  opac.kappa[ii][h+1][g+1],
                                                  opac.kappa[ii+1][h][g],
                                                  opac.kappa[ii+1][h][g+1],
                                                  opac.kappa[ii+1][h+1][g],
                                                  opac.kappa[ii+1][h+1][g+1],
                                                  temperature, pressure, atmos.lambda[i]+delta_lam + 0.0001e-6);
                                // MALSKY ADDITION
                                kappa_nu = fabs(kappa_nu);
                            }
                        }

                        /* Wind Only */
                        else if(DOPPLER==2){
                            u_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ew[o][c][j], atmos.vel_ew[o][c+1][j], atmos.vel_ew[o+1][c][j], atmos.vel_ew[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ns[o][c][j], atmos.vel_ns[o][c+1][j], atmos.vel_ns[o+1][c][j], atmos.vel_ns[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            w_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ve[o][c][j], atmos.vel_ve[o][c+1][j], atmos.vel_ve[o+1][c][j], atmos.vel_ve[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            
                            v_los = u_vel*sin(phi_lon_solid[l][m][j]*PI/180.0) + 
                            v_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*sin(theta_lat_solid[l][m][j]*PI/180.0) - 
                            w_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + 
                            R_VEL*cos((90.0-PHASE)*PI/180.0);

                            delta_lam = atmos.lambda[i]*v_los/CLIGHT;
                            Locate(NLAMBDA, atmos.lambda, atmos.lambda[i]+delta_lam, &ii);

                            if(temperature < 100.0 || atmos.lambda[i]+delta_lam >= atmos.lambda[NLAMBDA-1] || atmos.lambda[i]+delta_lam <= atmos.lambda[0])
                            {
                                kappa_nu = 0.0;
                            }
                            else
                            {
                                kappa_nu = lint3D(opac.T[g], opac.T[g+1], opac.P[h],
                                                  opac.P[h+1], atmos.lambda[ii], atmos.lambda[ii+1],
                                                  opac.kappa[ii][h][g],
                                                  opac.kappa[ii][h][g+1],
                                                  opac.kappa[ii][h+1][g],
                                                  opac.kappa[ii][h+1][g+1],
                                                  opac.kappa[ii+1][h][g],
                                                  opac.kappa[ii+1][h][g+1],
                                                  opac.kappa[ii+1][h+1][g],
                                                  opac.kappa[ii+1][h+1][g+1],
                                                  temperature, pressure, atmos.lambda[i]+delta_lam);
                            }
                        }

                        /* Rotation Only */
                        else if(DOPPLER==3){
                            u_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ew[o][c][j], atmos.vel_ew[o][c+1][j], atmos.vel_ew[o+1][c][j], atmos.vel_ew[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ns[o][c][j], atmos.vel_ns[o][c+1][j], atmos.vel_ns[o+1][c][j], atmos.vel_ns[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            w_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ve[o][c][j], atmos.vel_ve[o][c+1][j], atmos.vel_ve[o+1][c][j], atmos.vel_ve[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            
                            v_los = (cos(INPUT_INCLINATION)*(omega*(R_PLANET + atmos.alt[j])*sin(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + 
                            R_VEL*cos((90.0-PHASE)*PI/180.0)));

                            delta_lam = atmos.lambda[i]*v_los/CLIGHT;
                            Locate(NLAMBDA, atmos.lambda, atmos.lambda[i]+delta_lam, &ii);
                            
                            if(temperature < 100.0 || atmos.lambda[i]+delta_lam >= atmos.lambda[NLAMBDA-1] || atmos.lambda[i]+delta_lam <= atmos.lambda[0])
                            {
                                kappa_nu = 0.0;
                            }
                            else
                            {
                                kappa_nu = lint3D(opac.T[g], opac.T[g+1], opac.P[h],
                                                  opac.P[h+1], atmos.lambda[ii], atmos.lambda[ii+1],
                                                  opac.kappa[ii][h][g],
                                                  opac.kappa[ii][h][g+1],
                                                  opac.kappa[ii][h+1][g],
                                                  opac.kappa[ii][h+1][g+1],
                                                  opac.kappa[ii+1][h][g],
                                                  opac.kappa[ii+1][h][g+1],
                                                  opac.kappa[ii+1][h+1][g],
                                                  opac.kappa[ii+1][h+1][g+1],
                                                  temperature, pressure, atmos.lambda[i]+delta_lam);
                            }
                        }
                        
                        /* No Doppler */
                        else{
                            if(temperature < 100.0)
                            {
                                kappa_nu = 0.0;
                            }
                            else
                            {
                                kappa_nu = lint2D(opac.T[g], opac.T[g+1], opac.P[h],
                                                  opac.P[h+1], opac.kappa[i][h][g],
                                                  opac.kappa[i][h][g+1],
                                                  opac.kappa[i][h+1][g],
                                                  opac.kappa[i][h+1][g+1],
                                                  temperature, pressure);
                            }

                        }

                        /* add aerosol opacities to atmosphere, if clouds turned on */
                        if(CLOUDS==1)
                        {
                            kappa_nu += aero_lw_kappa_interp_1;                            
                            kappa_nu += aero_lw_kappa_interp_2;
                            kappa_nu += aero_lw_kappa_interp_3;
                            kappa_nu += aero_lw_kappa_interp_4;
                        }

                        dtau_em[l][m][j] = kappa_nu * dl[l][m][j];
                    }
                }
            }
        }
        
        
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
                {
                    for(j=0; j<NTAU; j++)
                    {
                        if(j==0)
                        {
                            tau_em[l][m][j] = dtau_em[l][m][j];
                        }
                        else if (j!=0)
                        {
                            tau_em[l][m][j] = tau_em[l][m][j-1] + dtau_em[l][m][j];
                        }
                    }
                }
            }
        }

        /*
        // Calculate the intensity of emergent rays at each latitude and longitude
        for(l=0; l<NLAT; l++){
            for(m=0; m<NLON; m++){
                if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
                {

                    intensity[l][m] = Planck(atmos.T_3d[l][m][NTAU-1], atmos.lambda[i]) * exp(-tau_em[l][m][NTAU-1]);
                    I_top[l][m] = Planck(temperature_3d[l][m][0], atmos.lambda[i]);
                }
            }
        }
 

        for(l=0; l<NLAT; l++){
            for(m=0; m<NLON; m++){
                if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE){
                    for(j=0; j<NTAU; j++)
                    {   
                        intensity[l][m] += Planck(temperature_3d[l][m][j], atmos.lambda[i]) * exp(-tau_em[l][m][j]) * dtau_em[l][m][j];
                    }
                
               }
            }
        }
        */

        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
                {
                    // Find the number of 0 elements in tau and temp
                    // This gets passes to the two stream function
                    // This is so the two stream arrays ignore 0 elements
                    kmin=0;
                    while(tau_em[l][m][kmin] < 1e-10 && kmin < NTAU-1)
                    {
                        kmin += 1;
                    }
                    
                    // Grab the incident sunlight fraction
                    // This should be defined in the file being passed in
                    // Incident fractions below 0 correspond to parts
                    // of the planet that solar radiation does not hit
                    incident_frac = 0;
                    if (atmos.incident_frac[l][m][NTAU-1] > 0)
                    {
                        incident_frac = atmos.incident_frac[l][m][NTAU-1];
                    }
                    
                    // This is where the magic happens
                    // We solve the two stream equations in https://ui.adsabs.harvard.edu/abs/1989JGR....9416287T/abstract
                    // Please don't break this
                    malsky_intensity[l][m] = two_stream(NTAU, kmin, 0.0, 0.0, atmos.T_3d[l][m], tau_em[l][m], \
                    CLIGHT / atmos.lambda[i], CLIGHT / atmos.lambda[i] - CLIGHT / atmos.lambda[i+1], TMI, incident_frac);        
                }
            }
        }
        
        
        /*Calculate the total flux received by us*/
        //FILE *fptr = fopen("../Testing/02W0_00g0.txt", "w"); 
        flux_pl[i] = 0.0;
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=450.0-PHASE && atmos.lon[m]<=630.0-PHASE)
                {                    
                    //fprintf(fptr, "%d, %d, %.8e, %.8e\n", l, m, intensity[l][m], malsky_intensity[l][m]);
                    flux_pl[i] += malsky_intensity[l][m] * fabs(SQ(cos(lat_rad[l]))*cos(lon_rad[m]-PI-PHASE*PI/180.0)*dlat_rad[l]*dlon_rad[m]);
                    //flux_pl[i] += intensity[l][m] * fabs(SQ(cos(lat_rad[l]))*cos(lon_rad[m]-PI-PHASE*PI/180.0)*dlat_rad[l]*dlon_rad[m]);
                }
            }
        }
        //fclose(fptr); 
        
        if(i % 100 == 0)
        {
            printf("%d out of %d lines (phase: %06.2f)\n", i, NLAMBDA, PHASE);
        }

        fprintf(file, "%10.8le\t%le\n", atmos.lambda[i], flux_pl[i] * PI/solid);

    }
    
    fclose(file);
    fclose(fp);

    return 0;
}
