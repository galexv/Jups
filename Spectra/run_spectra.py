#!/usr/bin/env python

from shutil import copyfile
from sys import exit
import os
import sys
import pandas as pd
import numpy as np
import run_grid


phases = [90.0]
inclinations = [0.3]
sytem_obliquity = 0
NTAU = 0

planet_name = 'low_grav_clear'


def add_columns(phases, inclinations):
    """For each phase and inclination, add some extra columns and double it
    The file names have to be pretty specific to be run in exotransmit
    phases (list): a list of all the phases to run
    inclinations (list): a list of all the inclinations to run
    """
    input_paths = []              
    output_paths = []
    inclination_strs = []
    phase_strs = []

    columns_to_add = ['aero_sw_tau_1', 'sw_asym_1', 'sw_pi0_1',
                      'aero_sw_tau_2', 'sw_asym_2', 'sw_pi0_2',
                      'aero_sw_tau_3', 'sw_asym_3', 'sw_pi0_3',
                      'aero_sw_tau_4', 'sw_asym_4', 'sw_pi0_4']

    for phase in phases:
        for inc in inclinations:
            phase = str(phase)
            inc = str(inc)
            data_file = 'DATA/init_' + planet_name + '_phase_{}_inc_{}.txt'.format(phase, inc)

            # Read the data file
            df = pd.read_csv(data_file, delim_whitespace=True, names=('lat', 'lon', 'level', 'alt', 'pres', 'temp', 'u', 'v', 'w', 'incident_frac'))
            df = df[(df['lon'] != 360)]

            # Add the 0s
            for column in columns_to_add:
                df[column] = 0

            # Double the data
            double = df.copy()
            double.lon = double.lon + 360.0

            # Get rid of the lon = 360 values
            # I don't know why Eliza's code needs this
            double = double[(double['lon'] != 360)]
            doubled = pd.concat([df, double])

            doubled = doubled.sort_values(by=['lat', 'lon'])
            numpy_df = doubled.to_numpy()

            final_path = 'DATA/Final_' + planet_name + '_phase_{}_inc_{}.txt'.format(phase, inc)
            input_paths.append(final_path)
            inclination_strs.append(inc)
            phase_strs.append(phase)
            np.savetxt(final_path, numpy_df,
            fmt='%-10E  %-10E  %i  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E  %-10E \t')
    
    return input_paths, inclination_strs, phase_strs


def run_exo(input_paths, inclination_strs, phase_strs, doppler_val):
    """
    This runs Eliza's code
    """
    inputs_file = 'input.h'
    output_paths = []

    # The output paths should be similar to the input paths
    # Minus the .dat file extension and saved to OUT/
    for file_path in input_paths:
        output_paths.append('OUT/Spec_' + str(doppler_val) + '_' + file_path[11:-4])

    # Each Run needs to have a specific input.h file
    # With the correct input and output paths
    for i in range(len(input_paths)):
        output_temp = output_paths[i]
        input_temp  = input_paths[i]
        
        # Copy the template for inputs
        try:
            copyfile('template_inputs.h', inputs_file)
        except IOError as e:
            print("Unable to copy file. %s" % e)
            exit(1)
        except:
            print("Unexpected error:", sys.exc_info())
            exit(1)
        
        # Read in the file
        with open(inputs_file, 'r') as file :
            filedata = file.read()

        # Replace the input and output paths
        filedata = filedata.replace("<<output_file>>", "\"" + output_temp + "\"")
        filedata = filedata.replace("<<input_file>>", "\"" + input_temp + "\"")
        filedata = filedata.replace("<<doppler>>", str(doppler_val))
        filedata = filedata.replace("<<inclination>>", inclination_strs[i])
        filedata = filedata.replace("<<phase>>", phase_strs[i])


        # Write the file out again
        with open(inputs_file, 'w') as file:
            file.write(filedata)
        
        # Run Eliza's code
        os.system('make clean')
        os.system('make rt_emission_aerosols.exe') 
        os.system('./rt_emission_aerosols.exe')


#run_grid.run_all_grid(planet_name, phases, inclinations, sytem_obliquity, NTAU)
input_paths, inclination_strs, phase_strs = add_columns(phases, inclinations)

# 0 is off
# 1 is everything
# 2 is Wind only
# 3 is rotation only
#dopplers = [0, 1, 2, 3]
dopplers = [0]
for doppler_val in dopplers:
    run_exo(input_paths, inclination_strs, phase_strs, doppler_val)




