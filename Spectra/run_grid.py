import pandas as pd
import numpy as np
from scipy.interpolate import Rbf


def run_all_grid(planet_name, phases,inclinations, sytem_obliquity):
    print ('Running the regridding')
    # This regrids all the stuff

    def df_to_txt(file, df):
        np.savetxt(file, df.values, fmt='%-8E  %-8E  %i  %-8E  %-8E  %-8E  %-8E  %-8E  %-8E %-8E\t')

    planet_file = '../Planets/' + planet_name + '.txt'

    for phase in phases:
        for inc in inclinations:
            df = pd.read_csv(planet_file,
                             delim_whitespace=True,
                             names=('lat', 'lon', 'level', 'alt', 'pres',
                                    'temp', 'u', 'v'), index_col='lat')

            # Get rid of any place where the longitude is equal to 360
            # Not needed by Eliza Code
            df = df[(df['lon'] != 360)].reset_index()

            def get_incident_flux(df, sytem_obliquity):
                """Add a new column corresponding to the fraction of intensity"""
                
                # Find the incident fraction depending on latitude and longitude
                df['incident_frac'] = np.sin(df.lat * np.pi / 180.0) * np.sin(sytem_obliquity) + \
                                      np.cos(df.lat * np.pi / 180.0) * np.cos(sytem_obliquity) * np.cos(df.lon * np.pi / 180.0)

                # Anywhere not on the circle illuminated set to 0
                #df.incident_frac = df.incident_frac.mask(df.incident_frac < 0.0, 0.0)
                return df

            def phase_rotation(df, phase):
                """ Rotate the planet a certain phase and rollover longitude"""
                # Rotate by the given phase
                df.lon = df.lon + phase

                # Rollover the longitude values if they are too big
                df.lon = df.lon.mask(df.lon >= 360.0, df.lon - 360.0)
                return df


            def wind_rot(df):
                """ Calculate the new u wind speed"""
                u = df.u
                v = df.v
                lon = df.lon
                lat = df.lat

                obs_theta_degree = 180.0 - phase
                obs_theta = obs_theta_degree * (np.pi / 180.0)
                obs_theta = 0

                # Convert to radians
                phi = lat * (np.pi / 180.0)
                theta = lon * (np.pi / 180.0)

                # Get the prime cordinates
                phi_prime = np.arcsin(np.cos(inc) * np.sin(phi) -
                                      np.cos(theta) * np.sin(inc) * np.cos(phi))

                theta_prime = np.arctan2(np.sin(theta) * np.cos(phi),
                                         np.cos(inc) * (np.cos(theta) * np.cos(phi)) + np.sin(inc) * np.sin(phi))

                # Update the lat and lon with the new rotation
                lat_prime = phi_prime * (180.0 / np.pi)
                lon_prime = theta_prime * (180.0 / np.pi)

                # I am not including any w terms in the vector rotation
                v_x = u * np.sin(theta - obs_theta) + v * np.cos(theta - obs_theta) * np.sin(phi)
                v_y = -u * np.cos(theta - obs_theta) + (v * np.sin(theta - obs_theta) * np.sin(phi))
                v_z = v * np.cos(phi)

                # All this vector math is in Deryl's google doc
                v_x_prime = -v_x * np.cos(inc) - (v_z * np.sin(inc))
                v_y_prime = -v_y
                v_z_prime = v_x * np.sin(inc) + (v_z * np.cos(inc))

                u_prime = (-v_x_prime) * np.sin(theta_prime) + v_y_prime * np.cos(theta_prime)
                v_prime = np.sin(phi_prime) * ((-v_x_prime * np.cos(theta_prime)) - v_y_prime * np.sin(theta_prime)) + v_z_prime * np.cos(phi_prime)

                # Make everything mod 360
                if lat_prime < -90.0:
                    lat_prime = 180.0 + lat_prime
                elif lat_prime > 90.0:
                    lat_prime = lat_prime - 180.0

                if lon_prime < 0.0:
                    lon_prime = 360.0 + lon_prime
                elif lon_prime >= 360.0:
                    lon_prime = lon_prime - 360.0
                 
                df.lat = lat_prime
                df.lon = lon_prime
                df.u = u_prime
                df.v = v_prime
                df['w'] = 0.0
                return df

            df = get_incident_flux(df, sytem_obliquity)
            df = phase_rotation(df, phase)
            df = df.apply(wind_rot, axis=1)

            running_df = pd.DataFrame(columns=['lon', 'lat', 'temp', 'pres', 'alt', 'u', 'v', 'w', 'incident_frac'])
            levels = np.linspace(1, NTAU, NTAU)
            for level in levels:
                sub_df = df[(df['level'] == level)].reset_index(drop=True)

                # Create a dataframe that has all the values
                # Used to find where the 0s should be
                full_df = sub_df.copy()

                # Create a dataframe that has all non zeros
                # Used to find where the 0s should be
                sub_df = sub_df[(sub_df['temp'] > 100)]

                x = np.array(list(sub_df.lon))
                y = np.array(list(sub_df.lat))

                x_test = np.array(list(full_df.lon))
                y_test = np.array(list(full_df.lat))

                t2 = np.round(np.linspace(-89, 89, 50), 3)
                t1 = np.round(np.linspace(0, 359, 50), 3)
                xx, yy = np.meshgrid(t1, t2)

                columns = t1.astype(np.str)
                rows = t2.astype(np.str)

                rbf = Rbf(x, y, np.array(list(sub_df.temp)), epsilon=2, function='linear', smooth=0.1)
                z1 = rbf(xx, yy)
                z1[z1 < 100] = 0

                rbf_test = Rbf(x_test, y_test, np.array(list(full_df.temp)), function='linear', smooth=1)
                z_test = rbf_test(xx, yy)
                df_temp1_test = pd.DataFrame(data=z_test, index=rows, columns=columns)
                temp_test = df_temp1_test.stack().reset_index().rename(columns={"level_0": "lat", "level_1": "lon", 0: "temp"})

                # Set the interpolation functions
                rbf = Rbf(x, y, np.array(list(sub_df.pres)), epsilon=2, function='linear', smooth=0.1)
                z2 = rbf(xx, yy)

                rbf = Rbf(x, y, np.array(list(sub_df.alt)), epsilon=2, function='linear', smooth=0.1)
                z3 = rbf(xx, yy)
                z3[z3 < 1e99] = min(sub_df.alt)

                rbf = Rbf(x, y, np.array(list(sub_df.u)), epsilon=2, function='linear', smooth=0.1)
                z4 = rbf(xx, yy)

                rbf = Rbf(x, y, np.array(list(sub_df.v)), epsilon=2, function='linear', smooth=0.1)
                z5 = rbf(xx, yy)

                rbf = Rbf(x, y, np.array(list(sub_df.incident_frac)), epsilon=2, function='linear', smooth=0.1)
                z6 = rbf(xx, yy)


                # Create more dataframe stuff
                df_temp1 = pd.DataFrame(data=z1, index=rows, columns=columns)
                temp = df_temp1.stack().reset_index().rename(columns={"level_0": "lat", "level_1": "lon", 0: "temp"})

                df_temp2 = pd.DataFrame(data=z2, index=rows, columns=columns)
                pres = df_temp2.stack().reset_index().rename(columns={"level_0": "lat", "level_1": "lon", 0: "pres"})

                df_temp3 = pd.DataFrame(data=z3, index=rows, columns=columns)
                alt = df_temp3.stack().reset_index().rename(columns={"level_0": "lat", "level_1": "lon", 0: "alt"})

                df_temp4 = pd.DataFrame(data=z4, index=rows, columns=columns)
                u = df_temp4.stack().reset_index().rename(columns={"level_0": "lat", "level_1": "lon", 0: "u"})

                df_temp5 = pd.DataFrame(data=z5, index=rows, columns=columns)
                v = df_temp5.stack().reset_index().rename(columns={"level_0": "lat", "level_1": "lon", 0: "v"})

                df_temp6 = pd.DataFrame(data=z6, index=rows, columns=columns)
                incident_frac = df_temp6.stack().reset_index().rename(columns={"level_0": "lat", "level_1": "lon", 0: "incident_frac"})


                # Merge dataframes
                big_df = pd.merge(temp, pres, how='left')
                big_df = pd.merge(big_df, alt, how='left')
                big_df = pd.merge(big_df, u, how='left')
                big_df = pd.merge(big_df, v, how='left')
                big_df = pd.merge(big_df, incident_frac, how='left')

                # More Code
                big_df['level'] = level

                # Filter out stuff where Exotransmit wants 0s
                big_df['temp'][temp_test['temp'] < 100] = 0

                # Merge the dataframes
                frames = [running_df, big_df]
                running_df = pd.concat(frames, sort=True)

            running_df['lat'] = running_df['lat'].astype(float)
            running_df['lon'] = running_df['lon'].astype(float)
            running_df['level'] = running_df['level'].astype(float)
            running_df['w'] = 0

            # Sort the data tables
            running_df = running_df[['lat', 'lon', 'level', 'alt', 'pres', 'temp', 'u', 'v', 'w', 'incident_frac']]
            running_df = running_df.sort_values(by=['lat', 'lon'], axis=0, ascending=[True, True])

            # Put this in the spectra folder
            # That is where it needs to be run from
            df_to_txt('../Spectra/DATA/init_' + planet_name + '_phase_{}_inc_{}.txt'.format(phase, inc), running_df)