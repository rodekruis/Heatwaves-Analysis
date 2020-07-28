## Scipt to detect heat wave events, according to 6 different definitions of heat wave, in ERA5 data,
# for one specific location. The six definitions are listed below. For each definition, the percentile 
# is calculated respect to the distribution of the temperatures in the whole 40-year period and the 
# condition needs to hold for 3 consecutive days:

# 1 - The maximum day temperature is higher than the 95th percentile of the maximum day temperatures.
 
# 2 - The maximum day temperature is higher than the 95th percentile of the maximum day temperatures 
# and the minimum night temperature is higher than the 95th percentile of the minimum day temperatures.

# 3 - The average day temperature is higher than the 95th percentile of the average day temperatures.

# 4 - Same as (1), but for the heat index, a combination of temperature and relative humidity.
# (according to https://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml)

# 5 - Same as (2), but for the heat index.

# 6 - Same as (3), but for the heat index.

# Imports
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from datetime import date, datetime, timedelta
from itertools import groupby
from operator import itemgetter
np.set_printoptions(threshold=np.inf)


# Definition of heat index (according to https://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml)

def heat_index(T_celsius, RH):
    
    # Convert Celsius to Farenheit
    T = (9*T_celsius/5)+32
    
    # Simple formula
    HI = 0.5*(T + 61.0 + ((T-68.0)*1.2) + (RH*0.094))
    
    # Corrections to the formula
    for i in range(len(HI)):
        
        hi = HI[i]
        rh = RH[i]
        t = T[i]
    
        if hi >= 80:

            hi = -42.379 + 2.04901523*t + 10.14333127*rh - 0.22475541*t*rh -             0.00683783*t*t - 0.05481717*rh*rh + 0.00122874*t*t*rh +             0.00085282*t*rh*rh - 0.00000199*t*t*rh*rh

            if rh < 13 and t >= 80 and t <= 112:

                hi -= ((13 - rh)/4)*np.sqrt((17-abs(t-95.))/17)

            if rh > 85 and t >= 80 and t <= 87:

                hi += ((rh-85)/10) * ((87-t)/5)
        
        HI[i] = hi
            
    return HI


# Define folders

data_folder_temp = '../data/temp/'
                     
data_folder_dew = '../data/dewpoint/'

figures_folder = '../figures/'


# Define city and its coordinates

city = 'Kampala'
lat_city = 1.37
lon_city = 32.29


# Create file lists for temperature and dewpoint at 2 meters

file_list_temp = []
file_list_dew = []

for file in os.listdir(data_folder_temp):
    # Exclude year 1981 because there seem to be some issues with that dataset
    if '.nc' in file and '1981' not in file:
        file_list_temp.append(file)
        
for file in os.listdir(data_folder_dew):
    # Exclude year 1981 because there seem to be some issues with that dataset
    if '.nc' in file and '1981' not in file:
        file_list_dew.append(file)

# Sort file lists (per year)
file_list_temp.sort()
file_list_dew.sort()

# Check whether the number of files for the two variables is the same
if len(file_list_temp) != len(file_list_dew):
    print('Lenghts of filelists are different!')


# Initialise lists

days = []   # days in format datetime

max_temps = [] # daily maximum temperatures
min_temps = [] # daily minimum temperatures
ave_temps = [] # daily average temperatures

max_HI = [] # daily maximum heat index
min_HI = [] # daily minimum heat index
ave_HI = [] # daily average heat index


for f in range(len(file_list_temp)):
    
    # Add the folder name to the filename
    file_name_temp = data_folder_temp+file_list_temp[f]
    file_name_dew = data_folder_dew+file_list_dew[f]
    
#     print(file_name_temp)
#     print(file_name_dew)

    # Check if the year of the two variables is the same
    if file_name_temp[-7:] != file_name_dew[-7:]:
        print(file_name_temp[-7:]+' != '+file_name_dew[-7:])

    # Open the netCDF file
    f_temp = Dataset(file_name_temp)

    # Define time, temperature, latitudes and longitudes
    time_temp = f_temp.variables['time'][:].data
    temp = f_temp.variables['t2m'][:].data
    lats_temp = f_temp.variables['latitude'][:].data
    lons_temp = f_temp.variables['longitude'][:].data

    # Define the start time of the dataset
    datetime_start_string_complete_temp = f_temp.variables['time'].units

    # Find the indices relative to the chosen latitude and longitude
    ind_lat_temp = np.where(abs(lats_temp-lat_city)==np.min(abs(lats_temp-lat_city)))[0][0]
    ind_lon_temp = np.where(abs(lons_temp-lon_city)==np.min(abs(lons_temp-lon_city)))[0][0]

    # Repeat the same for dewpoint
    f_dew = Dataset(file_name_dew)

    time_dew = f_dew.variables['time'][:].data
    dew = f_dew.variables['d2m'][:].data
    lats_dew = f_dew.variables['latitude'][:].data
    lons_dew = f_dew.variables['longitude'][:].data

    datetime_start_string_complete_dew = f_dew.variables['time'].units

    ind_lat_dew = np.where(abs(lats_dew-lat_city)==np.min(abs(lats_dew-lat_city)))[0][0]
    ind_lon_dew = np.where(abs(lons_dew-lon_city)==np.min(abs(lons_dew-lon_city)))[0][0]
    
    # Restrict the time to the minimum of the times considered by the two datasets
    time_min = np.min([len(time_dew),len(time_temp)])
    
    time_temp = time_temp[:time_min]
    time_dew = time_dew[:time_min]

    # Check whether the times are consistent between the two datasets
    if datetime_start_string_complete_temp != datetime_start_string_complete_dew or time_temp[0] != time_dew[0] or len(time_temp) != len(time_dew) or ind_lat_temp != ind_lat_dew or ind_lon_temp != ind_lon_dew:
        print('Something is inconsistent between the two datasets')
    else:
        # If the times are consinstent, define the time (it does not it matter from which dataset)
        time = time_dew
        # Detect the time respect to which the time is calculated (hours since ...)
        datetime_start_string = datetime_start_string_complete_temp.split(' ')[2:]
        date_start = [int(x) for x in datetime_start_string[0].split('-')]
        time_start = [int(float(x)) for x in datetime_start_string[1].split(':')]

    # Transform from Kelvin to degrees Celsius
    temp_city = temp[:time_min,ind_lat_temp,ind_lon_temp]-273.15
    dew_city = dew[:time_min,ind_lat_dew,ind_lon_dew]-273.15

    # Formula to calculate the relative humidity (RH) from dewpoint
    A = 7.625
    B = 243.04

    RH_city = 100*np.exp(A*B*(dew_city-temp_city)/((B+temp_city)*(B+dew_city)))

    # Calculate heat index for the chosen location
    HI_city = heat_index(temp_city, RH_city)

    # Number of days (since we have hourly data)
    number_days = int(len(time)/24)

#     print('Number days: '+str(number_days))

    # Extract daily minimum, maximum and average temperature and heat index
    for ind_day in range(number_days):

        # For each day, define the indices relative to the start and the end of the day
        ind_start = 24*ind_day

        ind_end = ind_start+24

        # Extract time, temperature and heat index for the chosen day
        hours_delta_day = time[ind_start:ind_end]

        temp_city_day = temp_city[ind_start:ind_end]
        HI_city_day = HI_city[ind_start:ind_end]

        # Time (in hours) relative to the first hour at which the data is calculated (for the chosen day)
        hours_delta = int(hours_delta_day[0])

        # Start date and time in datetime format (this could also be written out of the for-loop)
        start = datetime(date_start[0], date_start[1], date_start[2], hour=time_start[0], minute=time_start[1], second=time_start[2])

        # Difference (in datetime format) between the start and the chosen day
        delta = timedelta(hours=hours_delta)     

        # Time relative to the chosen day (in datetime format)
        day_time = start + delta      

        # From datetime format to day format
        day = day_time.date()

        # Append day and the variables to the lists
        days.append(day)

        max_temps.append(np.max(temp_city_day))
        min_temps.append(np.min(temp_city_day))
        ave_temps.append(np.mean(temp_city_day))

        max_HI.append(np.max(HI_city_day))
        min_HI.append(np.min(HI_city_day))
        ave_HI.append(np.mean(HI_city_day))


# Plot all the variables

# plt.figure()
# plt.plot(max_temps)

# plt.figure()
# plt.plot(min_temps)

# plt.figure()
# plt.plot(ave_temps)

# plt.figure()
# plt.plot(max_HI)

# plt.figure()
# plt.plot(min_HI)

# plt.figure()
# plt.plot(ave_HI)


# Initialise dictionaries

variables = {} # variables
ind_95 = {} # indices at which the value of the variable extends above the 95th percentile 
indexes = {} # indices relative to the occurrence of heat waves (heat waves day)
HW = {} # indices relative to the occurrence of heat waves (only one day per heat wave) 
HW_names = {} # name for the different definitions

# Insert variables calculated above in the dictionary
variables['0'] = max_temps
variables['1'] = min_temps
variables['2'] = ave_temps
variables['3'] = max_HI
variables['4'] = min_HI
variables['5'] = ave_HI

# Define names based on definitions
HW_names['0'] = 'T max'
HW_names['1'] = 'T max and min'
HW_names['2'] = 'T mean'
HW_names['3'] = 'HI max'
HW_names['4'] = 'HI max and min'
HW_names['5'] = 'HI mean'


for k in range(len(variables)):
    
    # Find indices for which the chosen variable is higher than the 95th percentile
    ind_95[str(k)] = np.where(variables[str(k)]>np.percentile(variables[str(k)], 95))[0]
    
    indexes[str(k)] = []
    
    if k == 1:
        
        # Maximum and minimum temperature together
        ind_95_0_1 = []
        
        for i in ind_95['0']:
            if i in ind_95['1']:
                ind_95_0_1.append(i)
        
        ind_95_k = ind_95_0_1
        
    elif k == 4:
        
        # Maximum and minimum heat index together
        ind_95_3_4 = []
        
        for i in ind_95['3']:
            if i in ind_95['4']:
                ind_95_3_4.append(i)
        
        ind_95_k = ind_95_3_4
                
    else:
                
        ind_95_k = ind_95[str(k)]
    
    # Divide the list in sublists of consecutive numbers
    for j, g in groupby(enumerate(ind_95_k), lambda ix : ix[0] - ix[1]):
        indexes[str(k)].append(list(map(itemgetter(1), g)))
        

# All the definitions at the same time
ind_95_tot = []
indexes[str(len(variables))] = []
        
for i in ind_95['0']:
    if all(i in ind_95[k] for k in [str(j) for j in range(1,6)]):
        ind_95_tot.append(i) 
        
for j, g in groupby(enumerate(ind_95_tot), lambda ix : ix[0] - ix[1]):
        indexes[str(len(variables))].append(list(map(itemgetter(1), g)))

# Restrict heat waves days to heat waves events (take only the ones for which the condition lasted for at least 3 days)         
for k in range(len(indexes)):
        
        HW_k = []
        
        for l in indexes[str(k)]:
            
            if len(l) >= 3:
        
                HW_k.append(l[-1])
            
        HW[str(k)] = np.array(HW_k)
        

# Create a bollean vector for each heat wave: for each time in the list 'days'
# - add 0 if there is no heat wave event
# - add 1 if there is a heat wave event

HW_days = {}

for k in range(len(HW)):
    
    # add 1 in the beginning, otherwise the function 'plt.scatter' does not work 
    HW_k = np.array([1])
    
    for d in range(len(days)):
        
        if d in HW[str(k)]:
            
            HW_k = np.append(HW_k, 1)
            
        else:
            
            HW_k = np.append(HW_k, np.nan)
    
    HW_days[str(k)] = HW_k 


## Plot

# Parameters
fontsize = 16
colors = ['r','b','y']*2+['k']
markers = ['o']*3+['*']*3+['d'] # '.', 'o', 'v', 's', 'P', '*', '+', 'x', 'd'
marker_size = 50
name_fig = 'HW_'+city

matplotlib.rc('xtick', labelsize=fontsize) 
matplotlib.rc('ytick', labelsize=fontsize)

# Define figure
fig = plt.figure(figsize=(20,6))

for k in range(len(HW)):
    
    # Add a fake day in the beginning (to take into account the '1' added in the beginning of each HW_days)
    plt.scatter([datetime(1900, 1, 1).date()]+days, (k+1)*HW_days[str(k)], c=colors[k], marker = markers[k], s=marker_size)
    # Number of HW events
    print(len(HW[str(k)]))
    
# Limits for y and x axis
plt.ylim(0,8)
plt.xlim(days[0],days[-1])
    
# Redefine y-labels
str_ticks_list = ['T max','T max and min','T average','HI max','HI max and min','HI average','All of them']
plt.yticks(range(1,8), str_ticks_list);

# Title
plt.title('Heat waves in '+city+' according to different definitions', fontweight = 'bold', fontsize=fontsize*5/4);

# Save figure
fig.savefig(figures_folder+name_fig+'.pdf', format='pdf', bbox_inches = 'tight', pad_inches = 0.3)
fig.savefig(figures_folder+name_fig+'.png', format='png', dpi=300, bbox_inches = 'tight', pad_inches = 0.3)

