#!/usr/bin/env python
from netCDF4 import Dataset
import numpy as np
import netCDF4 as nc
from PIL import *
from pylab import *
import math
import os
import matplotlib.pyplot as plt
import pdb
import argparse
import sys
from scipy.interpolate import interp1d
from iceberg_mass_comparisons import define_paths_array
from iceberg_mass_comparisons import find_max_min_years
from distributions_of_bergs import get_valid_data_from_traj_file
from distributions_of_bergs import add_constraint_to_data
from sea_ice_concentrations import plot_polar_field



def interpolat_bergs_onto_map(lat,lon,area,iceberg_lats,iceberg_lons,field,iceberg_mass,mass_weighted):
	M=area.shape
	field_gridded=np.zeros((M[0],M[1]))
	mass_gridded=np.zeros((M[0],M[1]))

	for k in range(len(field)):
		lat_ind=np.argmin(abs(lat-iceberg_lats[k]))
		lon_ind=np.argmin(abs(lon-iceberg_lons[k]))
		#print 'here', lon_ind,lat_ind

		if mass_weighted==False:
			field_gridded[lat_ind,lon_ind]=field_gridded[lat_ind,lon_ind]+field[k]
		if mass_weighted==True:
			mass_gridded[lat_ind,lon_ind]=mass_gridded[lat_ind,lon_ind]+iceberg_mass[k]
			field_gridded[lat_ind,lon_ind]=field_gridded[lat_ind,lon_ind]+(field[k]*iceberg_mass[k])

	return [field_gridded, mass_gridded]
			


###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
def main():
	#Clear screen
	#os.system('clear')


	#Defining possible paths
	all_paths=define_paths_array()

	#Flags
	plot_average_dist=1
	plot_full_time_series=0

	#Parameters and variables
	#start_year=1980
	end_year0=1985
	number_of_bins=100
	Number_of_years=25

	#Include contraints to use to filter the data
	# The constraint has the form [field_name, lower_bonud, upper_bound]
	constraint1=np.array(['lat',-85.,-10.])  # Latituds north of -60 in the southern hemisphere
	#constraint2=np.array(['lon',-60.,0.])  #Longitude of Weddel Sea    , lon appears to go from -270 to 90. 
	#constraint3=np.array(['lon',-120.,-60.])  #Longitude of Weddel Sea    , lon appears to go from -270 to 90. 
	constraints=[]
	constraints.append(constraint1)
	#constraints.append(constraint2)
	#constraints.append(constraint3)
	#contraints


	mass_scaling=np.array([2000, 200, 50, 20, 10, 5, 2, 1, 1, 1])



	#field_name='area';  x_min=1;x_max=1.e9  # Good for area including all Deltas
	field_name='area';  x_min=1;x_max=1.e7  # Good for area 
	#field_name='mass' ; #x_min=100;x_max=1.e12  # Good for mass
	#Defining a list of colors
	color_vec=np.array(['blue', 'red','purple','green', 'coral', 'cyan', 'magenta','orange', 'black', 'grey', 'yellow', 'orchid', 'blue', 'red','purple','green', 'coral', 'cyan' ])

	fig = plt.figure(1)



	input_file='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta1/19000101.icebergs_month.nc'
	area=nc.Dataset(input_file).variables['area'][:, :]
	lon=nc.Dataset(input_file).variables['xT'][:]
	lat=nc.Dataset(input_file).variables['yT'][:]


	count1=0
	#for k in range(13):#for k in np.array([3]):
	#for k in np.array([0 , 8]):
	for k in np.array([8 ]):
	#for k in np.array([0,1,2,5,8,9]):
		count1=count1+1
	#for k in range(0,8):
		input_folder=all_paths[k]


		#Make sure that data exists of the years allocated, othersize change it.
		(min_year, max_year)=find_max_min_years(input_folder,'.iceberg_trajectories.nc')
		end_year=min(max_year,end_year0)
		start_year=max(end_year-Number_of_years,min_year)
		#if max_year>=start_year:
		for year in range(start_year,end_year):
			#year=1945
			filename ='/' + str(year) + '0101.iceberg_trajectories.nc'
 	                input_file=input_folder + filename
			print input_file
			field_name1='uvel'
			#field_name1='lon'
			field_name2='vvel'
			iceberg_lats = get_valid_data_from_traj_file(input_file,'lat',subtract_orig=False)
			iceberg_lons = get_valid_data_from_traj_file(input_file,'lon',subtract_orig=False)
			iceberg_mass = get_valid_data_from_traj_file(input_file,'mass',subtract_orig=False)
			iceberg_field1 = get_valid_data_from_traj_file(input_file,field_name1,subtract_orig=False)
			iceberg_field2 = get_valid_data_from_traj_file(input_file,field_name2,subtract_orig=False)
			print 'Lendth of original field: ' , len(iceberg_field1)
			#Getting the distributions from the data

			#Handeling constraints
	                iceberg_lats=add_constraint_to_data(input_file,iceberg_lats,constraints)
	                iceberg_lons=add_constraint_to_data(input_file,iceberg_lons,constraints)
	                iceberg_mass=add_constraint_to_data(input_file,iceberg_mass,constraints)
	                iceberg_field1=add_constraint_to_data(input_file,iceberg_field1,constraints)
	                iceberg_field2=add_constraint_to_data(input_file,iceberg_field2,constraints)
	                print 'Lendth of field after constraints: ' , len(iceberg_field1)


			(field_gridded1,mass_gridded)=interpolat_bergs_onto_map(lat,lon,area,iceberg_lats,iceberg_lons,iceberg_field1,iceberg_mass,mass_weighted=True)
			(field_gridded2,mass_gridded)=interpolat_bergs_onto_map(lat,lon,area,iceberg_lats,iceberg_lons,iceberg_field2,iceberg_mass,mass_weighted=True)


			if year==start_year:
				Total_gridded1=field_gridded1
				Total_gridded2=field_gridded2
				Total_mass_gridded= mass_gridded
			else:
				Total_gridded1=Total_gridded1+field_gridded1
				Total_gridded2=Total_gridded2+field_gridded2
				Total_mass_gridded=Total_mass_gridded+ mass_gridded

		M=field_gridded1.shape
		norm_field1=np.zeros((M[0],M[1]))
		norm_field2=np.zeros((M[0],M[1]))
		for i in range(M[0]):
			for j in range(M[1]):
				if mass_gridded[i,j]>0:
					norm_field1[i,j]=field_gridded1[i,j]/mass_gridded[i,j]
					norm_field2[i,j]=field_gridded2[i,j]/mass_gridded[i,j]




		subplot(2,2,1)
		plot_polar_field(lat,lon,field_gridded1,pole='south',difference_on=1.,title='uvel',p_values=None,cscale=None,field=None)
		
		subplot(2,2,2)
		plot_polar_field(lat,lon,field_gridded2,pole='south',difference_on=1.,title='vvel',p_values=None,cscale=None,field=None)
		
		subplot(2,2,3)
		#plot_polar_field(lat,lon,field_gridded1/mass_gridded,pole='south',difference_on=1.,title='uvel',p_values=None,cscale=None,field=None)
		plot_polar_field(lat,lon,norm_field1,pole='south',difference_on=1.,title='uvel',p_values=None,cscale=None,field=None)
		
		subplot(2,2,4)
		#plot_polar_field(lat,lon,field_gridded2/mass_gridded,pole='south',difference_on=1.,title='vvel',p_values=None,cscale=None,field=None)
		plot_polar_field(lat,lon,norm_field2,pole='south',difference_on=1.,title='vvel',p_values=None,cscale=None,field=None)
			
			
			
			#Plotting the distributions
			#ax = fig.add_subplot(1,2,count1)
			#ax = fig.add_subplot(1,1,1)
			#plt.plot(field1,-field2,'o',color=color_vec[k])
			#plt.xlabel(field_name1)
			#plt.ylabel(field_name2)
			#plt.ylim([-75., -40]) #For mass
			#plt.xlim([10.e3, 7.4e11]) #For mass
			#ax.set_yscale('log')


	#plt.legend(loc='upper right', frameon=True)
	#fig = matplotlib.pyplot.gcf()
	#fig.set_size_inches(9,4.5)
	plt.show()


	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

