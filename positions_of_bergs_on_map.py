#!/usr/bin/env python
from netCDF4 import Dataset
import numpy as np
from PIL import *
from pylab import *
import math
import os
import scipy.io as sc
import matplotlib.pyplot as plt
import pdb
import argparse
import sys
from scipy.interpolate import interp1d
from sea_ice_concentrations import plot_polar_field
from sea_ice_concentrations import generate_data_file
from distributions_of_bergs import get_valid_data_from_traj_file
from mpl_toolkits.basemap import Basemap
from iceberg_mass_comparisons import find_max_min_years
from distributions_of_bergs import add_constraint_to_data


def save_traj_mat_file(lat,lon, Total_berg_lon, Total_berg_lat, mass0,pole,input_folder,start_year,end_year,months_str,constraint_name=None):
	field='traj'
	#Creating filename
	mat_filename='processed_data/'+ str.split(input_folder,'_bergs_')[-1]+'_'+field+'_'+str(start_year)+'_to_'+ str(end_year)+'_'+months_str+'_'+pole
	if constraint_name!=None:
		mat_filename=mat_filename + '_constraint_' + constraint_name
	mat_filename=mat_filename +'.mat' 
	
	sc.savemat(mat_filename, {'lat':lat , 'lon':lon , 'Total_berg_lon':Total_berg_lon , 'Total_berg_lat':Total_berg_lat , 'field':field, 'mass0':mass0 })
	print 'File: ' + mat_filename + ' saved.'



###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
def main():
	#Clear screen
	#os.system('clear')


	start_year=1915
	end_year0=2019
	Number_of_years=0
	input_folder='/ptmp/Alon.Stern/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta9/'
	second_folder=None
	field='CN'
	months_str='all'
	constraint_name=None
	months=range(0,11)
	file_type='ice_month'
	pole='north'
	save_traj_data=False

	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])



 	# The constraint has the form [field_name, lower_bonud, upper_bound]
	#constraint1=np.array(['lat',-65.,-10.])  # Latituds north of -60 in the southern hemisphere

	constraint_depth = {'Constraint_field_name': 'depth', 'lower_bound': 3000 , 'upper_bound': 8000, 'original_values': True}
 
	constraint_dist_from_calving = {'Constraint_field_name': 'distance_from_calving', 'lower_bound': 1000000 , 'upper_bound': 10000000000000, 'original_values': False}

	constraint2 = {'Constraint_field_name': 'lon', 'lower_bound': -150 , 'upper_bound': -65, 'original_values': True} #Longitude of AB
 	constraint3 = {'Constraint_field_name': 'lon', 'lower_bound': -210 , 'upper_bound': -65, 'original_values': True} #Longitude of AB_plus_Ross

 	constraint_SH = {'Constraint_field_name': 'lat', 'lower_bound': -90 , 'upper_bound': 0, 'original_values': True} #age
 	constraint_age = {'Constraint_field_name': 'age', 'lower_bound': 30 , 'upper_bound': 300, 'original_values': True} #age

	constraints=[]

	#mass=3.9e11 ;constraint_name='massD9'
	#mass=8.8e7 ;constraint_name='massD1'
	#constraint_m = {'Constraint_field_name': 'mass', 'lower_bound': mass-100 , 'upper_bound': mass+100, 'original_values': True}
	#constraints.append(constraint2)  ; constraint_name='AB'
	#constraints.append(constraint3)  ; constraint_name='AB_plus_Ross'
	#constraints.append(constraint4)
	#constraints.append(constraint_m)
	#constraints.append(constraint_depth)
	#constraints.append(constraint_dist_from_calving)
	#constraints.append(constraint_age)
	#constraints.append(constraint_SH)

	#plot_multiple_panels==1:
	constraint2 = {'Constraint_field_name': 'lon', 'lower_bound': -150 , 'upper_bound': -65, 'original_values': True} #Longitude of AB

	root_path='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg'

	(All_data, lat, lon,z) =generate_data_file(start_year, start_year, input_folder, field, months_str, file_type)
	data_mean=np.squeeze(np.mean(All_data,axis=0))
	data_mean=None

	#run_names=np.array(['freq','all_big', 'mass', 'all_small']) ; naming_flag='groups'
	#run_names=np.array(['Delta1', 'Delta2', 'Delta3', 'Delta6', 'Delta9','Delta10']) ; naming_flag='Delta'
	#run_names=np.array(['Delta1' ,'Delta10','Delta1_thick', 'Delta10_thick', 'Delta1b_thick', 'Delta10b_thick']) ; naming_flag='Thick'
	#run_names=np.array(['Delta6', 'Delta6','Delta9']) ; naming_flag='Delta'
	run_names=np.array(['tournadre']) ; naming_flag='Delta'
	#run_names=np.array(['Delta1','Delta3','Delta4','Delta5', 'Delta6','Delta9']) ; naming_flag='Delta'

	for k in range(1):
		input_folder=root_path + '_bergs_' + run_names[k]
		subplot(1,1,k)
		print input_folder	
		#Making sure the years work
		(min_year, max_year)=find_max_min_years(input_folder,'.iceberg_trajectories.nc')
		end_year=min(max_year,end_year0)
		start_year=max(end_year-Number_of_years,min_year)
		title1= run_names[k] + ' (' + str(start_year) + ' to ' + str(end_year) + ')'

		#Plot blank map.
		if pole=='south':
			projection = 'splaea'
			boundinglat=-45
		if pole=='north':
			projection = 'nplaea'
			boundinglat=45
		projection = 'nplaea'
		plot_polar_field(lat,lon,data_mean,pole,difference_on=0.,title=title1\
				,p_values=None,cscale=None,field=None,colorbar_on=True,return_data_map=False,plot_lat_lon_lines=True,boundinglat=boundinglat)
		m = Basemap(projection=projection, boundinglat=boundinglat, lon_0=180)
		print 'You are here'
		print projection

		
		count=0
		Total_berg_lat=np.array([])
		Total_berg_lon=np.array([])
		Total_berg_mass0=np.array([])
		for year in range(start_year, end_year+1):
			count=count+1
			filename ='/' + str(year) + '0101.iceberg_trajectories.nc'
			input_file=input_folder + filename
			print input_file
			all_bergs_lat = get_valid_data_from_traj_file(input_file,'lat')
			all_bergs_lon = get_valid_data_from_traj_file(input_file,'lon')
			all_bergs_mass0 =  get_valid_data_from_traj_file(input_file,'mass',subtract_orig=False,get_original_values=True) # Getting iceberg mass too.

			#Adding constraints to the data:
			print 'Lendth of original field: ' , len(all_bergs_lat)
			all_bergs_lat=add_constraint_to_data(input_file,all_bergs_lat,constraints)
			all_bergs_lon=add_constraint_to_data(input_file,all_bergs_lon,constraints)
			all_bergs_mass0=add_constraint_to_data(input_file,all_bergs_mass0,constraints)

			print 'Lendth of field after constraints: ' , len(all_bergs_lat)


			x, y = m(all_bergs_lon, all_bergs_lat)
			plot_polar_field(lat,lon,data_mean,pole,difference_on=0.,title=title1\
					,p_values=None,cscale=None,field=None,colorbar_on=True,return_data_map=False,plot_lat_lon_lines=True,boundinglat=boundinglat)
			m.scatter(x,y,3,marker='o',color=color_vec[k])

			#plt.plot(all_bergs_lon,all_bergs_lat,'o')
			
			Total_berg_lon=np.concatenate((Total_berg_lon,x),axis=0)
			Total_berg_lat=np.concatenate((Total_berg_lat,y),axis=0)
			Total_berg_mass0=np.concatenate((Total_berg_mass0,all_bergs_mass0),axis=0)
			
		
		if save_traj_data==True:
			save_traj_mat_file(lat,lon, Total_berg_lon, Total_berg_lat, Total_berg_mass0, pole,input_folder,start_year,end_year,months_str,constraint_name)
		#plot_polar_field(lat,lon,data_mean,pole,difference_on=0.,title=title1)
		#m.scatter(Total_berg_lon,Total_berg_lat,1,marker='o',color=color_vec[k])
			


	output_file= 'figures/traj_plot_' + str(start_year)+ '_to_' +  str(end_year) + '_with_' + run_names[k] + '.png'
	#plt.savefig(output_file, dpi=150, bbox_inches='tight', pad_inches=0.4)

	plt.show()



	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

