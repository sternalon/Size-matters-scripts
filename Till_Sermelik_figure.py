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


	start_year=2010
	end_year0=2010
	Number_of_years=3
	input_folder='/ptmp/Alon.Stern/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta9/'
	second_folder=None
	field='CN'
	months_str='all'
	constraint_name=None
	months=range(0,11)
	file_type='ice_month'
	pole='north'
	#pole='south'
	pole=None
	save_traj_data=False

	Title_list={'Delta1':'$L_{0}=60m$','Delta2':'$L_{0}=100m$','Delta3':'$L_{0}=200m$','Delta4':'$L_{0}=350m$','Delta5':'$L_{0}=500m$',\
		                        'Delta6':'$L_{0}=700m$','Delta7':'$L_{0}=900m$','Delta8':'$L_{0}=1200m$','Delta9':'$L_{0}=1600m$','Delta10':'$L_{0}=2200m$'}
	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	#cNorm = mpl.colors.LogNorm(vmin=8.8e7, vmax=7.4e11)
	#cNorm = mpl.colors.LogNorm(vmin=1.0e5, vmax=1.0e7)
	cNorm = mpl.colors.Normalize(vmin=0, vmax=0.1)
	#cNorm = mpl.colors.LogNorm(vmin=1.0e-3, vmax=1.0e0)
	
	plot_each_time=False


 	# The constraint has the form [field_name, lower_bonud, upper_bound]
	#constraint1=np.array(['lat',-65.,-10.])  # Latituds north of -60 in the southern hemisphere

	constraint_depth = {'Constraint_field_name': 'depth', 'lower_bound': 3000 , 'upper_bound': 8000, 'original_values': True}
	constraint_age = {'Constraint_field_name': 'age', 'lower_bound': 10 , 'upper_bound': 100, 'original_values': False}
 
	constraint_dist_from_calving = {'Constraint_field_name': 'distance_from_calving', 'lower_bound': 1000000 , 'upper_bound': 10000000000000, 'original_values': False}

	constraint_lon_Sermalik = {'Constraint_field_name': 'lon', 'lower_bound': -38 , 'upper_bound': -36, 'original_values': True} #Longitude of Sermelik
	constraint_lat_Sermalik = {'Constraint_field_name': 'lat', 'lower_bound': 65.0 , 'upper_bound': 66.0, 'original_values': True} #Longitude of Sermalik
	#constraint_mass0 = {'Constraint_field_name': 'mass', 'lower_bound': 3.2*(10**9) , 'upper_bound': 3.4*(10**9), 'original_values': True, 'subtract_original': False}
	#constraint_mass1 = {'Constraint_field_name': 'mass', 'lower_bound': 0*(10**9) , 'upper_bound': 1.0*(10**7), 'original_values': False, 'subtract_original': False}
	constraint_mass0 = {'Constraint_field_name': 'mass', 'lower_bound': 7.3*(10**11) , 'upper_bound': 7.6*(10**11), 'original_values': True, 'subtract_original': False}
	constraint_mass1 = {'Constraint_field_name': 'mass', 'lower_bound': 0*(10**9) , 'upper_bound': 1.0*(10**9), 'original_values': False, 'subtract_original': False}
	constraint_mass = {'Constraint_field_name': 'mass', 'lower_bound': 1.0*(10**10) , 'upper_bound': 4.0*(10**15), 'original_values': False, 'subtract_original': False}
	constraint_lon_WestGreenland = {'Constraint_field_name': 'lon', 'lower_bound': -45 , 'upper_bound': -20, 'original_values': True} 
	constraint_lat_SouthGreenland = {'Constraint_field_name': 'lat', 'lower_bound': 50.0 , 'upper_bound':66.0, 'original_values': True} 

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
	#constraints.append(constraint_lon_Sermalik)
	#constraints.append(constraint_lat_Sermalik)
	#constraints.append(constraint_SH)
	constraints.append(constraint_mass)
	constraints.append(constraint_mass0)
	#constraints.append(constraint_age)
	#constraints.append(constraint_lon_WestGreenland)
	constraints.append(constraint_lat_SouthGreenland)

	#plot_multiple_panels==1:
	constraint2 = {'Constraint_field_name': 'lon', 'lower_bound': -150 , 'upper_bound': -65, 'original_values': True} #Longitude of AB

	root_path='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg'

	(All_data, lat, lon,z) =generate_data_file(start_year, start_year, input_folder, field, months_str, file_type)
	data_mean=np.squeeze(np.mean(All_data,axis=0))
	data_mean=None

	#Deciding on the projection.
	if pole=='south':
		projection = 'splaea'
		boundinglat=-45
	if pole=='north':
		projection = 'nplaea'
		boundinglat=45
	if pole==None:
		projection = 'lcc'
		lat_1=80
		lat_2=60
		lat_0=63.5
		lon_0=-59.


	run_names=np.array(['rolling_tournadre_Burton_new','rolling_tournadre_Rolling_off'])
	Name_list=np.array(['Rolling', 'No Rolling'])

	nrows=1 ;ncols=2
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	for k in range(ncols*nrows):
		input_folder=root_path + '_bergs_' + run_names[k] + '/trajectories/'

		subplot(nrows,ncols,k+1)
		print input_folder	
		#Making sure the years work
		(min_year, max_year)=find_max_min_years(input_folder,'.iceberg_trajectories.nc')
		end_year=min(max_year,end_year0)
		start_year=max(end_year-Number_of_years,min_year)
		title1= run_names[k] + ' (' + str(start_year) + ' to ' + str(end_year) + ')'

		#plot_polar_field(lat,lon,data_mean,pole,difference_on=0.,title=title1\
		#		,p_values=None,cscale=None,field=None,colorbar_on=True,return_data_map=False,plot_lat_lon_lines=True,boundinglat=boundinglat)
		if pole==None:
			#m = Basemap(width=1750000,height=4300000, rsphere=(6378137.00,6356752.3142),resolution='l',area_thresh=1000.,\
			#		projection=projection,lat_1=lat_1,lat_2=lat_2,lat_0=lat_0,lon_0=lon_0)
			projection = 'lcc'
	                lat_1=80
        	        lat_2=60
                	lat_0=64.5
	                lon_0=-50.
        	        m = Basemap(width=2700000,height=4500000, rsphere=(6378137.00,6356752.3142),resolution='l',area_thresh=1000.,\
                                                                        projection=projection,lat_1=lat_1,lat_2=lat_2,lat_0=lat_0,lon_0=lon_0)

		else:
			m = Basemap(projection=projection, boundinglat=boundinglat, lon_0=180)
		m.drawcoastlines()
		#m.fillcontinents(color='grey',lake_color='white')

		
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
			#all_bergs_mass0 =  get_valid_data_from_traj_file(input_file,'mass',subtract_orig=False,get_original_values=False) # Getting iceberg mass too.
			#all_bergs_mass0 =  get_valid_data_from_traj_file(input_file,'age',subtract_orig=False,get_original_values=False) # Getting iceberg mass too.
			#all_bergs_mass0 =  get_valid_data_from_traj_file(input_file,'Me',subtract_orig=False,get_original_values=False) # Getting iceberg mass too.
			all_bergs_mass0 =  get_valid_data_from_traj_file(input_file,'cn',subtract_orig=False,get_original_values=False) # Getting iceberg mass too.

			#Adding constraints to the data:
			print 'Lendth of original field: ' , len(all_bergs_lat)
			all_bergs_lat=add_constraint_to_data(input_file,all_bergs_lat,constraints)
			all_bergs_lon=add_constraint_to_data(input_file,all_bergs_lon,constraints)
			all_bergs_mass0=add_constraint_to_data(input_file,all_bergs_mass0,constraints)

			print 'Lendth of field after constraints: ' , len(all_bergs_lat)


			x, y = m(all_bergs_lon, all_bergs_lat)
			if plot_each_time==True:
				plot_polar_field(lat,lon,data_mean,pole,difference_on=0.,title=title1\
						,p_values=None,cscale=None,field=None,colorbar_on=True,return_data_map=False,plot_lat_lon_lines=True,boundinglat=boundinglat)
				m.scatter(x,y,3,marker='.',color=color_vec[k])

			#plt.plot(all_bergs_lon,all_bergs_lat,'o')
			
			Total_berg_lon=np.concatenate((Total_berg_lon,x),axis=0)
			Total_berg_lat=np.concatenate((Total_berg_lat,y),axis=0)
			Total_berg_mass0=np.concatenate((Total_berg_mass0,all_bergs_mass0),axis=0)
		

		datamap=m.scatter(Total_berg_lon,Total_berg_lat, c=Total_berg_mass0, marker='.',cmap='jet',norm=cNorm, lw = 0.0)
		if run_names[k]=='tournadre':
			cbar=fig.colorbar(datamap)
        	        cbar.set_label('Calving Mass (kg)')
		else:
			plt.title(Name_list[k])
		
		if save_traj_data==True:
			save_traj_mat_file(lat,lon, Total_berg_lon, Total_berg_lat, Total_berg_mass0, pole,input_folder,start_year,end_year,months_str,constraint_name)
		#plot_polar_field(lat,lon,data_mean,pole,difference_on=0.,title=title1)
		#m.scatter(Total_berg_lon,Total_berg_lat,1,marker='o',color=color_vec[k])
			

	output_file= 'figures/Leigh_Rink_figure2.png'
	#output_file= 'figures/Northern_Trajectories_tournadre.png'
	#output_file= 'figures/Northern_Trajectories_Deltas.png'
	fig.set_size_inches(21.0, 8.5,forward=True)
	plt.savefig(output_file, dpi=150, bbox_inches='tight', pad_inches=0.4)

	plt.show()



	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

