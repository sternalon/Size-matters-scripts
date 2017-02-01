#!/usr/bin/env python
from netCDF4 import Dataset
import numpy as np
from PIL import *
from pylab import *
import math
import os
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



def construct_linear_sum(Big_data_matrix,weights):
	Total_weight=0.
	M=Big_data_matrix.shape
	Total_sum=np.zeros((M[1],M[2]))

	for k in range(M[0]):
		Total_sum=Total_sum+(np.squeeze(Big_data_matrix[k,:,:])*weights[k])
		Total_weight=Total_weight+weights[k]
	Total_sum=Total_sum/Total_weight
	return Total_sum







###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
def main():
	#Clear screen
	#os.system('clear')


	start_year=1980
	end_year0=1980
	field='melt'
	months_str='all'
	months=range(0,11)
	file_type='icebergs'
	pole='south'
	Number_of_years=20

	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	second_folder=None

	distribution=np.array([0.24, 0.12, 0.15, 0.18, 0.12, 0.07, 0.03, 0.03, 0.03, 0.02])
	initial_mass=np.array([8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11])   
	Total_mass=np.sum(distribution*initial_mass)

	#plot_multiple_panels==1:

	root_path='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg'
	#run_names=np.array(['freq','all_big', 'mass', 'all_small']) ; naming_flag='groups'
	#run_names=np.array(['Delta1', 'Delta2', 'Delta3', 'Delta6', 'Delta9','Delta10']) ; naming_flag='Delta'
	run_names=np.array(['Delta1', 'Delta2', 'Delta3','Delta4','Delta5','Delta6','Delta7', 'Delta8', 'Delta9','Delta10']) ; naming_flag='Delta'
	#run_names=np.array(['Delta1', 'Delta9']) ; naming_flag='Delta'
	Number_of_files=10

	multi_berg_run_name='freq'

	if multi_berg_run_name=='freq':
		weights=np.zeros((10,1))
		for i in range(10):
			weights[i]=distribution[i]*initial_mass[i]/Total_mass
		print sum(weights)

	count=-1
	for k in range(Number_of_files):
		count=count+1
		input_folder=root_path + '_bergs_' + run_names[k]
		(min_year, max_year)=find_max_min_years(input_folder,'.' + file_type + '_month.nc')
		end_year=max_year
		start_year=max(end_year-Number_of_years+1,min_year)

		(All_data, lat, lon) =generate_data_file(start_year, end_year, input_folder, field,months,file_type)
		data_mean=np.squeeze(np.mean(All_data,axis=0))

		#Set up the matrix on the first time around
		if count==0:
			M=data_mean.shape
			Big_data_matrix=np.zeros((Number_of_files,M[0],M[1]))
		
		Big_data_matrix[count,:,:]=data_mean

	linear_sum_data=construct_linear_sum(Big_data_matrix,weights)

	

	#Getting the data for the multi_berg_run
	input_folder=root_path + '_bergs_' + multi_berg_run_name
	(min_year, max_year)=find_max_min_years(input_folder,'.' + file_type + '_month.nc')
	end_year=max_year
	start_year=max(end_year-Number_of_years+1,min_year)
	(All_data, lat, lon) =generate_data_file(start_year, end_year, input_folder, field,months,file_type)
	data_mean=np.squeeze(np.mean(All_data,axis=0))

	relative_error=(data_mean-linear_sum_data)#/data_mean
	#ralative_error=relative_error*np.where(data_mean!=0)

	#Plotting the linear sum	
	subplot(1,3,1)
	plot_polar_field(lat,lon,linear_sum_data,pole,difference_on=0.,title='linear combination',p_values=None,cscale=None,field=field)


	subplot(1,3,2)
	plot_polar_field(lat,lon,data_mean,pole,difference_on=0.,title=multi_berg_run_name,p_values=None,cscale=None,field=field)


	subplot(1,3,3)
	plot_polar_field(lat,lon,-relative_error,pole,difference_on=0.,title='relative_error',p_values=None,cscale=None,field=field)

	output_file='linearity_test_neg' + multi_berg_run_name + '.png'
	plt.savefig(output_file, dpi=150, bbox_inches='tight', pad_inches=0.4)





	plt.show()



	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

