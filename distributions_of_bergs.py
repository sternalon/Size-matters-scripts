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
import scipy.io as sc
from iceberg_mass_comparisons import find_max_min_years

def get_valid_data_from_traj_file(filename,field_name,subtract_orig=False,get_original_values=False):

	#Importing the data file
	nc = Dataset(filename, mode='r')

	#Use the length field to help filter the data
	length = nc.variables['length'][:]
	js=np.where(length>0)[0]

	#Special case for field_name==area
	getting_area=0.
	getting_LW_max=0.
	if field_name=='area':
		field_name='width' ; getting_area=1.
	if field_name=='LW_max':
		field_name='width' ; getting_LW_max=1.

	#Special case for iceberg age
	if field_name=='age':
		nc = Dataset(filename, mode='r')
		
		field = nc.variables['day'][:]
		field1=field[js]
		field2=field[js-1]  #This only works when just one position is recorded at a time
		field_day=field1-field2
		
		field = nc.variables['year'][:]
		field1=field[js]
		field2=field[js-1]  #This only works when just one position is recorded at a time
		field_year=field1-field2
		#field_year[np.where(field_year>1)[0]]=(field_year[np.where(field_year>1)[0]])
		#special fix for bug in the model:  (earlier simulation before bug fix)
		field_year[np.where(field_year>1)[0]]=(field_year[np.where(field_year>1)[0]]/2)

		field=(field_year)+(field_day/365.25)
	
	elif field_name=='distance_from_calving':
		subtract_orig=True
		R_earth=6360.*1000
		lat_field = get_valid_data_from_traj_file(filename,'lat',subtract_orig='False',get_original_values='False')
		lat_dist = get_valid_data_from_traj_file(filename,'lat',subtract_orig,get_original_values)
		lon_dist = get_valid_data_from_traj_file(filename,'lon',subtract_orig,get_original_values)
		
		#Converting to distance
		y_dist=lat_dist*R_earth*(np.pi/180)
		x_dist=lon_dist*R_earth*(np.pi/180)*cos(lat_field*np.pi/180)
		field=sqrt(  ((x_dist)**2)  +  ((y_dist)**2)   )
	
	elif field_name=='speed':
		uvel = get_valid_data_from_traj_file(filename,'uvel',subtract_orig,get_original_values)
		vvel = get_valid_data_from_traj_file(filename,'vvel',subtract_orig,get_original_values)
		field=sqrt(uvel**2 +vvel**2)

	elif field_name=='speed_o':
		uo = get_valid_data_from_traj_file(filename,'uo',subtract_orig,get_original_values)
		vo = get_valid_data_from_traj_file(filename,'vo',subtract_orig,get_original_values)
		field=sqrt(uo**2 +vo**2)

	elif field_name=='speed_a':
		ua = get_valid_data_from_traj_file(filename,'ua',subtract_orig,get_original_values)
		va = get_valid_data_from_traj_file(filename,'va',subtract_orig,get_original_values)
		field=sqrt(ua**2 +va**2)

	elif field_name=='Mb': #melt due to basal melt
		uo = get_valid_data_from_traj_file(filename,'uo',subtract_orig,get_original_values)
		vo = get_valid_data_from_traj_file(filename,'vo',subtract_orig,get_original_values)
		uvel = get_valid_data_from_traj_file(filename,'uvel',subtract_orig,get_original_values)
		vvel = get_valid_data_from_traj_file(filename,'vvel',subtract_orig,get_original_values)
		dvo = np.sqrt( ((uvel-uo)**2 + (vvel-vo)**2))
		L = get_valid_data_from_traj_file(filename,'length',subtract_orig,get_original_values)
		SST = get_valid_data_from_traj_file(filename,'sst',subtract_orig,get_original_values)
		Mb=0.58*(dvo**0.8)*(SST+4.0)/(L**0.2)
		Mb[np.where(Mb<0.0)]=0.0
		field=Mb
	
	elif field_name=='side_decay':
		Me = get_valid_data_from_traj_file(filename,'Me',subtract_orig,get_original_values)
		L = get_valid_data_from_traj_file(filename,'length',subtract_orig,get_original_values)
		W = get_valid_data_from_traj_file(filename,'width',subtract_orig,get_original_values)
		H = get_valid_data_from_traj_file(filename,'thickness',subtract_orig,get_original_values)
		field= H*(L+W)*Me


	elif field_name=='Ss':  #Sea state
		uo = get_valid_data_from_traj_file(filename,'uo',subtract_orig,get_original_values)
		vo = get_valid_data_from_traj_file(filename,'vo',subtract_orig,get_original_values)
		ua = get_valid_data_from_traj_file(filename,'ua',subtract_orig,get_original_values)
		va = get_valid_data_from_traj_file(filename,'va',subtract_orig,get_original_values)
		dva = np.sqrt( ((ua-uo)**2 + (va-vo)**2))
		field=1.5*(dva**0.5)+(0.1*dva)

	elif field_name=='Me':  #melt due to side erosion
		Ss = get_valid_data_from_traj_file(filename,'Ss',subtract_orig,get_original_values)
		IC = get_valid_data_from_traj_file(filename,'cn',subtract_orig,get_original_values)
		SST = get_valid_data_from_traj_file(filename,'sst',subtract_orig,get_original_values)
		Me= (1./12.)*(SST+2.)*Ss*(1+cos(np.pi*(IC**3))) 
		Me[np.where(Me<0.0)]=0.0
		field=Me

	elif field_name=='Mv':  #melt due to buoyant convection
		SST = get_valid_data_from_traj_file(filename,'sst',subtract_orig,get_original_values)
		Mv= (7.62e-3*SST)+(1.29e-3*(SST**2))
		Mv[np.where(Mv<0.0)]=0.0
		field=Mv

	else:
		#Importing the field
		field = nc.variables[field_name][:]

		if subtract_orig==True:
			field1=field[js]
			field2=field[js-1]  #This only works when just one position is recorded at a time

			field=field1-field2
			if field_name=='lon':
				field[np.where(field>90)]=field[np.where(field>90)]-360

		elif get_original_values==True:  
			field=field[js-1]  #This only works when just one position is recorded at a time

		else:
			field=field[js]

		#Special case for Area again
		if getting_area==1.: 
			field=field*length[js]
		if getting_LW_max==1.: 
			temp=length[js]
			for k in range(len(field)):
				field[k]=max(field[k],temp[k])
	
	return field


def smooth(y, box_pts):
	box = np.ones(box_pts)/box_pts
	y_smooth = np.convolve(y, box, mode='same')
	y_smooth[0:box_pts]=0.
	y_smooth[len(y)-box_pts:len(y)]=0.
	y_smooth[-1]=0.
	return y_smooth

def remove_neg(f):
	f=0.5*(abs(f) + f)
	return f

def create_size_dist(field,number_of_bins,x1_min=1,x1_max=1.e7):
	#Input the Cumulative distribution
	x=sort(field)

	#Create a histogram on bins defined by x1
	#x1 = np.linspace(x1_min, x1_max, num=number_of_bins, endpoint=True)
	x1 = np.logspace(np.log10(x1_min), np.log10(x1_max), num=number_of_bins, endpoint=True) 

	hist1, x1 =np.histogram(x, bins=x1)

	#Form the probabilty distribution
	x2=x1[:-1] + diff(x1)/2
	prob_dist=hist1/(float(len(x))*diff(x1))  #To get a probabilty dist

	#print 'prob sum is ',np.sum(prob_dist*diff(x1))

	Cu_dist1=np.cumsum(prob_dist*diff(x1))

	return [prob_dist, Cu_dist1, x2]


def add_constraint_to_data(input_file,field,constraints):
	valid_data_boolean=np.ones((len(field)))
	if constraints!=None:
		for k in range(len(constraints)):
			current_constraint=constraints[k]
			constraint_name=current_constraint['Constraint_field_name']
			lower_bound=current_constraint['lower_bound']
			upper_bound=current_constraint['upper_bound']
			get_original_values=current_constraint['original_values']

			if 'subtract_original' in current_constraint.keys():
				subtract_orig=current_constraint['subtract_original']
			else:
				subtract_orig=False


			if constraint_name=='depth':
				#lat_field = get_valid_data_from_traj_file(input_file,'lat',subtract_orig,get_original_values)
				#lon_field = get_valid_data_from_traj_file(input_file,'lon',subtract_orig,get_original_values)
				
				lat_field = get_valid_data_from_traj_file(input_file,'lat')
				lon_field = get_valid_data_from_traj_file(input_file,'lon')
				constraint_field=lat_field*0.
				#print 'here', field[0],lat_field[0]

				#Load topography from a file
				latlon_filename='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta9/19000101.ice_month.nc'
				nc = Dataset(latlon_filename, mode='r')
				lat = nc.variables['yT'][:] 
				lon = nc.variables['xT'][:]
				
				topog_filename='/home/aas/Iceberg_Project/iceberg_scripts/python_scripts/data/topog.nc'
				nc = Dataset(topog_filename, mode='r')
				depth = nc.variables['depth'][:,:] #ny by nx

				#Using the ocean file takes really long.
				#ocean_filename='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta9/19000101.ocean_month_z.nc'
				#nc = Dataset(ocean_filename, mode='r')
				#lat = nc.variables['yh'][:] 
				#lon = nc.variables['xq'][:]
				#depth = nc.variables['depth_ocean'][:,:] #ny by nx



				#Find depth of closest lat/lon point
				for k in range(len(lat_field)):
					#current_lat=lat_field[k]
					lat_ind=(abs(lat-lat_field[k])).argmin()
					lon_ind=(abs(lon-lon_field[k])).argmin()

					#print k, lon_field[k]
					#if lon_ind==388:
					#	print k, lon_field[k]
					#	print lon
					#	print len(lon)
					#print (abs(lat-lat_field[k]))

					#if lon_ind>=360:
					#	lon_ind=lon_ind-360
					
					#print 'lat', lat_ind, lat[lat_ind]
					#print 'lon', lon_ind
					#print lon[lon_ind]

					#index_min=values.argmin()
					#lat_ind=lat.index(test)
					#lat_ind=lat.index(min(abs(lat-current_lat)))
					#lon_ind=lon.index(min(abs(lon-current_lon)))
					#print constraint_field[k]
					constraint_field[k]=depth[lat_ind,lon_ind]
					#print constraint_field[k]

			elif constraint_name=='distance_from_calving':
				subtract_orig=True
				R_earth=6360.*1000
				lat_field = get_valid_data_from_traj_file(input_file,'lat',subtract_orig='False',get_original_values='False')
				lat_dist = get_valid_data_from_traj_file(input_file,'lat',subtract_orig,get_original_values)
				lon_dist = get_valid_data_from_traj_file(input_file,'lon',subtract_orig,get_original_values)
				
				#Converting to distance
				y_dist=lat_dist*R_earth*(np.pi/180)
				x_dist=lon_dist*R_earth*(np.pi/180)*cos(lat_field*np.pi/180)
				constraint_field=sqrt(  ((x_dist)**2)  +  ((y_dist)**2)   )


			else:
				constraint_field = get_valid_data_from_traj_file(input_file,constraint_name,subtract_orig,get_original_values)
			
			
			#if 'absolute_values' in  current_constraint.values():
			#	if current_constraint['absolute_values']==Ture:
			#		constraint_field=abs(constraint_field)


			#print constraint_name,upper_bound, lower_bound
			#print max(constraint_field), min(constraint_field)

			valid_data_boolean=valid_data_boolean*((constraint_field>lower_bound) * (constraint_field<upper_bound))



	js=np.where(valid_data_boolean)[0]
	#print constraint_name, max(constraint_field[js]), min(constraint_field[js])
	
	#Use the constraint field to help constrain the data
	field=field[js]
	
	return field

def get_number_of_icebergs(input_folder, field_name,start_year, end_year, \
		number_of_bins,x_min,x_max,constraints,split_distribution_by_mass,save_mat_file):
	print 'Generating ' + field_name + ' data from folder: ' + input_folder

	number_matrix_all=np.zeros((end_year-start_year+1,1))
	number_matrix_constraints=np.zeros((end_year-start_year+1,1))
	number_matrix_splits=np.zeros((end_year-start_year+1,10))
	if split_distribution_by_mass==True:
		initial_mass=np.array([8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11]) 
		mass_scaling=np.array([2000, 200, 50, 20, 10, 5, 2, 1, 1, 1])


	constraints_orig=constraints
	count=-1	
	for year in range(start_year, end_year+1):
		count=count+1
		filename ='/' + str(year) + '0101.iceberg_trajectories.nc'
		input_file=input_folder + filename
		print input_file
		field = get_valid_data_from_traj_file(input_file,field_name)
		print 'Lendth of original field: ' , len(field)
		number_matrix_all[count]=len(field)

		if split_distribution_by_mass==True:
			#for mass in range(len(initial_mass)):
			mass_count=-1
			for k in range(len(initial_mass)):
				if k==0:
					constraint_m = {'Constraint_field_name': 'mass', 'lower_bound': 0 , 'upper_bound': initial_mass[k], 'original_values': False}
				else:
					constraint_m = {'Constraint_field_name': 'mass', 'lower_bound': initial_mass[k-1] , 'upper_bound': initial_mass[k], 'original_values': False}
				mass_constraints=[]
				#mass_constraints=constraints
				mass_constraints.append(constraint_m)
				field_m=add_constraint_to_data(input_file,field,mass_constraints)
				print 'Lendth of field in mass class ', initial_mass[k], ' is : ' , len(field_m)
				number_matrix_splits[count,k]=len(field_m)
			number_matrix_constraints[count,0]=sum(number_matrix_splits[count,:])
		else:
			#Handeling constraints
			field=add_constraint_to_data(input_file,field,constraints)
			print 'Lendth of field after constraints: ' , len(field)
			number_matrix_constraints[count,0]=len(field)

	
	return [number_matrix_all, number_matrix_constraints,number_matrix_splits]


def  generate_dist_timeseries_matrix(input_folder,field_name, start_year, end_year, number_of_bins,x_min,x_max,constraints=None,split_distribution_by_mass=False,save_mat_file=False):
	print 'Generating ' + field_name + ' data from folder: ' + input_folder

	if split_distribution_by_mass==True:
		initial_mass=np.array([8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11]) 
		mass_scaling=np.array([2000, 200, 50, 20, 10, 5, 2, 1, 1, 1])
		#initial_mass=np.array([8.8e7,  3.9e11, 7.4e11]) 
		#initial_mass=np.array([ 3.9e11]) 
		#mass_scaling=np.array([2000, 200, 50, 20, 10, 5, 2, 1, 1, 1])
		split_prob_matrix=np.zeros((end_year-start_year+1,number_of_bins-1,len(initial_mass)))
		split_Cu_matrix=np.zeros((end_year-start_year+1,number_of_bins-1,len(initial_mass)))
	else:
		split_prob_matrix=None
		split_Cu_matrix=None

		


	prob_dist_matrix=np.zeros((end_year-start_year+1,number_of_bins-1))
	Cu_dist_matrix=np.zeros((end_year-start_year+1,number_of_bins-1))
	x_matrix=np.zeros((end_year-start_year+1,number_of_bins-1))


	constraints_orig=constraints
	count=-1	
	for year in range(start_year, end_year+1):
		count=count+1
		filename ='/' + str(year) + '0101.iceberg_trajectories.nc'
		input_file=input_folder + filename
		print input_file
		field = get_valid_data_from_traj_file(input_file,field_name)
		print 'Lendth of original field: ' , len(field)
		
		

		if split_distribution_by_mass==True:
			Number_list=[]
			prob_list=[]
			Cu_list=[]

			

			#for mass in range(len(initial_mass)):
			for mass in initial_mass:

				constraint_m = {'Constraint_field_name': 'mass', 'lower_bound': mass-100 , 'upper_bound': mass+100, 'original_values': True}
				mass_constraints=[]
				mass_constraints=constraints
				mass_constraints.append(constraint_m)

				#Handeling constraints
				field_m=add_constraint_to_data(input_file,field,mass_constraints)
				print 'Lendth of field in mass class ', mass, ' is : ' , len(field_m)


				(prob_dist_m, Cu_dist_m, x)=create_size_dist(field_m,number_of_bins,x_min,x_max)
				if len(field_m)== 0:
					prob_dist_m[:]=0
					Cu_dist_m[:]=0

				print 'A)prob sum is ',np.sum(prob_dist_m[:-1]*diff(x))
				print 'A) Cu max is',max(Cu_dist_m)

				prob_list.append(prob_dist_m)
				Cu_list.append(Cu_dist_m)
				Number_list.append(len(field_m))
				mass_constraints.pop()

			
			#Normalizing
			scaling_total=0
			prob_dist=x*0
			Cu_dist=x*0
			scaling_on=True
			for k in range(len(prob_list)):
				if scaling_on==False:
					prob_dist=prob_dist+(prob_list[k]*Number_list[k])
					Cu_dist=Cu_dist+(Cu_list[k]*Number_list[k])
					scaling_total=scaling_total+Number_list[k]
					split_prob_matrix[count,:,k]=prob_list[k]
					split_Cu_matrix[count,:,k]=Cu_list[k]
				else:
					prob_dist=prob_dist+(prob_list[k]*Number_list[k]*mass_scaling[k])
					Cu_dist=Cu_dist+(Cu_list[k]*Number_list[k]*mass_scaling[k])
					scaling_total=scaling_total+(Number_list[k]*mass_scaling[k])
					split_prob_matrix[count,:,k]=prob_list[k]
					split_Cu_matrix[count,:,k]=Cu_list[k]
			
			prob_dist=prob_dist/scaling_total
			Cu_dist=Cu_dist/scaling_total

			
			#dx1=diff(x)
			print 'prob sum is ',np.sum(prob_dist[:-1]*diff(x))
			print 'Cu max is',max(Cu_dist)




		else:
			#Handeling constraints
			field=add_constraint_to_data(input_file,field,constraints)
			print 'Lendth of field after constraints: ' , len(field)

			#(prob_dist, Cu_dist, x)=create_size_dist(field,number_of_bins,x_max,x_min)
			(prob_dist, Cu_dist, x)=create_size_dist(field,number_of_bins,x_min,x_max)
			#(prob_dist,Cu_dist, x)=other_size_dist(field,number_of_bins)
		
		prob_dist_matrix[count,:]=prob_dist
		Cu_dist_matrix[count,:]=Cu_dist
		x_matrix[count,:]=x

	if save_mat_file==True:
		save_prob_dist_mat_file(split_distribution_by_mass,prob_dist_matrix, Cu_dist_matrix, x_matrix ,split_prob_matrix,\
				split_Cu_matrix,input_folder,start_year, end_year,constraints,field_name )
	
	return [prob_dist_matrix, Cu_dist_matrix, x_matrix ,split_prob_matrix, split_Cu_matrix]


def save_prob_dist_mat_file(split_distribution_by_mass,prob_dist_matrix, Cu_dist_matrix, x_matrix ,split_prob_matrix, split_Cu_matrix,input_folder,start_year, end_year,constraints,field_name ):
	if split_distribution_by_mass==False:
		split_prob_matrix=0
		split_Cu_matrix=0
   
	exp_name=str.split(input_folder,'_bergs_')[-1]

	#Creating filename
	mat_filename='processed_data/'+exp_name+'_size_distribution_'+str(start_year)+'_to_'+ str(end_year)

	if field_name!='Area':
		mat_filename=mat_filename+'_'+field_name

	if constraints!=None:
		if len(constraints)>0:
			mat_filename=mat_filename+'_constraints'
			for k in range(len(constraints)):
				current_constraint=constraints[k]
				constraint_name=current_constraint['Constraint_field_name']
				lower_bound=current_constraint['lower_bound']
				upper_bound=current_constraint['upper_bound']
				original_values=current_constraint['original_values']
				mat_filename=mat_filename+ '_'+constraint_name + '_' + str(abs(upper_bound)) +'_to_'+ str(abs(lower_bound))

	mat_filename=mat_filename +'.mat'

	sc.savemat(mat_filename, {'prob_dist_matrix':prob_dist_matrix,'Cu_dist_matrix':Cu_dist_matrix,'x_matrix':x_matrix,\
		'split_prob_matrix':split_prob_matrix, 'split_Cu_matrix':split_Cu_matrix }) 
	
	print 'File: ' + mat_filename + ' saved.'




def plot_power_laws_on_figure(x,plot_1_2_line, plot_2_3_line, plot_3_2_line, plot_4_3_line):
	if plot_1_2_line==True:
		y=x**(-0.5)
		total=1e6
		y=y/total
		#plt.plot(x,y, linewidth=1.5, linestyle="-",color='black',label='$y=x^{-0.5}$')

	if plot_2_3_line==True:
		y=x**(-0.667)
		total=1e4
		y=y/total
		plt.plot(x,y, linewidth=1.5, linestyle="--",color='black',label='$y=x^{-0.667}$')

	if plot_3_2_line==True:
		y=x**(-3./2.)
		total=0.0577
		y=y/total
		plt.plot(x,y, linewidth=1.5, linestyle=":",color='black',label='$y=x^{-3/2}$')
	
	#if plot_4_5_line==True:
	#	y=x**(-4./5.)
	#	total=1e3
	#	y=y/total
	#	plt.plot(x,y, linewidth=1.5, linestyle=":",color='black',label='$y=x^{-4/5}$')
	
	if plot_4_3_line==True:
		y=x**(-4./3.)
		total=0.0877
		total=1e2
		y=y/total
		plt.plot(x,y, linewidth=1.5, linestyle="-",color='black',label='$y=x^{-4/3}$')

###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
def main():
	#Clear screen
	#os.system('clear')



	root_path='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg'
	#run_names=np.array(['tournadre','freq','all_big', 'mass', 'all_small']) ; naming_flag='groups'
	#run_names=np.array(['freq','all_big', 'mass', 'all_small']) ; naming_flag='groups'
	run_names=np.array(['Delta1', 'Delta2', 'Delta3', 'Delta4', 'Delta5','Delta6', 'Delta7', 'Delta8', 'Delta9', 'Delta10']) ; naming_flag='Delta'
	run_names=np.array(['rolling_Delta1_Burton_new','rolling_Delta5_Burton_new','rolling_Delta9_Burton_new',\
			'rolling_tournadre_Burton_breaking_aspect_fixed','rolling_tournadre_Burton_breaking_aspect_fixed'])
	run_names=np.array(['rolling_Delta9_Burton_new', 'rolling_tournadre_Burton_new','Delta9','tournadre',\
			'rolling_tournadre_Burton_breaking_D_fixed','rolling_tournadre_Burton_breaking_aspect_fixed'])
	run_names=np.array(['tournadre_C2'])

	#Flags
	find_probability_dist=False
	plot_average_dist=0
	plot_full_time_series=0
	plot_cumulative_ditribution=False
	plot_prob_ditribution=False
	save_mat_file=False
	split_distribution_by_mass=True   # Splits the distribution up into mass classes

	plot_iceberg_numbers=True  #Only use when all others are false

	#Parameters and variables
	start_year=1900 
	#end_year0=2019
	end_year0=2050
	number_of_bins=100
	Number_of_years=150
	#Number_of_years=25
	#Number_of_years=15
	#Number_of_years=10
	col_num=1


	#Include contraints to use to filter the data
	# The constraint has the form [field_name, lower_bonud, upper_bound]
	#constraint1=np.array(['lat',-65.,-10.,'not_original'])  # Latituds north of -60 in the southern hemisphere

	constraint_mass = {'Constraint_field_name': 'mass', 'lower_bound': 3.8e11 , 'upper_bound': 3.95e11, 'original_values': True}
	constraint_small_mass = {'Constraint_field_name': 'mass', 'lower_bound': 0.0e8 , 'upper_bound': 4.0e8, 'original_values': False}
	constraint_big_mass = {'Constraint_field_name': 'mass', 'lower_bound': 1.0e11 , 'upper_bound': 1.0e15, 'original_values': False}
	constraint_lat = {'Constraint_field_name': 'lat', 'lower_bound': -70 , 'upper_bound': -10, 'original_values': False}
	constraint_depth = {'Constraint_field_name': 'depth', 'lower_bound': 1000 , 'upper_bound': 8000, 'original_values': False}



	constraint_dist_from_calving = {'Constraint_field_name': 'distance_from_calving', 'lower_bound': 50*1000 , 'upper_bound': 10000000000000, 'original_values': False}


	#constraint2=np.array(['lon',-60.,0.,'not_original'])  #Longitude of Weddel Sea    , lon appears to go from -270 to 90. 
	constraints=[]
	constraints.append(constraint_small_mass)
	#constraints.append(constraint_mass)
	#constraints.append(constraint_lat)
	#constraints.append(constraint_depth)
	#constraints.append(constraint_dist_from_calving)
	#contraints


	mass_scaling=np.array([2000, 200, 50, 20, 10, 5, 2, 1, 1, 1])
	initial_mass=np.array([8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11]) 



	#field_name='area';  x_min=1;x_max=1.e9  # Good for area including all Deltas
	field_name='area';  x_min=1;x_max=1.e9  # Good for area 
	field_name='mass' ; x_min=100;x_max=1.e12  # Good for mass
	#Defining a list of colors
	color_vec=np.array(['blue', 'red','purple','green', 'coral', 'cyan', 'magenta','orange', 'black', 'grey', 'yellow', 'orchid', 'blue', 'red','purple','green', 'coral', 'cyan' ])

	fig = plt.figure(1)
	#for k in np.array([0,1,2,3,4]):
	#for k in np.array([0,1,2,3]):
	#for k in np.array([4,5,8]):
	#for k in np.array([6,9]):
	for k in np.array([0]):
	#for k in range(9):
		input_folder=root_path + '_bergs_' + run_names[k]

		#Make sure that data exists of the years allocated, othersize change it.
		(min_year, max_year)=find_max_min_years(input_folder,'.iceberg_trajectories.nc')
		end_year=min(max_year,end_year0)
		start_year=max(end_year-Number_of_years,min_year)
		
		if max_year>=start_year:

			#Getting the distributions from the data
			print input_folder
			print start_year
			if find_probability_dist is True:
				(prob_dist_matrix, Cu_dist_matrix, x_matrix, split_prob_matrix, split_Cu_matrix) = generate_dist_timeseries_matrix(input_folder, field_name,start_year, end_year, \
						number_of_bins,x_min,x_max,constraints,split_distribution_by_mass,save_mat_file)
			else:
				(Num_bergs_all, Num_bergs_constraints,Num_bergs_splits) = get_number_of_icebergs(input_folder, field_name,start_year, end_year, \
					number_of_bins,x_min,x_max,constraints,split_distribution_by_mass,save_mat_file)



			#Plotting the distributions
			#plot_distributions(start_year, end_year, plot_full_time_series, plot_average

			if plot_full_time_series==1:
				count=-1
				for year in range(start_year, end_year+1):
					count=count+1
					Cu_dist=Cu_dist_matrix[count,:]
					x=x_matrix[count,:]
					prob_dist=prob_dist_matrix[count,:]
					subplot(1,2,1)
					plt.plot(x,Cu_dist,color=color_vec[k], linewidth=1.5, linestyle="-")
					subplot(1,1,1)
					plt.plot(x,prob_dist,color=color_vec[k], linewidth=1.5, linestyle="-")

			if plot_average_dist==1:
				label= str.split(input_folder,'_')[-1]
				Cu_dist=np.squeeze(np.mean(Cu_dist_matrix, axis=0))
				prob_dist=np.squeeze(np.mean(prob_dist_matrix, axis=0))
				x=np.squeeze(np.mean(x_matrix, axis=0))
				if plot_cumulative_ditribution==True:
					ax = fig.add_subplot(1,2,1)
					plt.plot(x,Cu_dist, color=color_vec[k], linewidth=5.,label=label)
					#plt.plot(x,Cu_dist, color='black', linewidth=5.,label=label)
					col_num=2
				if plot_prob_ditribution==True:
					ax = fig.add_subplot(1,col_num,col_num)
					plt.plot(x,prob_dist, color=color_vec[k], linewidth=5.,label=label)
					#plt.plot(x,prob_dist, color='black', linewidth=5.,label=label)

				#Plotting the split up distributions
				if split_distribution_by_mass==True:
					print 'you are here'
					M=split_prob_matrix.shape
					for mass_num in range(M[2]):
						label= 'D'+str(mass_num+1)
						Cu_dist=np.squeeze(np.mean(split_Cu_matrix[:,:,mass_num], axis=0))
						prob_dist=np.squeeze(np.mean(split_prob_matrix[:,:,mass_num], axis=0))
						x=np.squeeze(np.mean(x_matrix, axis=0))
						if plot_cumulative_ditribution==True:
							ax = fig.add_subplot(1,2,1)
							plt.plot(x,Cu_dist, color=color_vec[mass_num], linewidth=1.,label=label)
						if plot_prob_ditribution==True:
							ax = fig.add_subplot(1,col_num,col_num)
							plt.plot(x,prob_dist, color=color_vec[mass_num], linewidth=1.,label=label)

	if plot_cumulative_ditribution==True:
		ax = fig.add_subplot(1,2,1)
		ax.set_xscale('log')
		ax.set_ylim([0, 1])
		plt.xlabel(field_name)
		plt.ylabel('Cumulative Probability')
		#plt.legend(loc='upper left', frameon=True)

	if plot_prob_ditribution==True:
		ax = fig.add_subplot(1,col_num,col_num)
		plt.ylim([1e-7, 1e-2]) #For area
		plt.ylim([1e-10, 1e-2]) #For area, with constraints
		plt.ylim([1e-13, 1e-3]) #For area, with constraints
		plt.xlim([1e1, 1e9]) #For area, with constraints
		#plt.ylim([1e-9, 1e-2]) #For area, all deltas
		plt.xlim([1e1, 1e13]) #For mass
		plt.ylim([1e-13, 1e-4]) #For mass
		ax.set_xscale('log')
		ax.set_yscale('log')
		plt.xlabel(field_name)
		#plt.plot(x,x**(-1.5), linewidth=1.5, linestyle=":")
		plot_1_2_line=True; plot_2_3_line=True; plot_3_2_line=True; plot_4_3_line=True 
		plot_power_laws_on_figure(x,plot_1_2_line, plot_2_3_line, plot_3_2_line, plot_4_3_line,)
		#plt.plot(x,1e4*x**(-1.5), linewidth=1.5, linestyle=":",color='black',label="x^{-1.5}")
		#plt.plot(x,1e4*x**(-2/3), linewidth=1.5, linestyle="--",color='black',label="x^{-0.66}")
		#plt.plot(x,1e4*x**(-1/2), linewidth=1.5, linestyle=":",color='black',label="x^{-0.5}")
		plt.ylabel('Probabilty')
		plt.legend(loc='upper right', frameon=True)
	
	if plot_iceberg_numbers==True:
		ax = fig.add_subplot(2,1,1)
		plt.plot(Num_bergs_all,'k')
		plt.ylabel('Number of icebergs')
		plt.xlabel('Year')
		#plt.plot(Num_bergs_constraints,'r')
		ax = fig.add_subplot(2,1,2)
		for k in range(10):
			if k>0:
				label=str(initial_mass[k-1]/(1*10**9))+ 'Gt to ' + str(initial_mass[k]/(1*10**9)) + 'Gt' 
			else:
				label=str(0.0)+ 'Gt to ' + str(initial_mass[k]/(1*10**9)) +'Gt' 
			plt.plot(Num_bergs_splits[:,k],color_vec[k],label=label)
		plt.ylabel('Number of icebergs')
		plt.xlabel('Year')
		#ax = fig.add_subplot(2,1,2)
		#plt.plot(Num_berg_spl)
		#ax.set_xscale('log')
		#ax.set_ylim([0, 1])
		#plt.xlabel(field_name)
		plt.legend(loc='lower right', frameon=True)


	#plt.legend(loc='upper right', frameon=True)
	fig = matplotlib.pyplot.gcf()
	fig.set_size_inches(9,4.5)
	plt.show()


	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

