#!/usr/bin/env python
from netCDF4 import Dataset
import numpy as np
from PIL import *
from pylab import *
import math
import os
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pdb
import argparse
import sys
from scipy.interpolate import interp1d
import scipy.io as sc
import scipy.stats as sp
from iceberg_mass_comparisons import define_paths_array
from iceberg_mass_comparisons import find_max_min_years
from distributions_of_bergs import get_valid_data_from_traj_file
from distributions_of_bergs import add_constraint_to_data
from sea_ice_concentrations import MidpointNormalize



def plot_age_distribution(values):
	ages = np.logspace(np.log10(0.1), np.log10(12), num=20, endpoint=True)
	#hist1, ages =np.histogram(values, bins=ages)
	#ages=ages[:-1]
	fig = plt.figure(1)
	ax = fig.add_subplot(1,1,1)

	#Fitting the distribution
	param=sp.lognorm.fit(values)
	ages=np.linspace(0,12,50)
	pdf_fitted = sp.lognorm.pdf(ages, param[0], loc=param[1], scale=param[2])
	#plt.plot(ages,pdf_fitted,'r-',ages,hist1,'g-')
	plt.plot(ages,pdf_fitted,'r-')
	#shape, loc, scale = sc.stats.lognorm.fit(hist, floc=0)

	mu=np.log(param[2])
	sigma=param[0]
	#plt.plot(ages,hist1)
	print param[0],param[1],param[2]
	print mu, sigma
	#plt.hist(ages,hist1)
	plt.hist(values,bins=30,normed=True)
	plt.show()



def construct_occurance_matrix(x,y,Matrix,x_vec, y_vec):
	#age is y   (field1)
	#mass is x  (field2)
	for i in range(len(x_vec)-1):
		for j in range(len(y_vec)-1):
			x_valid=(x>x_vec[i])*(x<x_vec[i+1])
			y_valid=(y>y_vec[j])*(y<y_vec[j+1])
			values=y[np.where(x_valid*y_valid)]
			values=values[np.where(np.isnan(values)==0.)]
			Matrix[i,j]=Matrix[i,j]+len(values)

	#Handeling the end points
	return Matrix
			
		
	

def expected_values(x,y,number_of_bins,min_bin_value=None):
	#finding expected value of x in bins of y
	plot_distribution=False

	#age is y   (field1)
	#mass is x  (field2)

	if min_bin_value is not None:	
		x1_min=max(10**5,min(x)) ; x1_max=max(x)   #Use with mass vs age
		#x1_min=10**5 ; x1_max=10**12
	else:
		x1_min=min(x) ; x1_max=max(x)
	iceberg_age=0
	iceberg_mass=0
	#number_of_bins=300
	x1 = np.logspace(np.log10(x1_min), np.log10(x1_max), num=number_of_bins, endpoint=True)
	x1[0]=0.
	#x1 = np.linspace(x1_min, x1_max, num=number_of_bins, endpoint=True)
	
	x_bins=x1[:-1]
	
	y_ave=x1[:-1]*0
	y_std=x1[:-1]*0
	y_median=x1[:-1]*0
	y_mode=x1[:-1]*0
	sigma=x1[:-1]*0
	loc=x1[:-1]*0
	scale=x1[:-1]*0
	mu=x1[:-1]*0
	count=-1
	for k in x_bins:
		count=count+1
		valid=(x>x1[count])*(x<x1[count+1])
		#valid=(x>x1[k])*(x<x1[k+1])

		values=y[np.where(valid)]
		values=values[np.where(np.isnan(values)==0.)]
		y_ave[count]=np.mean(values)
		y_std[count]=np.std(values)
		
		#Log normal fitting
		param=sp.lognorm.fit(values)
		sigma[count]=param[0] #Shape
		loc[count]=param[1] #Loc
		scale[count]=param[2] #Scale
		mu[count]=np.log(param[2])
		
		#median and mode
		y_median[count]=np.median(values)
		#values2=np.round(values,1)
		#print values2
		#y_mode[count]=max(set(values2), key=values2.count)
		#y_mode[count]=sp.mode(values2[:])

		if plot_distribution==True:
			Bin_number=30
			if count==Bin_number:
				plot_age_distribution(values)
	return [x_bins,y_ave,y_std,sigma,scale,loc,mu,y_median,y_mode]


def derivative(x,y):
	N=len(x)
	first_ind=np.insert(np.arange(1,N),N-2 ,N-1 )
        second_ind=np.insert(np.arange(0,N-1),0 ,0 )
	dy_dx=(y[first_ind]-y[second_ind])/(x[first_ind]-x[second_ind])
	
	return dy_dx

def fit_derivative(x,y): #y is age, x in mass
	max_age=10
	x1=x[np.where((y<max_age) *(~np.isnan(y)))]
	y1=y[np.where((y<max_age) *(~np.isnan(y)))]

	fit_with_interp=True
	if fit_with_interp==True:
		y2=np.linspace(0.05,max(y1),30)
		set_interp = interp1d(y1, x1, kind='linear')
		x2 = set_interp(y2)
		x_new=x2 ; y_new=y2
	else:
		x_new=x1 ; y_new=y1

	#Fitting x=p(y)
	N=10 #polynomial degree
	p=np.polyfit(y_new,x_new,N)
	p_der=np.polyder(p)
	x_fit=np.polyval(p,y_new)
	#y_fit

	dy_dx=np.polyval(p_der,y_new)  #dm/dt
	#dy_dx=derivative(y_new,x_new)

	return [dy_dx, x_new, y_new, x_fit]

def find_average_over_time_without_nan(matrix,ave_type='median'):
        M= matrix.shape
        Total=np.zeros([1,M[1]])
	for j in range(M[1]):
		temp=np.squeeze(matrix[:,j])
		temp=temp[np.where(np.isnan(temp)==0)]
		if ave_type=='median':
			Total[0,j]=np.median(temp)
		elif ave_type=='mean':
			Total[0,j]=np.mean(temp)
	Total=np.squeeze(Total)
        return Total


def load_mass_vs_age_from_mat_file(filename,load_type='lognormal'):

        mat_contents=sc.loadmat(filename)
        Total_x=np.squeeze(mat_contents['Total_x'][:][0])

        if load_type=='occurance':
                occurance_x=np.squeeze(mat_contents['occurance_x'][:])
                occurance_y=np.squeeze(mat_contents['occurance_y'][:])
                occurance_matrix=mat_contents['occurance_matrix'][:]
                print occurance_x.shape
                print occurance_y.shape
                print occurance_matrix.shape
                return [ occurance_x, occurance_y, occurance_matrix]
        if load_type=='median':
                Total_median=np.squeeze(mat_contents['Total_median'][:][0])
                Total_mean=np.squeeze(mat_contents['Total_mean'][:][0])
                return [Total_median,  Total_mean,  Total_x ]
        if load_type=='lognormal':
                Total_sigma=mat_contents['Total_sigma'][0]
                Total_scale=mat_contents['Totla_scale'][0]
                return [Total_scale, Total_sigma,Total_x]

        if load_type=='normal':
                #mean_matrix=mat_contents['mean_matrix']
                std_matrix=mat_contents['std_matrix'][0]
                x_matrix=mat_contents['std_matrix'][0]
                #return [mean_matrix, std_matrix]
                Total_mean=mat_contents['Total_mean'][0]
                Total_std=mat_contents['Total_std'][0]
                return [Total_mean, Total_std,Total_x]




def save_mat_file(mean_matrix, std_matrix, x_matrix, Total_mean, Total_std, Total_x, start_year,end_year,Delta_name,field_name1,field_name2,\
		sigma_matrix, scale_matrix, loc_matrix, mu_matrix, Total_sigma, Total_loc, Total_scale, Total_mu,Total_median,constraints,
		occurance_x,occurance_y,occurance_matrix):
	mat_filename= Delta_name + '_' + field_name1 + '_vs_' + field_name2 + '_' + str(start_year) +'_to_' + str(end_year) + '_expected_values'
	
	if constraints!=None:
		if len(constraints)>0:
			mat_filename=mat_filename+'_constraints'
			for k in range(len(constraints)):
				current_constraint=constraints[k]
				constraint_name=current_constraint['Constraint_field_name']
				lower_bound=current_constraint['lower_bound']
				upper_bound=current_constraint['upper_bound']
				mat_filename=mat_filename+ '_'+constraint_name + '_' + str(abs(lower_bound)) +'_to_'+ str(abs(upper_bound))
	
	mat_filename='processed_data/' + mat_filename +'_scatter.mat'

	sc.savemat(mat_filename, {'mean_matrix':mean_matrix,'std_matrix':std_matrix,'x_matrix':x_matrix,\
			'Total_mean':Total_mean,'Total_std':Total_std,'Total_x':Total_x ,\
			'sigma_matrix':sigma_matrix, 'scale_matrix':scale_matrix, 'loc_matrix':loc_matrix, 'mu_matrix':mu_matrix,\
			'Total_sigma':Total_sigma, 'Total_loc': Total_loc, 'Totla_scale':Total_scale, 'Total_mu':Total_mu, 'Total_median':Total_median,\
			'occurance_x':occurance_x ,'occurance_y':occurance_y,'occurance_matrix':occurance_matrix}) 
	print 'File: ' + mat_filename + ' saved.'

###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
def main():
	#Clear screen
	#os.system('clear')

	#Defining possible paths
	all_paths=define_paths_array()

	#Flags
	plot_full_time_series=0
	plot_scatter_plot=False
	plot_mean_fields_again=True
	plot_occurance_matrix=False
	plot_derivatives_and_fits=False
	plot_rolling_equation=False
	find_expected_values=True
	save_the_data=False
	load_the_data=True

	#Parameters and variables
	#start_year=1980
	end_year0=2010
	#end_year0=1999
	Number_of_years=3
	#Number_of_years=0

	#Which fields do you want to compare?
	#field_name1='area'
	#field_name1='age'
	#field_name1='side_decay'
	#field_name1='Me'
	#field_name1='sst'
	#field_name1='Ss'
	field_name1='cn'
	#field_name1='speed'
	#field_name1='Mb'
	#field_name1='mass'
	#field_name1='length'
	#field_name1='width'
	#field_name1='thickness'
	#field_name1='LW_max'
	#field_name1='distance_from_calving'
	
	field_name2='mass'
	#field_name2='area'
	#field_name2='age'
	#field_name2='thickness'
	#field_name2='length'
	#field_name2='width'

	initial_mass_vec=np.array([ 8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11])
	#Include contraints to use to filter the data
	
	constraint2 = {'Constraint_field_name': 'lon', 'lower_bound': -150 , 'upper_bound': -65, 'original_values': True} #Longitude of AB
	constraint_SH = {'Constraint_field_name': 'lat', 'lower_bound': -150 , 'upper_bound': -5, 'original_values': True} #Longitude of AB
	constraint_NH = {'Constraint_field_name': 'lat', 'lower_bound': 5 , 'upper_bound': 100, 'original_values': True} #Longitude of AB
	constraint_birth_year = {'Constraint_field_name': 'year', 'lower_bound': 1899.9 , 'upper_bound': 1900.1, 'original_values': True}
	constraint_depth = {'Constraint_field_name': 'depth', 'lower_bound': 500 , 'upper_bound': 8000, 'original_values': False}
	constraint_dist_from_calving = {'Constraint_field_name': 'distance_from_calving', 'lower_bound': 1000 , 'upper_bound': 10000000000000, 'original_values': False}
	#constraint_mass = {'Constraint_field_name': 'mass', 'lower_bound': -10**16 , 'upper_bound': -0.5*(10**8), 'original_values': False, 'subtract_original': True}
	constraint_mass = {'Constraint_field_name': 'mass', 'lower_bound': 1.0*(10**10) , 'upper_bound': 4.0*(10**15), 'original_values': False, 'subtract_original': False}
	constraint_mass0 = {'Constraint_field_name': 'mass', 'lower_bound': 7.3*(10**11) , 'upper_bound': 7.6*(10**11), 'original_values': True, 'subtract_original': False}
	#constraint_mass0 = {'Constraint_field_name': 'mass', 'lower_bound': 3.2*(10**9) , 'upper_bound': 3.4*(10**9), 'original_values': True, 'subtract_original': False}
	constraint_day = {'Constraint_field_name': 'day', 'lower_bound': 92 , 'upper_bound': 94., 'original_values': True, 'subtract_original': False}
	constraint_vvel = {'Constraint_field_name': 'vvel', 'lower_bound': 0.0001 , 'upper_bound': 100, 'original_values': False, 'subtract_original': False}
	constraint_uvel = {'Constraint_field_name': 'uvel', 'lower_bound': -100 , 'upper_bound': -0.00100, 'original_values': False, 'subtract_original': False}
	constraint_age = {'Constraint_field_name': 'age', 'lower_bound': 0 , 'upper_bound': 5, 'original_values': False} 
        constraint_lon_WestGreenland = {'Constraint_field_name': 'lon', 'lower_bound': -45 , 'upper_bound': -20, 'original_values': True}
        constraint_lat_SouthGreenland = {'Constraint_field_name': 'lat', 'lower_bound': 50.0 , 'upper_bound':66.0, 'original_values': True}

	#constraints.append(constraint2)


	#if field_name1=='age' or field_name2=='age':
	#	print 'applying age constraint'
	#	day_constraint = {'Constraint_field_name': 'day', 'lower_bound': 0 , 'upper_bound': 258, 'original_values': True} #Longitude of AB
	constraints=[]	
	#constraints.append(constraint_dist_from_calving)
	#constraints.append(constraint_depth)
	constraints.append(constraint_mass0)
	#constraints.append(constraint_day)
	#constraints.append(constraint_age)
	#constraints.append(constraint_SH)
	#constraints.append(constraint_birth_year)
	#constraints.append(constraint_NH)
	#constraints.append(constraint_vvel)
	#constraints.append(constraint_uvel)
	#	constraints.append(day_constraint)
	#constraints.append(constraint_lon_WestGreenland)
	constraints.append(constraint_lat_SouthGreenland)



	mass_scaling=np.array([2000, 200, 50, 20, 10, 5, 2, 1, 1, 1])



	#field_name='area';  x_min=1;x_max=1.e9  # Good for area including all Deltas
	field_name='mass';  x_min=1;x_max=1.e7  # Good for area 
	#field_name='mass' ; #x_min=100;x_max=1.e12  # Good for mass
	#Defining a list of colors
	color_vec=np.array(['blue', 'red','purple','green', 'coral', 'cyan', 'magenta','orange', 'black', 'grey', 'yellow', 'orchid', 'blue', 'red','purple','green', 'coral', 'cyan' ])
	root_path='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_'
	
	#Delta_list=np.array(['Delta9', 'rolling_tournadre_WM', 'rolling_Delta9_Burton_new'])
	#Delta_list=np.array(['rolling_Delta1_Burton_new','rolling_Delta5_Burton_new','rolling_Delta9_Burton_new'])
	#Delta_list=np.array(['rolling_tournadre_Burton_breaking_aspect_fixed','rolling_tournadre_Burton_breaking_D_fixed','rolling_tournadre_Burton_new'])
	Delta_list=np.array(['tournadre','rolling_tournadre_WM','rolling_tournadre_Burton_new'])
	#Delta_list=np.array(['rolling_tournadre_Burton_new/trajectories/','rolling_tournadre_Rolling_off/trajectories/'])
	Delta_list=np.array(['rolling_tournadre_Burton_new','rolling_tournadre_Rolling_off'])
	#Delta_list=np.array(['Delta1','Delta9'])
	#Delta_list=np.array(['Delta2','Delta3', 'Delta4','Delta5','Delta6'])

	#Rolling_file='processed_data/rolling_tournadre_Burton_new_Me_vs_mass_1960_to_2010_expected_values_constraints_mass_7.3e+11_to_7.6e+11_lat_50.0_to_75.0_scatter.mat'
	#No_Rolling_file='processed_data/rolling_tournadre_Rolling_off_Me_vs_mass_1960_to_2010_expected_values_constraints_mass_7.3e+11_to_7.6e+11_lat_50.0_to_75.0_scatter.mat'
	
	No_Rolling_file='processed_data/rolling_tournadre_Rolling_off_Me_vs_mass_1960_to_2010_expected_values_constraints_mass_7.3e+11_to_7.6e+11_lon_45_to_20_lat_50.0_to_75.0_scatter.mat'
	Rolling_file='processed_data/rolling_tournadre_Burton_new_Me_vs_mass_1960_to_2010_expected_values_constraints_mass_7.3e+11_to_7.6e+11_lon_45_to_20_lat_50.0_to_75.0_scatter.mat'
	load_filename_list=np.array([Rolling_file,No_Rolling_file])
	
	#Name_list=np.array(['Martin and Adcroft','Weeks and Mellor','Burton'])
	Name_list=np.array(['Rolling', 'No Rolling'])

	fig = plt.figure(1)
	#for k in range(13):#for k in np.array([3]):
	count1=0
	print 'field1 = ', field_name1, ', field2 = ' , field_name2
	if load_the_data is True:
		for k in np.array([0,1]):
		#for k in np.array([0,1,2,3,4]):
		#for k in np.array([8,5,3,1]):
		#for k in np.array([9,5,1]):
		#for k in np.array([4,5,6,9]):
		#for k in np.array([7,8,10]):
		#for k in np.array([10,9,8,7,6,5,4,3,2,1]):
		#for k in range(1,9):
			#ax = fig.add_subplot(2,1,k+1)
			ax = fig.add_subplot(1,1,1)
			count1=count1+1
			#input_folder=all_paths[k]
			#Delta_name='Delta'+str(k)
			Delta_name=Delta_list[k]
			input_folder=root_path+Delta_name 
			input_folder = input_folder + '/trajectories/'
			print input_folder
			
			#Make sure that data exists of the years allocated, othersize change it.
			(min_year, max_year)=find_max_min_years(input_folder,'.iceberg_trajectories.nc')
			end_year=min(max_year,end_year0)
			start_year=max(end_year-Number_of_years,min_year)

			number_of_bins=100
			mean_matrix=np.zeros((end_year-start_year+1,number_of_bins-1))
			std_matrix=np.zeros((end_year-start_year+1,number_of_bins-1))
			x_matrix=np.zeros((end_year-start_year+1,number_of_bins-1))
			
			sigma_matrix=np.zeros((end_year-start_year+1,number_of_bins-1))
			scale_matrix=np.zeros((end_year-start_year+1,number_of_bins-1))
			loc_matrix=np.zeros((end_year-start_year+1,number_of_bins-1))
			mu_matrix=np.zeros((end_year-start_year+1,number_of_bins-1))
			median_matrix=np.zeros((end_year-start_year+1,number_of_bins-1))
			mode_matrix=np.zeros((end_year-start_year+1,number_of_bins-1))

			number_of_age_bins=40
			occurance_y= np.logspace(np.log10(10**5), np.log10(10**11), num=number_of_bins, endpoint=True)
			occurance_x = np.linspace(0, 3, num=number_of_age_bins, endpoint=True)
			#occurance_x = np.linspace(0, 40, num=number_of_age_bins, endpoint=True)
			occurance_matrix=np.zeros((len(occurance_x),len(occurance_y)))

			year_count=-1
			for year in range(start_year,end_year+1):
				year_count=year_count+1
				#year=1945
				filename ='/' + str(year) + '0101.iceberg_trajectories.nc'
				input_file=input_folder + filename
				print input_file
				
				field1 = get_valid_data_from_traj_file(input_file,field_name1,subtract_orig=False)
				field2 = get_valid_data_from_traj_file(input_file,field_name2,subtract_orig=False)

				print 'Lendth of field: ' , len(field1)
				#Getting the distributions from the data

				#Handeling constraints
				field1=add_constraint_to_data(input_file,field1,constraints)
				field2=add_constraint_to_data(input_file,field2,constraints)
				print 'Lendth of field after constraints: ' , len(field1)
				
				#field2=field2*0.81
				

				#if create_2D_histogram==True:


				
				if find_expected_values==True:
					(x_bins,y_ave,y_std,sigma,scale,loc,mu,y_median,y_mode)=expected_values(field2,field1,number_of_bins, min_bin_value=True)
					mean_matrix[year_count,:]=y_ave
					std_matrix[year_count,:]=y_std
					x_matrix[year_count,:]=x_bins
					
					sigma_matrix[year_count,:]=sigma
					scale_matrix[year_count,:]=scale
					loc_matrix[year_count,:]=loc
					mu_matrix[year_count,:]=mu

					median_matrix[year_count,:]=y_median
					mode_matrix[year_count,:]=y_mode
					

					#Dealing with occurance
					Occurance_matrix=construct_occurance_matrix(field1,field2,occurance_matrix, occurance_x, occurance_y)


				if plot_scatter_plot==True:
					#Plotting the distributions
					#ax = fig.add_subplot(3,1,k+1)
					#ax = fig.add_subplot(2,1,count1)
					label=Name_list[k]
					#label='Delta' + str(k)
					plt.plot(field1,field2,'.',color=color_vec[k],label=label)
					plt.xlabel(field_name1)
					plt.ylabel(field_name2)
					#plt.ylim([-75., -40]) #For mass
					#plt.xlim([0, 30]) #For mass
					#ax.set_yscale('log')
					#ax.set_xscale('log')
					plt.title(Name_list[k])
				

				if plot_derivatives_and_fits==True:
					#dy_dx=derivative(x_bins,y_ave)
					#Note:  x_bins is mass, y_ave is age
					(dy_dx, x_new, y_new, x_fit)=fit_derivative(x_bins,y_ave)

					ax = fig.add_subplot(1,1,1)
					label='Delta' + str(k+1)

					plt.plot(y_ave,x_bins,color=color_vec[k],label=label,linewidth=2.)
					plt.plot(y_ave,x_bins,'*',color=color_vec[k])
					plt.plot(y_new,x_new,color='red',label=label,linewidth=4.)
					plt.plot(y_new,x_fit,color='magenta',label=label,linewidth=2.)


					plt.plot(y_ave+y_std,x_bins,color='green')
					plt.plot(y_ave-y_std,x_bins,color='green')
					plt.xlabel(field_name1)
					plt.ylabel(field_name2)
					#ax.set_yscale('log')
					#ax.set_xscale('log')
					
					ax = fig.add_subplot(1,2,2)
					plt.plot(x_bins,dy_dx,color=color_vec[k],label=label,linewidth=4.)
					plt.plot(x_bins,dy_dx,'*')
					plt.plot(x_new,dy_dx)
					plt.ylabel('dy/dx')
					plt.xlabel('Mass')
					#ax.set_xscale('log')
					#ax.set_yscale('log')


			if plot_rolling_equation==True:
				D=np.linspace(0.1,300,100)
				L=sqrt(0.92*(D**2)+(58.32*D))
				plt.plot(L,D,linestyle='--')
				#plt.plot(L,L,linestyle=':')
				plt.plot(L,0.81*L,linestyle=':')
				#plt.plot(L,L/1.2,linestyle=':')
				#plt.plot(L,L-(,linestyle=':')

			
			#Normalizing
			print np.sum(occurance_matrix)
			occurance_matrix=occurance_matrix/(np.sum(occurance_matrix))
			
			if plot_occurance_matrix is True:
				vmax=np.max(occurance_matrix)
				#cNorm = MidpointNormalize(vmin=0, vmax=5,midpoint=0)
				cNorm = mpl.colors.Normalize(vmin=0, vmax=vmax)
				cNorm = mpl.colors.LogNorm(vmin=10**-5, vmax=0.1)
				print len(occurance_x),len(occurance_y)
				plt.pcolor(occurance_x, occurance_y, occurance_matrix.transpose(),cmap='jet', norm=cNorm)
				plt.colorbar()
				plt.xlabel(field_name1)
				plt.ylabel(field_name2)
				ax.set_yscale('log')
				ax.set_xscale('log')


			if plot_mean_fields_again==True:
				Total_median=find_average_over_time_without_nan(median_matrix,'median')
				Total_mean=find_average_over_time_without_nan(mean_matrix,'mean')
				Total_std=find_average_over_time_without_nan(np.sqrt(std_matrix**2),'mean')
				Total_x=find_average_over_time_without_nan(x_matrix,'mean')
				Total_sigma=find_average_over_time_without_nan(sigma_matrix,'mean')
				Total_scale=find_average_over_time_without_nan(scale_matrix,'mean')
				Total_loc=find_average_over_time_without_nan(loc_matrix,'mean')
				Total_mu=find_average_over_time_without_nan(mu_matrix,'mean')
				if save_the_data==True:
					save_mat_file(mean_matrix, std_matrix, x_matrix, Total_mean, Total_std, Total_x, start_year,end_year,Delta_name,field_name1,field_name2,\
							sigma_matrix, scale_matrix, loc_matrix, mu_matrix, Total_sigma, Total_loc, Total_scale, Total_mu,Total_median,constraints,\
							occurance_x,occurance_y,occurance_matrix)

				#ax = fig.add_subplot(1,1,1)
				#label='Delta' + str(k+1)
				label=Delta_list[k]

				#These should be uncommented
				#plt.plot(Total_mean,Total_x,color='black',label=label,linewidth=4)
				#plt.plot(Total_median,Total_x,color='green',label=label,linewidth=4)
				#plt.plot(Total_mean+Total_std,Total_x,color='magenta',linewidth=4)
				#plt.plot(Total_mean-Total_std,Total_x,color='magenta',linewidth=4)

				#plt.plot(Total_median,Total_x,color=color_vec[k],label=label,linewidth=4)
				print Total_median
				plt.title(Name_list[k])
				M=median_matrix.shape
				#for loop1 in range(M[0]):
				#	plt.plot(np.squeeze(median_matrix[loop1,:]),Total_x,color=color_vec[loop1+3],label=label,linewidth=1)
				plt.plot(Total_median,Total_x,color=color_vec[k],label=label,linewidth=4)
				plt.xlabel(field_name1)
				plt.ylabel(field_name2)
				ax.set_yscale('log')
				#ax.set_xscale('log')
	else:
		if plot_mean_fields_again==True:
				ax = fig.add_subplot(1,1,1)
				for k in range(2):
					print k	
					filename=load_filename_list[k]
					(Total_median, Total_mean,  Total_x) =load_mass_vs_age_from_mat_file(filename,load_type='median')
					plt.title(Name_list[k])
					print Total_median.shape, Total_x.shape
					plt.plot(Total_median,Total_x,color=color_vec[k],linewidth=4)
					plt.xlabel(field_name1)
					plt.ylabel(field_name2)
					ax.set_yscale('log')


		if plot_occurance_matrix is True:
			for k in range(2):
				ax = fig.add_subplot(1,1,1)
				filename=load_filename_list[k]
				(occurance_x, occurance_y, occurance_matrix) =load_mass_vs_age_from_mat_file(filename,load_type='occurance')
				dp=occurance_x[1]-occurance_x[0]
				max_occurance=np.zeros((len(occurance_y)))
				mode_occurance=np.zeros((len(occurance_y)))
				for j in range(len(occurance_y)):
					max_occurance[j]=np.sum(  occurance_matrix[:,j]*occurance_x /np.sum(occurance_matrix[:,j])  )
				#plt.subplot(2,1,k+1)
				vmax=np.max(occurance_matrix[:,:])
				#cNorm = MidpointNormalize(vmin=0, vmax=5,midpoint=0)
				#cNorm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)
				cNorm = mpl.colors.LogNorm(vmin=10**-5, vmax=0.1)
				#if k==0:
				#       max_occurance0=max_occurance
				if k==1:
					plt.pcolor(occurance_x, occurance_y, (np.squeeze(occurance_matrix[:,:])).transpose(),cmap='jet', norm=cNorm)
					plt.colorbar()
					plt.title(Name_list[k])
				plt.plot(max_occurance,occurance_y,color=color_vec[k],linewidth=2)
				#plt.plot(mode_occurance,occurance_y,color=color_vec[k],linewidth=2)
				#title=Name_list[k]
				plt.plot(max_occurance,occurance_y,color=color_vec[k],linewidth=1, label=Name_list[k])
				#plt.plot(mode_occurance,occurance_y,color=color_vec[k],linewidth=1, label=Name_list[k])
				ax.set_yscale('log')

		

	plt.legend(loc='upper right', frameon=True)

	#plt.legend(loc='lower right', frameon=True,prop={'size':12})
	#plt.legend(loc='lower right', frameon=True,prop={'size':12})
	#plt.legend(loc='upper right', frameon=True)
	#fig = matplotlib.pyplot.gcf()
	fig.set_size_inches(9,4.5)
	print 'field1 = ', field_name1, ', field2 = ' , field_name2
	plt.show()


	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

