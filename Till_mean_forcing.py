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




def expected_values(x,y,number_of_bins,min_bin_value=None,use_log_axis=False):
	#finding expected value of x in bins of y
	plot_distribution=False

	#age is y   (field1)
	#mass is x  (field2)

	if min_bin_value is not None:	
		x1_min=max(10**5,min(x)) ; x1_max=max(x)   #Use with mass vs age
	else:
		x1_min=min(x) ; x1_max=max(x)
	iceberg_age=0
	iceberg_mass=0
	#number_of_bins=300
	if use_log_axis is True:
		x1 = np.logspace(np.log10(x1_min), np.log10(x1_max), num=number_of_bins, endpoint=True)
	else:
		x1 = np.linspace(x1_min, x1_max, num=number_of_bins, endpoint=True)
	
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

		#if np.isnan(y_median[count]):
		#	print 'STOP', valid
		#	print 'x1', x1, x
		#	print values,x1[count],x1[count+1]
		#	halt

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

def filter_curve(age, mass,k,mass_ind):
	mass=mass[np.where(np.isnan(age)==0)]
	age=age[np.where(np.isnan(age)==0)]
	tol=0.3
	#tol=np.max(age)/1.
	#print tol
	filtered_age=age
	filtered_mass=mass
	#for i in range(len(age)-1):
	#	if abs(age[i] - age[i-1])>tol:#   and  abs(age[i] - age[i-1])>tol :
	#		filtered_age[i]=filtered_age[i-1]
	#		filtered_mass[i]=filtered_mass[i+1]
	#for i in range(len(age)-1,30,-1):

	N=0
	#if k==0 and mass_ind>2  and mass_ind<7:
	#	N=13
	filtered_age[0:N]=0.0
	filtered_mass[0:N]=0.0


	return [filtered_age , filtered_mass]

def fit_curve(age,mass): #x is age, y in mass
	
	#Parameters used in fitting
	max_age=14
	number_of_points=40
	N=13 #polynomial degree

        #Filtering the data
	mass=mass[np.where((age<max_age) *(~np.isnan(age)))]
	age=age[np.where((age<max_age) *(~np.isnan(age)))]
	mass=np.log(mass)
	
	if len(mass)>0:
		#Creating new age vector
		new_age=np.linspace(np.min(age),np.max(age),number_of_points)


		#fit_with_interp=False
		#if fit_with_interp==True:
		#	y2=np.linspace(0.05,max(y1),30)
		#	set_interp = interp1d(y1, x1, kind='linear')
		#	x2 = set_interp(y2)
		#	x_new=x2 ; y_new=y2
		#else:
		#	x_new=x1 ; y_new=y1

		#Fitting mass=p(age)
		p=np.polyfit(age,mass,N)
		#p=np.polyfit(y_new,x_new,N)
		new_mass=np.polyval(p,new_age)
	else:
		new_mass=0; new_age =0;

	new_mass=np.exp(new_mass)


	return [new_age,new_mass]

def load_data_from_mat_file(filename, field_name):
	print 'loading ' , filename , ' ...'
	mean_field_name='mean_'+field_name
	Total_mean_field_name='Total_mean_'+field_name

        mat_contents=sc.loadmat(filename)
	print 'Keys so far are: ' , sc.loadmat(filename).keys()

	field=mat_contents[field_name][:]
	mean_field=mat_contents[mean_field_name][:]
	Total_mean_field=mat_contents[Total_mean_field_name][:]
        
        return [field, mean_field, Total_mean_field]

def save_mat_file(Delta_name,field1_matrix,field_name1,constraints,mean_field1,Total_mean_field1,start_year, end_year):
	mat_filename= Delta_name + '_Mean_forcing_'+field_name1 + '_' + str(start_year) +'_to_' + str(end_year) 
	
	if constraints!=None:
		if len(constraints)>0:
			mat_filename=mat_filename+'_constraints'
			for k in range(len(constraints)):
				current_constraint=constraints[k]
				constraint_name=current_constraint['Constraint_field_name']
				lower_bound=current_constraint['lower_bound']
				upper_bound=current_constraint['upper_bound']
				mat_filename=mat_filename+ '_'+constraint_name + '_' + str(abs(lower_bound)) +'_to_'+ str(abs(upper_bound))

	mat_filename='processed_data/' + mat_filename +'.mat'
	
	#if os.path.isfile(mat_filename):
	#	with open(mat_filename,'ab') as f:
	#		sc.savemat(f, {field_name1:field1_matrix, 'mean_'+field_name1:mean_field1, 'Total_mean_'+field_name1:Total_mean_field1}) #append
	#		print 'Appending to file: ' + mat_filename + ', saved.'
	#else:
	sc.savemat(mat_filename, {field_name1:field1_matrix, 'mean_'+field_name1:mean_field1, 'Total_mean_'+field_name1:Total_mean_field1})
	print 'Creating file: ' + mat_filename + ', saved.'
	print 'Keys so far are: ' , sc.loadmat(mat_filename).keys() 
	print 'File saving complete'

###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
def main():
	#Clear screen
	#os.system('clear')

	#Defining possible paths
	all_paths=define_paths_array()

	#Flags
	save_the_data=True
	load_the_data=True

	#Parameters and variables
	#start_year=1980
	#end_year0=2010
	end_year0=2000
	Number_of_years=50
	#Number_of_years=0

	#Which fields do you want to compare?
	#field_name1='sst'
	#field_name1='ua'
	#field_name1='va'
	#field_name1='uo'
	#field_name1='vo'
	#field_name1='ui'
	#field_name1='vi'
	#field_name1='cn'
	#field_name1='speed_o'
	#field_name1='speed_i'
	#field_name1='speed_a'
	#field_name1='Mb'
	#field_name1='Me'
	field_name1='Mv'
	

	initial_mass=np.array([ 8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11])
	#Include contraints to use to filter the data
	constraint_SH = {'Constraint_field_name': 'lat', 'lower_bound': -150 , 'upper_bound': -5, 'original_values': True} 
	constraint_NH = {'Constraint_field_name': 'lat', 'lower_bound': 5 , 'upper_bound': 150, 'original_values': True} 
	

	#if field_name1=='age' or field_name2=='age':
	#	print 'applying age constraint'
	#	day_constraint = {'Constraint_field_name': 'day', 'lower_bound': 0 , 'upper_bound': 258, 'original_values': True} #Longitude of AB
	constraints=[]	
	#constraints.append(constraint_NH)
	constraints.append(constraint_SH)


	rho_ice=850.0
	rho_sw=1025.0

	mass_ind_list=np.array([0,1,2,3,4,5,6,7,8,9])
	#mass_ind_list=np.array([9])


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
	Delta_list=np.array(['rolling_tournadre_Burton_new','rolling_tournadre_Rolling_off'])
	#Delta_list=np.array(['rolling_tournadre_Burton_new/trajectories/','rolling_tournadre_Rolling_off/trajectories/'])
	#Delta_list=np.array(['Delta1','Delta9'])
	#Delta_list=np.array(['Delta2','Delta3', 'Delta4','Delta5','Delta6'])
	
	#Name_list=np.array(['Martin and Adcroft','Weeks and Mellor','Burton'])
	Name_list=np.array(['Rolling', 'No Rolling'])

	filename='processed_data/rolling_tournadre_Burton_new_Mean_forcing_'+ field_name1+'_1950_to_2000.mat'


	#load_filename_list=np.array(['processed_data/rolling_tournadre_Rolling_off_age_vs_mass_1950_to_2000_expected_values.mat',\
	#	'processed_data/rolling_tournadre_Burton_new_age_vs_mass_1950_to_2000_expected_values.mat'])
	load_filename_list=np.array(['processed_data/rolling_tournadre_Rolling_off_age_vs_mass_1950_to_2000_expected_values_logaxis.mat',\
		'processed_data/rolling_tournadre_Burton_new_age_vs_mass_1950_to_2000_expected_values_logaxis.mat'])

	fig = plt.figure(1)
	ax = fig.add_subplot(1,1,1)
	if load_the_data==True:
		count1=0
		for k in np.array([1]):
		#for k in np.array([1]):
			count1=count1+1
			#input_folder=all_paths[k]
			#Delta_name='Delta'+str(k)
			Delta_name=Delta_list[k]
			input_folder=root_path+Delta_name +'/trajectories/'

			#Make sure that data exists of the years allocated, othersize change it.
			(min_year, max_year)=find_max_min_years(input_folder,'.iceberg_trajectories.nc')
			end_year=min(max_year,end_year0)
			start_year=max(end_year-Number_of_years,min_year)

			number_of_bins=500
			Number_of_classes=10

			#Initializing matricies
			field1_matrix=np.zeros((Number_of_classes,end_year-start_year+1))

			year_count=-1
			for year in range(start_year,end_year+1):
				year_count=year_count+1
				#year=1945
				filename ='/' + str(year) + '0101.iceberg_trajectories.nc'
				input_file=input_folder + filename
				print input_file
				
				field1_full = get_valid_data_from_traj_file(input_file,field_name1,subtract_orig=False)

				print 'Lendth of field: ' , len(field1_full)
				#Getting the distributions from the data


				#Handling all the different mass classes
				for mass_ind in mass_ind_list:
					mass=initial_mass[mass_ind]
					constraint_m = {'Constraint_field_name': 'mass', 'lower_bound': mass-100 , 'upper_bound': mass+100, 'original_values': True}
					mass_constraints=[]
					mass_constraints=constraints
					mass_constraints.append(constraint_m)
					#Handeling constraints
					field1=0. ; field2=0.
					field1=add_constraint_to_data(input_file,field1_full,mass_constraints)
					print 'Lendth of field after constraints using mass0= ',   mass ,' is ' , len(field1)
					mass_constraints.pop()

					field1_matrix[mass_ind,year_count]=np.mean(field1)
		mean_field1=np.mean(field1_matrix,axis=1)
		Total_mean_field1=np.mean(mean_field1)
		if save_the_data is True:
			save_mat_file(Delta_name,field1_matrix,field_name1,constraints,mean_field1,Total_mean_field1,start_year, end_year)
					
	else:
	       (field1_matrix, mean_field1, Total_mean_field1)=load_data_from_mat_file(filename,field_name1)	


	#print field1_matrix[:,:]				
	print mean_field1
	print Total_mean_field1

	#plt.legend(loc='upper right', frameon=True,prop={'size':12})
	#plt.legend(loc='upper right', frameon=True)
	#fig = matplotlib.pyplot.gcf()
	#fig.set_size_inches(9,4.5)
	#plt.show()


	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

