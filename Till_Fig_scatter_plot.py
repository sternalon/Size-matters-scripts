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




def expected_values(x,y,number_of_bins,min_bin_value=None):
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
	#x1 = np.logspace(np.log10(x1_min), np.log10(x1_max), num=number_of_bins, endpoint=True)
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

def save_mat_file(mean_matrix, std_matrix, x_matrix, Total_mean, Total_std, Total_x, start_year,end_year,Delta_name,field_name1,field_name2,\
		sigma_matrix, scale_matrix, loc_matrix, mu_matrix, Total_sigma, Total_loc, Total_scale, Total_mu,Total_median,constraints):
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
	
	mat_filename='processed_data/' + mat_filename +'.mat'

	sc.savemat(mat_filename, {'mean_matrix':mean_matrix,'std_matrix':std_matrix,'x_matrix':x_matrix,\
			'Total_mean':Total_mean,'Total_std':Total_std,'Total_x':Total_x ,\
			'sigma_matrix':sigma_matrix, 'scale_matrix':scale_matrix, 'loc_matrix':loc_matrix, 'mu_matrix':mu_matrix,\
			'Total_sigma':Total_sigma, 'Total_loc': Total_loc, 'Totla_scale':Total_scale, 'Total_mu':Total_mu, 'Total_median':Total_median}) 
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
	plot_average_dist=1
	plot_full_time_series=0
	plot_scatter_plot=True
	plot_mean_fields=False
	plot_mean_fields_again=False
	plot_derivatives_and_fits=False
	find_expected_values=False
	save_the_data=False
	plot_rolling_equation=False
	save_figure=True

	#Parameters and variables
	#start_year=1980
	#end_year0=2010
	end_year0=2000
	Number_of_years=0
	#Number_of_years=0

	#Which fields do you want to compare?
	#field_name1='area'
	#field_name1='age'
	#field_name1='width'
	#field_name1='length'
	field_name1='width'
	#field_name1='thickness'
	#field_name1='LW_max'
	#field_name1='distance_from_calving'
	
	#field_name2='area'
	#field_name2='age'
	#field_name2='thickness'
	#field_name2='length'
	#field_name2='width'
	#field_name2='LW_max'
	#field_name2='mass'

	initial_mass_vec=np.array([ 8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11])
	#Include contraints to use to filter the data
	
	constraint2 = {'Constraint_field_name': 'lon', 'lower_bound': -150 , 'upper_bound': -65, 'original_values': True} #Longitude of AB
	constraint_SH = {'Constraint_field_name': 'lat', 'lower_bound': -150 , 'upper_bound': -5, 'original_values': True} #Longitude of AB
	constraint_depth = {'Constraint_field_name': 'depth', 'lower_bound': 500 , 'upper_bound': 8000, 'original_values': False}
	constraint_dist_from_calving = {'Constraint_field_name': 'distance_from_calving', 'lower_bound': 1000 , 'upper_bound': 10000000000000, 'original_values': False}
	#constraint_mass = {'Constraint_field_name': 'mass', 'lower_bound': -10**16 , 'upper_bound': -0.5*(10**8), 'original_values': False, 'subtract_original': True}
	constraint_mass = {'Constraint_field_name': 'mass', 'lower_bound': 3.8*(10**11) , 'upper_bound': 4.0*(10**11), 'original_values': True, 'subtract_original': False}
	constraint_mass9 = {'Constraint_field_name': 'mass', 'lower_bound': 7.3*(10**11) , 'upper_bound': 7.6*(10**11), 'original_values': True, 'subtract_original': False}
	constraint_mass8 = {'Constraint_field_name': 'mass', 'lower_bound': 2.1*(10**11) , 'upper_bound': 2.3*(10**11), 'original_values': True, 'subtract_original': False}
	constraint_mass6 = {'Constraint_field_name': 'mass', 'lower_bound': 7.4*(10**10) , 'upper_bound': 7.6*(10**10), 'original_values': True, 'subtract_original': False}
	constraint_mass3 = {'Constraint_field_name': 'mass', 'lower_bound': 3.2*(10**9) , 'upper_bound': 3.4*(10**9), 'original_values': True, 'subtract_original': False}
	constraint_day = {'Constraint_field_name': 'day', 'lower_bound': 92 , 'upper_bound': 94., 'original_values': True, 'subtract_original': False}
	constraint_vvel = {'Constraint_field_name': 'vvel', 'lower_bound': 0.0001 , 'upper_bound': 100, 'original_values': False, 'subtract_original': False}
	constraint_uvel = {'Constraint_field_name': 'uvel', 'lower_bound': -100 , 'upper_bound': -0.00100, 'original_values': False, 'subtract_original': False}
	constraint_age = {'Constraint_field_name': 'age', 'lower_bound': 0 , 'upper_bound': 5, 'original_values': False} 
	#constraints.append(constraint2)


	#if field_name1=='age' or field_name2=='age':
	#	print 'applying age constraint'
	#	day_constraint = {'Constraint_field_name': 'day', 'lower_bound': 0 , 'upper_bound': 258, 'original_values': True} #Longitude of AB
	constraints=[]	
	#constraints.append(constraint_dist_from_calving)
	#constraints.append(constraint_depth)
	#constraints.append(constraint_mass0)
	#constraints.append(constraint_day)
	#constraints.append(constraint_age)
	#constraints.append(constraint_SH)
	#constraints.append(constraint_vvel)
	#constraints.append(constraint_uvel)
	#	constraints.append(day_constraint)


	mass_scaling=np.array([2000, 200, 50, 20, 10, 5, 2, 1, 1, 1])



	#field_name='area';  x_min=1;x_max=1.e9  # Good for area including all Deltas
	field_name='mass';  x_min=1;x_max=1.e7  # Good for area 
	#field_name='mass' ; #x_min=100;x_max=1.e12  # Good for mass
	#Defining a list of colors
	#color_vec=np.array(['blue', 'red','purple','green', 'coral', 'cyan', 'magenta','orange', 'black', 'grey', 'yellow', 'orchid', 'blue', 'red','purple','green', 'coral', 'cyan' ])
	color_vec=np.array(['red', 'salmon','mistyrose','green', 'coral', 'cyan', 'magenta','orange', 'black', 'grey', 'yellow', 'orchid', 'blue', 'red','purple','green', 'coral', 'cyan' ])
	color_vec=np.array(['darkred','red', 'salmon','mistyrose','green', 'coral', 'cyan', 'magenta','orange', 'black', 'grey', 'yellow', 'orchid', 'blue', 'red','purple','green', 'coral', 'cyan' ])
	root_path='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_'
	
	#Delta_list=np.array(['Delta9', 'rolling_tournadre_WM', 'rolling_Delta9_Burton_new'])
	#Delta_list=np.array(['rolling_Delta1_Burton_new','rolling_Delta5_Burton_new','rolling_Delta9_Burton_new'])
	#Delta_list=np.array(['rolling_tournadre_Burton_breaking_aspect_fixed','rolling_tournadre_Burton_breaking_D_fixed','rolling_tournadre_Burton_new'])
	#Delta_list=np.array(['rolling_tournadre_Rolling_off','tournadre_C2','rolling_tournadre_WM','rolling_tournadre_Burton_new'])
	#Delta_list=np.array(['rolling_tournadre_Rolling_off','rolling_tournadre_WM','rolling_tournadre_Burton_new'])
	Delta_list=np.array(['rolling_tournadre_Rolling_off','rolling_tournadre_WM_corrected','rolling_tournadre_Burton_new'])
	#Delta_list=np.array(['Delta1','Delta9'])
	#Delta_list=np.array(['Delta2','Delta3', 'Delta4','Delta5','Delta6'])
	
	#Name_list=np.array(['No Rolling', 'Martin and Adcroft','Weeks and Mellor','Burton'])
	Name_list=np.array(['No Rolling', 'Weeks and Mellor \n(corrected)','Burton'])
	#mass_string=np.array(['3.3e9 kg','2.2e11 kg'])
	mass_string=np.array(['7.5e10 kg','3.3e9 kg'])

	fig = plt.figure(1)
	#for k in range(13):#for k in np.array([3]):
	for outer_loop in range(1):
		if outer_loop==0:
			field_name2='thickness'
			#field_name1='length'
		if outer_loop==1:
			field_name2='length'
		#for loop in range(1):
		for loop in np.array([1]):
			constraints=[]	
			if loop==0:
				constraints.append(constraint_mass6)
			#if loop==2:
			#	constraints.append(constraint_mass8)
			if loop==1:
				constraints.append(constraint_mass3)

			count1=0
			for k in np.array([0,1,2]):
			#for k in np.array([1]):
				count1=count1+1
				#input_folder=all_paths[k]
				#Delta_name='Delta'+str(k)
				Delta_name=Delta_list[k]
				input_folder=root_path+Delta_name 

				#Make sure that data exists of the years allocated, othersize change it.
				print input_folder
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
					

					if plot_scatter_plot==True:
						#Plotting the distributions
						#ax = fig.add_subplot(3,2,2*(k)+1+outer_loop)
						#ax = fig.add_subplot(3,1,1*(k)+1+outer_loop)
						ax = fig.add_subplot(1,1,1)
						#ax = fig.add_subplot(3,2,k+1)
						#ax = fig.add_subplot(2,1,count1)
						#label=Name_list[k]
						#label='Calving mass =' + mass_string[loop]
						label=Name_list[k]
						#plt.plot(field1,field2,'o',color=color_vec[loop],label=label)
						if year_count==0:
							#plt.plot(field2,field1/field2,'o',color=color_vec[loop],label=label)
							#plt.plot(field2,field1/field2,'o',color=color_vec[k],label=label)
							#plt.plot(field2,field1/field2,color=color_vec[k],label=label,linestyle='circ')
							#plt.scatter(field2, field1/field2, s=20, facecolors='none', edgecolors=color_vec[k],label=label)
							plt.scatter(field2, field1/field2, s=4, facecolors=color_vec[k], edgecolors=color_vec[k],label=label)
						else:
							#plt.plot(field2,field1/field2,'o',color=color_vec[k])
							plt.scatter(field2, field1/field2, s=4, facecolors=color_vec[k], edgecolors=color_vec[k])
							#plt.plot(field2,field1/field2,'o',color=color_vec[loop])
						
						ind=np.argmax(field2)
						C1=field1[ind]
						C2=np.max(field2)

						print C1
						H=np.linspace(0.001,C2,1000)
						plt.plot(H,C1/H,color='grey',linewidth=3)#,linestyle=':')
						if Name_list[k]=='Burton':
							rho_ice=850.0; rho_sw=1025.0
							B=np.sqrt(6.*(rho_ice/rho_sw)*(1-(rho_ice/rho_sw)))
							#print B
							plt.plot(H,B+(H*0.0),color='blue',linewidth=3,linestyle='--')
							#plt.plot(H,1.0+(H*0.0),color='k',linewidth=3,linestyle=':')
						if Name_list[k][0]=='W':
							rho_ice=850.0; rho_sw=1025.0 ; Delta=6.0  ; alpha=rho_ice/rho_sw
							#B=np.sqrt(0.92*(H**2)-(58.32*H) )
							#plt.plot(H,B/H,color='k',linewidth=3)
							B=sqrt((6.0*alpha*(1.-alpha)*(H**2))-(12*Delta*alpha*H) ) 
							plt.plot(H,B/H,color='b',linewidth=3)
							B=sqrt((6.0*alpha*(1.-alpha)*(H**2))+(12*Delta*alpha*H) ) 
							#plt.plot(H,B/H,color='b',linewidth=3,linestyle='-.')
							plt.plot(H,B/H,color='cornflowerblue',linewidth=3)

						if k+1==3:
							#plt.xlabel(field_name2)
							plt.xlabel('H (m)')
							plt.legend(loc='upper right', frameon=True,prop={'size':12})
						#plt.xlabel(field_name2)
						plt.ylabel('W / H')
						#plt.ylabel(field_name1+ '/'+ field_name2 )
						#plt.ylim([0., 1.5]) #For mass
						plt.ylim([0., 2.2]) #For mass
						plt.xlim([0, 140]) #For mass
						#plt.xlim([0, 135]) #For mass
						#ax.set_yscale('log')
						#ax.set_xscale('log')
						#plt.title(Name_list[k])
						#plt.title('Mass Class 3')
						plt.plot(C2,C1/C2,color='grey', marker='*',markersize=18)
						print C2, C1/C2
			


		if plot_rolling_equation==True:
			D=np.linspace(0.1,300,100)
			L=sqrt(0.92*(D**2)+(58.32*D))
			plt.plot(L,D,linestyle='--')
			#plt.plot(L,L,linestyle=':')
			plt.plot(L,0.81*L,linestyle=':')
			#plt.plot(L,L/1.2,linestyle=':')
			#plt.plot(L,L-(,linestyle=':')
			
			
		if save_the_data==True:
			save_mat_file(mean_matrix, std_matrix, x_matrix, Total_mean, Total_std, Total_x, start_year,end_year,Delta_name,field_name1,field_name2,\
					sigma_matrix, scale_matrix, loc_matrix, mu_matrix, Total_sigma, Total_loc, Total_scale, Total_mu,Total_median,constraints)

	#plt.legend(loc='lower right', frameon=True,prop={'size':12})
	#plt.legend(loc='upper right', frameon=True)
	#fig = matplotlib.pyplot.gcf()
	fig.set_size_inches(9,4.5)
	#fig.set_size_inches(12.0, 10.5,forward=True)
	if save_figure==True:
		figure_name='Till_paper_figures/Till_Fig_scatter_plot2.png'
		plt.savefig(figure_name,dpi=300,bbox_inches='tight')
		print 'Saving file : ', figure_name

	plt.show()


	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

