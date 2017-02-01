#!/usr/bin/env python

import sys
import argparse
import netCDF4 as nc
import numpy as np
import scipy.io as sc
from scipy import stats
import matplotlib
import matplotlib as mpl
from pylab import *
from matplotlib import ticker
#matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from sea_ice_concentrations import plot_polar_field
#from sea_ice_concentrations import define_longitudes
from sea_ice_concentrations import create_output_file_and_save


def load_distribution_from_mat_file(filename,load_split_matrix):
        mat_contents=sc.loadmat(filename)
	#print mat_contents

	prob_dist_matrix=mat_contents['prob_dist_matrix']
	x_matrix=mat_contents['x_matrix']

	if load_split_matrix==True:
		split_prob_matrix=mat_contents['split_prob_matrix']
		return [prob_dist_matrix,x_matrix, split_prob_matrix]
	else:
		return [prob_dist_matrix,x_matrix]
	

	#if split_prob_matrix==0:
	#	split_prob_matrix=None


def load_mass_vs_age_from_mat_file(filename,load_type='lognormal'):

        mat_contents=sc.loadmat(filename)
	Total_x=mat_contents['Total_x'][0]
	
	if load_type=='lognormal':
		#sigma_matrix=mat_contents['sigma_matrix']
		#scale_matrix=mat_contents['scale_matrix']
		#return [scale_matrix, sigma_matrix]
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


	

def produce_distribution_panel(filename,nrows,ncols,position,field_name,Delta_list,Area_list,color_vec,letter_labels):


	print 'Plotting ' + field_name + ' fields for ' + Delta_list[0] 
	(prob_dist_matrix,x_matrix, split_prob_matrix)=load_distribution_from_mat_file(filename,load_split_matrix=True)
	
	prob_dist=np.squeeze(np.mean(prob_dist_matrix, axis=0))
	x=np.squeeze(np.mean(x_matrix, axis=0))
	
	ax=subplot(nrows,ncols,position)
	label= 'All Iceberg Areas'
	plt.plot(x,prob_dist, color='black', linewidth=4.,label=label)

	#Plotting the split up distributions
	M=split_prob_matrix.shape
	for mass_num in range(M[2]):
		#label= 'D'+str(mass_num+1)
		label=Area_list[mass_num]
		prob_dist=np.squeeze(np.mean(split_prob_matrix[:,:,mass_num], axis=0))
		x=np.squeeze(np.mean(x_matrix, axis=0))
		plt.plot(x,prob_dist, color=color_vec[mass_num], linewidth=2.,label=label)



	Gladstone_A=np.array([0.0026, 0.0072, 0.0292, 0.1210, 0.1788, 0.3529, 0.5647, 1.0353, 1.8353, 3.4824])  #In km^2
	Gladstone_A=Gladstone_A*(10**(6))  #Convert to m^2
	Tournadre_dist=np.array([0.3407, 0.2701, 0.2051, 0.1034, 0.0175, 0.0235,0.0122,0.0120, 0.0085,0.0069])
	dx_A=np.zeros([10,1]) 
	Gladstone_A2=np.zeros([10,1]) 
	h_Tournadre=np.zeros([10,1]) 
	diff_G=diff(Gladstone_A)
	for k in range(10):
		if k==0:
			dx_A[0]=Gladstone_A[0]
			Gladstone_A2[0]=Gladstone_A[0]
		else:
			dx_A[k]=diff_G[k-1]
			Gladstone_A2[k]=Gladstone_A[k]-(dx_A[k]/2)

		h_Tournadre[k]=Tournadre_dist[k]/dx_A[k]
	y=x**(-1.5)
	#total=0.5*(Gladstone_A[0]**(-0.5)-Gladstone_A[9]**(-0.5))
	total=0.0577 
	y=y/total

	plt.scatter(Gladstone_A, h_Tournadre, marker='x',linewidth=4,color='red',label='Calving distribution', zorder=10)
	#plt.scatter(Gladstone_A2, h_Tournadre, marker='x',linewidth=4,color='magenta',label='Calving distribution', zorder=11)
	plt.plot(x,y, linewidth=1.5, linestyle=":",color='black',label='$y=x^{-1.5}$')




	#ax = fig.add_subplot(1,2,2)
	text(1,1,letter_labels[position-1], ha='right', va='bottom',transform=ax.transAxes)
	plt.ylim([1e-7, 1e-2]) #For area
	plt.ylim([1e-10, 1e-2]) #For area, with constraints
	plt.ylim([1e-9, 1e-2]) #For area, with constraints
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.xlabel('Area (m$^2$)')
	plt.ylabel('Probabilty')


def plot_power_laws_on_figure(plot_1_2_line, plot_2_3_line, plot_3_2_line):
	if plot_1_2_line==True:
		y=x**(-0.5)
		total=1e6
		y=y/total
		plt.plot(x,y, linewidth=1.5, linestyle=":",color='black',label='$y=x^{-1.5}$')

	if plot_2_3_line==True:
		y=x**(-0.667)
		total=1e8
		y=y/total
		plt.plot(x,y, linewidth=1.5, linestyle=":",color='black',label='$y=x^{-1.5}$')

	if plot_3_2_line==True:
		y=x**(-1.5)
		total=0.0577
		y=y/total
		plt.plot(x,y, linewidth=1.5, linestyle=":",color='black',label='$y=x^{-1.5}$')


def produce_delta_distribution_panel(datapath,nrows,ncols,position,field_name,Delta_list,Area_list,color_vec,letter_labels,mass_vs_area):
	
	year_range_list={'Delta1':'1959_to_2019','Delta2':'1939_to_1999','Delta3':'1959_to_2019','Delta4':'1959_to_2019','Delta5':'1959_to_2019',\
			'Delta6':'1959_to_2019','Delta7':'1930_to_1955','Delta8':'1912_to_1927','Delta9':'1959_to_2019','Delta10':'1923_to_1948'}
	#		'Delta6':'1959_to_2019','Delta7':'1900_to_1955','Delta8':'1900_to_1927','Delta9':'1959_to_2019','Delta10':'1900_to_1948'}

	#Plotting no constraints data
	#for num in np.array([0,2,3,4,5,8]):
	for num in np.array([0,1,2,3,4,5,6,8,9]):  #Delta8 missing
	#for num in range(10):
		print 'Plotting ' + field_name + ' fields for ' + Delta_list[num] 
		file_extension='_size_distribution_' + year_range_list[Delta_list[num]]+'.mat'
		#file_extension='_size_distribution_' + year_range_list[Delta_list[num]]+ '_area_constraints_distance_from_calving_10000000000000_to_50000.mat'
		#file_extension='_size_distribution_' + year_range_list[Delta_list[num]]+ '_mass_constraints_distance_from_calving_10000000000000_to_30000.mat'
		if mass_vs_area=='mass':
			file_extension='_size_distribution_' + year_range_list[Delta_list[num]]+ '_mass.mat'

		filename=datapath + Delta_list[num] + file_extension

		(prob_dist_matrix,x_matrix)=load_distribution_from_mat_file(filename,load_split_matrix=False)
		prob_dist=np.squeeze(np.mean(prob_dist_matrix, axis=0))
		x=np.squeeze(np.mean(x_matrix, axis=0))
		
		ax=subplot(nrows,ncols,position)
		label= 'All Iceberg Areas'
		plt.plot(x,prob_dist, color=color_vec[num], linewidth=4.)
	
	if mass_vs_area=='area':
		#Plotting with depth constraints data
		for num in range(10):
			print 'Plotting ' + field_name + ' fields for ' + Delta_list[num] 
			#file_extension='_size_distribution_1959_to_2019_constraints_depth_8000_to_1000.mat'
			file_extension='_size_distribution_' + year_range_list[Delta_list[num]]+'_constraints_depth_8000_to_1000.mat'
			filename=datapath + Delta_list[num] + file_extension

			(prob_dist_matrix,x_matrix)=load_distribution_from_mat_file(filename,load_split_matrix=False)
			prob_dist=np.squeeze(np.mean(prob_dist_matrix, axis=0))
			x=np.squeeze(np.mean(x_matrix, axis=0))
			
			ax=subplot(nrows,ncols,position)
			label= 'All Iceberg Areas'
			plt.plot(x,prob_dist, color=color_vec[num], linewidth=4.,linestyle=":")

	plot_power_laws_on_figure(plot_1_2_line=False,plot_2_3_line=False,plot_3_2_line=False)

	plot_tipping_point=True
	if plot_tipping_point==True:
		Area_list=np.array([0.0026, 0.0072,0.029,0.12,0.18,0.35,0.56,1.0,1.8,3.5])
		Depth_list=np.array([40,67,133,175,250,250,250,250,250,250])
		tipping_L=sqrt(0.92*Depth_list**2 +58.32*Depth_list)
		tipping_area=(tipping_L*tipping_L/1.5)
		tipping_mass=(tipping_L*tipping_L/1.5)*Depth_list*900
		print tipping_L
		print tipping_area
		for num in range(5):
			if mass_vs_area=='area':
				plt.plot(np.array([tipping_area[num],tipping_area[num]]),np.array([10e-10,10e-2]), color=color_vec[num], linewidth=0.5, linestyle="--")
			if mass_vs_area=='mass':
				plt.plot(np.array([tipping_mass[num],tipping_mass[num]]),np.array([10e-14,10e-2]), color=color_vec[num], linewidth=0.5, linestyle="--")


	#ax = fig.add_subplot(1,2,2)
	text(1,1,letter_labels[position-1], ha='right', va='bottom',transform=ax.transAxes)
	
	if mass_vs_area=='area':
		plt.ylim([1e-9, 1e-2]) #For area, with constraints
		plt.xlabel('Area (m$^2$)')
	if mass_vs_area=='mass':
		plt.ylim([1e-14, 1e-4]) #For mass
		plt.xlabel('Mass (kg)')

	ax.set_xscale('log')
	ax.set_yscale('log')
	#plt.xlabel('Area (m$^2$)')
	plt.ylabel('Probabilty')





def produce_mass_vs_age_panel(datapath,nrows,ncols,position,field_name,Delta_list,Area_list,color_vec,letter_labels,fig):
	
	year_range_list={'Delta1':'1959_to_2019','Delta2':'1939_to_1999','Delta3':'1959_to_2019','Delta4':'1959_to_2019','Delta5':'1959_to_2019',\
			'Delta6':'1959_to_2019','Delta7':'1930_to_1955','Delta8':'1912_to_1927','Delta9':'1959_to_2019','Delta10':'1923_to_1948'}
	#		'Delta6':'1959_to_2019','Delta7':'1900_to_1955','Delta8':'1900_to_1927','Delta9':'1959_to_2019','Delta10':'1900_to_1948'}

	#Plotting no constraints data
	#for num in np.array([9]):
	Experiment_numbers=np.array([0,1,2,3,4,5,8])
	# pr	xperiment_numbers=np.array([0,1,2,3,4,5,8])
	life_times=np.zeros([len(Experiment_numbers),1])
	Area_value=np.zeros([len(Experiment_numbers),1])
	count=-1
	for num in Experiment_numbers:  #Delta8 missing
		count=count+1
	#for num in range(10):
		print 'Plotting ' + field_name + ' fields for ' + Delta_list[num] 
		#file_extension='_size_distribution_1959_to_2019.mat'
		file_extension='_age_vs_mass_' + year_range_list[Delta_list[num]]+'_expected_values.mat'
		filename=datapath + Delta_list[num] + file_extension

		(mean_age,std_age,mass_bins)=load_mass_vs_age_from_mat_file(filename,load_type='normal')
		#(mean_age,Total_std,mass_bins)=load_mass_vs_age_from_mat_file(filename,load_type='lognormal')
		mass_bins=mass_bins[3:len(mass_bins)-1]
		mean_age=mean_age[3:len(mean_age)-1]
		std_age=std_age[3:len(std_age)-1]
		
		#Normalizing by the initial mass
		mass_bins=mass_bins/max(mass_bins)

		#mean_age=mean_age/max(mean_age)
		life_times[count]=np.max(mean_age)
		
		#for k in range(len(mean_age)-1):
		#	Area_value[count]=Area_value[count]+(mass_bins[k])*(-mean_age[k+1]+mean_age[k])
		


		ax=subplot(nrows,ncols,position)
		label= 'All Iceberg Areas'
		plt.plot(mean_age, mass_bins, color=color_vec[num], linewidth=4.)
		#plt.plot(life_times[count]-mean_age, mass_bins, color=color_vec[num], linewidth=4.)
		
		#plt.plot(mean_age, mass_bins, 'o',color=color_vec[num])
		#plt.plot(std_age, mass_bins, color=color_vec[num], linewidth=4.)
		#plt.plot(std_age, mass_bins, 'o',color=color_vec[num])
	
	text(1,1,letter_labels[position-1], ha='right', va='bottom',transform=ax.transAxes)
	#print 1./Area_value
	#ax = fig.add_subplot(1,2,2)
	#plt.ylim([1e-7, 1e-2]) #For area
	#plt.ylim([1e-10, 1e-2]) #For area, with constraints
	#plt.ylim([1e-10, 1e1]) #For area, with constraints
	#ax.set_xscale('log')
	#ax.set_yscale('log')
	plt.xlabel('Age (years)')
	plt.ylabel('Fraction of initial mass')

		
	#MI and FC from time_series figure script:
	MI=np.array([ 1451.60485153,  1968.97211227,  2852.14401112,  4029.6269745,  4606.9498816, 5205.2215391, 6849.74128436])
	FC=np.array([2502.29828072,  2519.82106124,  2543.77607754,  2493.8283016,  2541.3579632,   2509.7450503,   2541.88664836])
	initial_mass_vec=np.array([ 8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11])
	Good_masses_vec=np.array([ 8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 3.9e11])
	print 'life times:', life_times
	print 'MI/FC:', MI/FC

	plot_Life_vs_M_total=False
	if plot_Life_vs_M_total==True:
		ax=fig.add_axes([0.63,0.72,0.15,0.15])
		##plt.plot((MI/FC),life_times)
		plt.plot(life_times,(MI/FC))
		for num in range(len(life_times)):
			plt.plot(life_times[num],(MI[num]/FC[num]),'o',color=color_vec[Experiment_numbers[num]])
		plt.ylabel('M$_{T}$ / F$_{C}$ (years)')
		plt.xlabel('Average Life time (years)')


	plot_mass_vs_age=True
	if plot_mass_vs_age==True:
		ax=fig.add_axes([0.63,0.72,0.15,0.15])
		##plt.plot((MI/FC),life_times)
		#plt.plot(life_times,(MI/FC))
		plt.plot(Good_masses_vec,life_times)
		for num in range(len(life_times)):
			plt.plot(Good_masses_vec[num],life_times[num],'o',color=color_vec[Experiment_numbers[num]])
			#ax.set_xscale('log')
		#plt.xlabel('(Iceberg Mass)/(Calving Flux) (years)')
		plt.ylabel('Lifetime (years)')
		plt.xlabel('Calving Mass (kg)')



##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################

def main():

	parser = argparse.ArgumentParser()
   	args = parser.parse_args()
	
	
	#Flags
	save_figure=True
	significant_vals_only=False
	plot_split_Tournadre_dist=True
	plot_split_Tournadre_dist_depth=True
	plot_Delta_distributions=True
	plot_mass_vs_age=True
	
	plot_all_panels_without_data=False   #False makes the real figure

	#Parameters;
	pole='south'
	title=None
	cscale=None
	Tournadre_list=np.array(['tournadre'])
	Delta_list=np.array(['Delta1', 'Delta2', 'Delta3', 'Delta4', 'Delta5','Delta6', 'Delta7', 'Delta8', 'Delta9', 'Delta10'])
	letter_labels=np.array(['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'])
	Area_list=np.array(['Area$_{1}$=0.0026km$^2$', 'Area$_{2}$=0.0072km$^2$','Area$_{3}$=0.029km$^2$',\
		'Area$_{4}$=0.12km$^2$','Area$_{5}$=0.18km$^2$', 'Area$_{6}$=0.35km$^2$','Area$_{7}$=0.56km$^2$', 'Area$_{8}$=1.0km$^2$','Area$_{9}$=1.8km$^2$', 'Area$_{10}$=3.5km$^2$'])
	Title_list=np.array(['TOURNADRE distribution'])
	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta',  'orange', 'coral', 'yellow', 'orchid',   'orange', 'coral', 'yellow', 'orchid' ])
	y_title_list=np.array(['Iceberg Melt'])
	colorbar_unit_list=np.array(['melt (kg/m$^2$)'])	


	nrows=2
	ncols=2
	
	datapath='/home/Alon.Stern/Iceberg_Project/iceberg_scripts/python_scripts/size_matters_paper/processed_data/'
        
	#Setting up the figure
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	ax=subplot(nrows,ncols,2)
	fig.delaxes(ax)
	
	if plot_Delta_distributions==True:
		position=1
		field_name='size distribution for Deltas'
		mass_vs_area='mass'
		produce_delta_distribution_panel(datapath,nrows,ncols,position,field_name,Delta_list,Area_list,color_vec,letter_labels,mass_vs_area)
		
	if plot_mass_vs_age==True:
		position=2
		field_name='mass vs time for Deltas'
		produce_mass_vs_age_panel(datapath,nrows,ncols,position,field_name,Delta_list,Area_list,color_vec,letter_labels,fig)
	

	if plot_split_Tournadre_dist==True:
		position=3
		field_name='size distribution no constrains'
		label= Tournadre_list[0]
		file_extension='_size_distribution_1959_to_2019.mat'
		filename=datapath + Tournadre_list[0] + file_extension
		produce_distribution_panel(filename,nrows,ncols,position,field_name,Tournadre_list,Area_list,color_vec,letter_labels)

	if plot_split_Tournadre_dist_depth==True:
		position=4
		field_name='size distribution depth constrains'
		file_extension='_size_distribution_1959_to_2019_constraints_lat_10_to_70.mat'
		file_extension='_size_distribution_1959_to_2019_constraints_depth_8000_to_1000.mat'
		filename=datapath + Tournadre_list[0] + file_extension
		produce_distribution_panel(filename,nrows,ncols,position,field_name,Tournadre_list,Area_list,color_vec,letter_labels)
		fig.subplots_adjust(right=0.7)
		plt.legend(loc=(1.1, 0.0), frameon=True,prop={'size':12})
		

	subplots_adjust(left=None, bottom=None, right=0.8, top=None, wspace=None, hspace=None)
	#fig.tight_layout()
	fig.set_size_inches(20.5, 20.5,forward=True)
	if save_figure==True:
		plt.savefig('paper_figures/Fig_distributions_Tournadre_and_Deltas_with_age_vs_mass.png')
	plt.show()


if __name__ == '__main__':
	main()
	#sys.exit(main())

