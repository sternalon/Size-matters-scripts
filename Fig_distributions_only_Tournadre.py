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


def load_data_from_mat_file(filename):
        mat_contents=sc.loadmat(filename)
	lat=mat_contents['lat'][0,:]
	lon=mat_contents['lon'][0,:]
	data=mat_contents['data'][:,:]
	field=mat_contents['field']
	title=mat_contents['title']
	p_values=mat_contents['p_values']
	if p_values==0:
		p_values=None
	return [lat, lon, data,p_values,field ]


def load_distribution_from_mat_file(filename):
        mat_contents=sc.loadmat(filename)
	#print mat_contents

	prob_dist_matrix=mat_contents['prob_dist_matrix']
	x_matrix=mat_contents['x_matrix']
	split_prob_matrix=mat_contents['split_prob_matrix']
	

	#if split_prob_matrix==0:
	#	split_prob_matrix=None
	return [prob_dist_matrix,x_matrix, split_prob_matrix]


def produce_distribution_panel(filename,nrows,ncols,position,field_name,Delta_list,Area_list,color_vec,letter_labels):

	print 'Plotting ' + field_name + ' fields for ' + Delta_list[0] 
	(prob_dist_matrix,x_matrix, split_prob_matrix)=load_distribution_from_mat_file(filename)
	
	prob_dist=np.squeeze(np.mean(prob_dist_matrix, axis=0))
	x=np.squeeze(np.mean(x_matrix, axis=0))
	
	ax=subplot(nrows,ncols,position)
	label= 'All Areas'
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

	plt.scatter(Gladstone_A, h_Tournadre, marker='x',linewidth=4,color='red',label='Calving\ndistribution', zorder=10)
	#plt.scatter(Gladstone_A2, h_Tournadre, marker='x',linewidth=4,color='magenta',label='Calving distribution', zorder=11)
	plt.plot(x,y, linewidth=1.5, linestyle=":",color='black',label='$p=A^{-1.5}$')


	text(1,1,letter_labels[position-1], ha='right', va='bottom',transform=ax.transAxes,fontsize=15)

	#ax = fig.add_subplot(1,2,2)
	#plt.ylim([1e-7, 1e-2]) #For area
	#plt.ylim([1e-10, 1e-2]) #For area, with constraints
	plt.ylim([1e-9, 1e-2]) #For area, with constraints
	plt.ylim([10**(-8.5), 10**(-2.8)]) #For area, with constraints
	plt.xlim([5, 1e7]) #For area, with constraints
	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.xlabel('Area (m$^2$)')
	if position==1:
		plt.ylabel('Probabilty')
	if position==2:
		ax.set_yticks([])

	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +  ax.get_xticklabels() + ax.get_yticklabels()):
		item.set_fontsize(20)

##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################

def main():

	parser = argparse.ArgumentParser()
   	args = parser.parse_args()
	
	
	#Flags
	save_figure=True
	significant_vals_only=False
	plot_melt_fields=False
	plot_split_Tournadre_dist=True
	plot_split_Tournadre_dist_depth=True
	
	plot_all_panels_without_data=False   #False makes the real figure

	#Parameters;
	pole='south'
	title=None
	cscale=None
	Delta_list=np.array(['tournadre'])
	#Area_list=np.array(['Area$_{1}$=0.0026km$^2$', 'Area$_{2}$=0.0072km$^2$','Area$_{3}$=0.029km$^2$',\
	#	'Area$_{4}$=0.12km$^2$','Area$_{5}$=0.18km$^2$', 'Area$_{6}$=0.35km$^2$','Area$_{7}$=0.56km$^2$', 'Area$_{8}$=1.0km$^2$','Area$_{9}$=1.8km$^2$', 'Area$_{10}$=3.5km$^2$'])
	Area_list=np.array(['0.0026km$^2$', '0.0072km$^2$','0.029km$^2$',\
		'0.12km$^2$','0.18km$^2$', '0.35km$^2$','0.56km$^2$', '1.0km$^2$','1.8km$^2$', '3.5km$^2$'])
	Title_list=np.array(['TOURNADRE distribution'])
	letter_labels=np.array(['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'])
	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta',  'orange', 'coral', 'yellow', 'orchid',   'orange', 'coral', 'yellow', 'orchid' ])
	y_title_list=np.array(['Iceberg Melt'])
	colorbar_unit_list=np.array(['melt (kg/m$^2$)'])	


	nrows=1
	ncols=2
	
	datapath='/home/Alon.Stern/Iceberg_Project/iceberg_scripts/python_scripts/size_matters_paper/processed_data/'
        
	#Setting up the figure
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	ax=subplot(nrows,ncols,2)
	fig.delaxes(ax)
	
	
	#Plotting Trajectory Fields

	#Plotting Melt Fields
	if plot_melt_fields==True:
		cscale=None
		row_num=1
		difference_on=0
		field_name='melt'
		file_extension='_melt_1959_to_2019_all_south_xy.mat'
		title=None
		filename=datapath + Delta_list[0] + file_extension
		(lat ,lon, data, p_values,field)=load_data_from_mat_file(filename)
		print 'Plotting ' + field_name + ' fields for ' + Delta_list[row_num-1] 
		ax=subplot(nrows,ncols,row_num)
		(cscale,datamap)=plot_polar_field(lat,lon,data,pole,difference_on,title,p_values,cscale,field,colorbar_on=False,return_data_map=True)
		#Creating colorbar
		cbar_ax = fig.add_axes([0.52, 0.62, 0.04, 0.3])
		cbar=fig.colorbar(datamap, cax=cbar_ax)
		cbar.set_label(colorbar_unit_list[row_num-1], rotation=90)

	if plot_split_Tournadre_dist==True:
		position=1
		field_name='size distribution no constrains'
		label= Delta_list[0]
		file_extension='_size_distribution_1959_to_2019.mat'
		filename=datapath + Delta_list[0] + file_extension
		produce_distribution_panel(filename,nrows,ncols,position,field_name,Delta_list,Area_list,color_vec,letter_labels)

	if plot_split_Tournadre_dist_depth==True:
		position=2
		field_name='size distribution no constrains'
		file_extension='_size_distribution_1959_to_2019_constraints_lat_10_to_70.mat'
		file_extension='_size_distribution_1959_to_2019_constraints_depth_8000_to_1000.mat'
		filename=datapath + Delta_list[0] + file_extension
		produce_distribution_panel(filename,nrows,ncols,position,field_name,Delta_list,Area_list,color_vec,letter_labels)
		fig.subplots_adjust(right=0.65)
		plt.legend(loc=(1.1, 0.0), frameon=True,prop={'size':20})

	subplots_adjust(left=None, bottom=None, right=0.8, top=None, wspace=None, hspace=None)
	#fig.tight_layout()
	fig.set_size_inches(20.5, 10.5,forward=True)
	if save_figure==True:
		plt.savefig('paper_figures/Fig_distributions_only_Tournadre.png',dpi=300,bbox_inches='tight')
	plt.show()


if __name__ == '__main__':
	main()
	#sys.exit(main())

