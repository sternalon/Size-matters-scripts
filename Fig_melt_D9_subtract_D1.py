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





##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################

def main():

	parser = argparse.ArgumentParser()
   	args = parser.parse_args()
	
	
	#Flags
	save_figure=True
	significant_vals_only=False
	plot_melt_fields=True

	plot_all_panels_without_data=False   #False makes the real figure
	

	#Parameters;
	pole='south'
	title=None
	cscale=None
	Delta_list=np.array(['Delta1', 'Delta9'])
	#Title_list=np.array(['Delta1 (Area=0.0026km^2)', 'Delta6 (Area=0.029km^2)', 'Delta9 (Area=1.8km^2)'])
	#Title_list=np.array(['Area=0.0026km^2', 'Area=1.8km^2'])
	Title_list=np.array(['DELTA1 - DELTA9', 'DELTA9 -DELTA1'])
	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	letter_labels=np.array(['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'])
	y_title_list=np.array(['Iceberg Melt'])
	colorbar_unit_list=np.array(['melt (m/year)'])
	nrows=1
	ncols=2
	
	datapath='/home/Alon.Stern/Iceberg_Project/iceberg_scripts/python_scripts/size_matters_paper/processed_data/'
        
	#Setting up the figure
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	
	
	#Plotting Trajectory Fields

	#Plotting Melt Fields
	if plot_melt_fields==True:
		cscale=None
		row_num=1
		difference_on=0
		field_name='melt'
		file_extension='_melt_1959_to_2019_all_south_xy.mat'

		title=None

		filename1=datapath + Delta_list[0] + file_extension
		filename2=datapath + Delta_list[1] + file_extension
		(lat ,lon, data1, p_values,field)=load_data_from_mat_file(filename1)
		(lat ,lon, data2, p_values,field)=load_data_from_mat_file(filename2)

		#Converting melt rates to m/yr
		data1=data1*(60*60*24*356.25) #Multiply to get units /year rather than /s
		data1=data1/850 # Multiply to get units of m/year
		data2=data2*(60*60*24*356.25) #Multiply to get units /year rather than /s
		data2=data2/850 # Multiply to get units of m/year


		#Panel 1
		print 'Plotting ' + field_name + ' fields for ' + Delta_list[0] + ' minus '+ Delta_list[1]
		ax=subplot(nrows,ncols,1)
		(cscale,datamap)=plot_polar_field(lat,lon,data1-data2,pole,difference_on,title,p_values,cscale,field,colorbar_on=False,return_data_map=True)
		plt.ylabel(y_title_list[row_num-1], fontsize=18, labelpad=25)
		ax.set_title(Title_list[0], fontsize=18,y=1.13)
		text(1,1,letter_labels[0], ha='right', va='bottom',transform=ax.transAxes)

		#Panel 2
		print 'Plotting ' + field_name + ' fields for ' + Delta_list[1] + ' minus '+ Delta_list[0]
		ax=subplot(nrows,ncols,2)
		(cscale,datamap)=plot_polar_field(lat,lon,data2-data1,pole,difference_on,title,p_values,cscale,field,colorbar_on=False,return_data_map=True)
		ax.set_title(Title_list[1], fontsize=18,y=1.13)
		text(1,1,letter_labels[1], ha='right', va='bottom',transform=ax.transAxes)
		
		#Creating colorbar
		fig.subplots_adjust(right=0.8)
		cbar_ax = fig.add_axes([0.85, 0.3, 0.05, 0.4])
		cbar=fig.colorbar(datamap, cax=cbar_ax)
		cbar.set_label(colorbar_unit_list[row_num-1], rotation=90)
		if field!='melt':
			tick_locator = ticker.MaxNLocator(nbins=5)
			cbar.locator = tick_locator
			cbar.update_ticks()


	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
	#fig.set_size_inches(10.5, 20.5,forward=True)
	if save_figure==True:
		plt.savefig('paper_figures/Fig_melt_D9_subtract_D1.png')
	plt.show()


if __name__ == '__main__':
	main()
	#sys.exit(main())

