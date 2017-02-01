#!/usr/bin/env python

import sys
import argparse
import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import scipy.io as sc
from scipy import stats
import matplotlib
import matplotlib as mpl
from pylab import *
#matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import Normalize
from sea_ice_concentrations import plot_polar_field
from sea_ice_concentrations import plot_vertical_section
#from sea_ice_concentrations import define_longitudes
from sea_ice_concentrations import create_output_file_and_save


def load_traj_data_from_mat_file(filename):
        mat_contents=sc.loadmat(filename)
	lat=mat_contents['lat'][0,:]
	lon=mat_contents['lon'][0,:]
	Total_berg_lat=mat_contents['Total_berg_lat'][0,:]
	Total_berg_lon=mat_contents['Total_berg_lon'][0,:]
	field=mat_contents['field']

	return [lat ,lon,Total_berg_lon,Total_berg_lat]


def load_data_from_mat_file(filename):
        mat_contents=sc.loadmat(filename)
	lat=mat_contents['lat'][0,:]
	lon=mat_contents['lon'][0,:]
	data=mat_contents['data'][:,:]
	data=ma.masked_values(data,1.e20)  #Recreating the mask over the data.
	field=mat_contents['field']
	title=mat_contents['title']
	lat_lon_bounds=mat_contents['lat_lon_bounds'][0]
	layer_interface=mat_contents['layer_interface'][0]
	p_values=mat_contents['p_values']
	if p_values==0:
		p_values=None
	return [lat, lon, data,p_values,field,layer_interface,lat_lon_bounds ]

def produce_figure_variable_row(fig,axes,nrows,ncols,file_type_list, title_list,colorbar_unit_list,y_title_list,datapath, file_extension,cscale, row_num, \
		difference_on,pole,field_name,variable_list,axes_direction='xy',letter_labels=None):
			difference_on=0
			title=None
			for k in range(len(variable_list)):
			#for k in np.array([0,1,2]):
				print 'Plotting ' + field_name + ' fields for ' + file_type_list[variable_list[k]]
				ax=subplot(nrows,ncols,(row_num-1)*ncols+(k+1))
				filename=datapath +'Control_' + file_type_list[variable_list[k]] + file_extension
				#print filename
				(lat ,lon, data, p_values,field,layer_interface,lat_lon_bounds)=load_data_from_mat_file(filename)
				if axes_direction=='xy':
					(cscale,datamap)=plot_polar_field(lat,lon,data,pole,difference_on,title,p_values,cscale,field,colorbar_on=False,return_data_map=True)
				else:
					file_type='ocean_month_z'
					datamap=plot_vertical_section(lon,lat,data,difference_on,title,p_values,cscale,field,layer_interface,axes_direction,lat_lon_bounds,file_type,\
							colorbar_on=False,return_data_map=True)
				
				if row_num==1:
					#ax.set_title(title_list[k], fontsize=18, y=1.13)
					ax.set_title(title_list[variable_list[k]], fontsize=18, y=1.13)
				
				#if k==0:
				#	plt.ylabel(y_title_list[row_num-1], fontsize=18, labelpad=25)
				if k>0:
					ax.set_yticks([])
				
				#if k==0 and row_num<2:
					#ax.annotate(y_title_list[row_num-1], xy=(-0.6, 0.5), xycoords='axes fraction', fontsize=21,rotation=90,verticalalignment='center')
					#ax.annotate('Salinity', xy=(-0.45, 0.5), xycoords='axes fraction', fontsize=21,rotation=90,verticalalignment='center')
				#	plt.ylabel(y_title_list[i], fontsize=18, labelpad=25)
				#if k==0 and row_num>1:
					#plt.ylabel('depth (m)', fontsize=14)
					#ax.annotate(y_title_list[row_num-1], xy=(-0.5, 0.5), xycoords='axes fraction', fontsize=21,rotation=90,verticalalignment='center')
				
				#if row_num!=nrows:
				#	ax.set_xticks([])
				#else:
				#	plt.xlabel('latitude (deg)', fontsize=14)
				#	ax.set_xticks([-40 ,-60, -80])

				#if row_num==nrows:
						
				#cbar=fig.colorbar(datamap,shrink=0.5)
				#cbar=fig.colorbar(datamap)
				#Creating colorbar
				#fig.subplots_adjust(bottom=0.15)
				text(1,1,letter_labels[(row_num-1)*ncols+(k+1)-1], ha='right', va='bottom',transform=ax.transAxes,fontsize=15)


				#cbar_ax = fig.add_axes([0.13+(k*0.28),0.07, 0.20, 0.03])
				#cbar=fig.colorbar(datamap, cax=cbar_ax, orientation="horizontal")
				cbar=fig.colorbar(datamap, orientation="horizontal")
				cbar.set_label(colorbar_unit_list[k])
				tick_locator = ticker.MaxNLocator(nbins=5)
				cbar.locator = tick_locator
				cbar.update_ticks()



##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################

def main():

	parser = argparse.ArgumentParser()
   	args = parser.parse_args()
	
	
	#Flags
	save_figure=True
	significant_vals_only=False
	plot_Sea_Surface_fields=True
	plot_AB_plus_Ross_Section=True
	plot_EastA_plus_Adele_Section=True
	plot_Weddell_Section=True
	plot_depth_int_fields=False

	plot_all_panels_without_data=False   #False makes the real figure
	

	#Parameters;
	pole='south'
	title=None
	cscale=None
	#variable_list=np.array(['temp', 'salt', 'u'])
	variable_list=np.array(['temp', 'salt', 'u','rho','CN','HI'])
	#yitle_list=np.array(['Delta1 (Area=0.0026km^2)', 'Delta6 (Area=0.029km^2)', 'Delta9 (Area=1.8km^2)'])
	#title_list = {'Delta1': 'Area=0.0026km^2', 'Delta3':'Area=0.029km^2', 'Delta4':'Area=0.12km^2' , 'Delta6': 'Area=0.35km^2' , 'Delta9': 'Area=1.8km^2'}
	title_list = {'temp': 'Temperature', 'salt':'Salinity', 'u':'Zonal Velocity','rho':'Density','CN':'Sea Ice Concentration','HI':'Sea Ice Thickness'}
	ocean_list = {'temp': 'temp', 'salt':'salt', 'u':'u','rho':'rho','CN':'CN','HI':'HI'}
	scale_list = {'temp': 'temp', 'salt':'salt', 'u':'u','rho':'SSD','CN':'CN','HI':'HI'}
	sea_surface_list = {'temp': 'SST', 'salt':'SSS', 'u':'UO','rho':'SSD','CN':'CN','HI':'HI'}
	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	letter_labels=np.array(['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'])
	#y_title_list=np.array(['Sea Surface Salinity', 'AB+Ross', 'East Antarctia', 'Weddell','Depth Integrated Salinity'])
	y_title_list=np.array(['Sea Surface', 'AB+Ross', 'East Antarctia', 'Weddell'])
	#y_title_list=np.array(['SSS','Column Salinity', 'AB+Ross', 'East Antarctia', 'Weddell'])
	colorbar_unit_list=np.array(['Temperature (C)', 'Salinity (psu)','Zonal Velocity (m/s)','Density (kg/m$^{3}$)','Sea Ice Concentration (non dim)', 'Sea Ice Thickness (m)'])
	nrows=2
	ncols=3#len(variable_list)#3
	temperature_scale=20
	salinity_scale=30
	velocity_scale=0.2
	
	datapath='/home/Alon.Stern/Iceberg_Project/iceberg_scripts/python_scripts/size_matters_paper/processed_data/'
        
	#Setting up the figure
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	
	#This is used to plot the figure without the data for playing this anotations.
	if plot_all_panels_without_data==True:
		save_figure=False
	else:

		#Plotting Sea Surface Variables
		if plot_Sea_Surface_fields==True:
			cscale=None
			#cscale=salinity_scale
			row_num=1
			difference_on=0
			axes_direction='xy'
			field_name='Sea Surface Variables'
			file_extension= '_1959_to_2019_all_south_xy.mat'
			produce_figure_variable_row(fig,axes,nrows,ncols,sea_surface_list, title_list,colorbar_unit_list,y_title_list, datapath, file_extension,cscale, \
					row_num,difference_on, pole,field_name,variable_list,axes_direction,letter_labels)

		#Plotting Depth Integrated Variables
		if plot_depth_int_fields==True:
			#cscale=40
			cscale=None
			row_num=2
			difference_on=0
			axes_direction='xy'
			field_name='Depth Integrated'
			file_extension='_1959_to_2019_all_south_xy.mat'
			produce_figure_variable_row(fig,axes,nrows,ncols,ocean_list, title_list,colorbar_unit_list,y_title_list, datapath, file_extension,cscale, \
					row_num,difference_on, pole,field_name,variable_list,axes_direction,letter_labels)

			

	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
	fig.set_size_inches(13.5, 13.5,forward=True)
	if save_figure==True:
		plt.savefig('paper_figures/Fig_control_mean_CN_HI_S_T_U_rho.png',dpi=300,bbox_inches='tight')

	plt.show()


if __name__ == '__main__':
	main()
	#sys.exit(main())

