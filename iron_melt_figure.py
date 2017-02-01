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
	
	mass0=mat_contents['mass0'][0,:]

	field=mat_contents['field']

	return [lat ,lon,Total_berg_lon,Total_berg_lat,mass0]


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
	#if p_values==0:
	p_values=None
	return [lat, lon, data,p_values,field,layer_interface,lat_lon_bounds ]

def produce_figure_panel(fig,axes,nrows,ncols,Delta_list, title_list,colorbar_unit_list,y_title_list,datapath, file_extension,cscale, panel_num, \
		difference_on,pole,field_name,field_type_list,axes_direction='xy'):
			title=None
			
			print 'Plotting ' + field_name + ' fields for ' + Delta_list[0]
			ax=subplot(nrows,ncols,panel_num)
			filename=datapath + Delta_list[0] + file_extension
			#print filename
			(lat ,lon, data, p_values,field,layer_interface,lat_lon_bounds)=load_data_from_mat_file(filename)

			if field=='melt':
				data=data*(60*60*24*356.25) #Multiply to get units /year rather than /s
				data=data/850 # Multiply to get units of m/year
				data=data*(0.00069)
				field='iron'

			if axes_direction=='xy':
				(cscale,datamap)=plot_polar_field(lat,lon,data,pole,difference_on,title,p_values,cscale,field,colorbar_on=False,return_data_map=True)
			else:
				file_type='ocean_month_z'
				datamap=plot_vertical_section(lon,lat,data,difference_on,title,p_values,cscale,field,layer_interface,axes_direction,lat_lon_bounds,file_type,\
						colorbar_on=False,return_data_map=True)
			
			#if panel_num==1:
			#	ax.set_title(title_list[Delta_list[k]], fontsize=18, y=1.13)
			
			#Creating colorbar
			if panel_num>0:
				#cbar_ax=fig.add_axes([0.41+((panel_num-2)*0.28),0.15,0.2,0.05])
				#cbar=fig.colorbar(datamap, cax=cbar_ax, orientation="horizontal")
				cbar=fig.colorbar(datamap, orientation="horizontal")
				cbar.set_label(colorbar_unit_list[panel_num-1])
			#tick_locator = ticker.MaxNLocator(nbins=5)
			#tick_locator = ticker.MaxNLocator(nbins=5)
			#cbar.locator = tick_locator
			#cbar.update_ticks()


			#if panel_num==1:
			#	ax.subplots_adjust(bottom=0.5)



##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################

def main():

	parser = argparse.ArgumentParser()
   	args = parser.parse_args()
	
	
	#Flags
	save_figure=True
	significant_vals_only=False
	plot_traj_fields=False
	plot_melt_fields=True
	plot_CN_fields=False

	plot_all_panels_without_data=False   #False makes the real figure
	

	#Parameters;
	pole='south'
	title=None
	cscale=None
	Delta_list=np.array(['tournadre'])
	#yitle_list=np.array(['Delta1 (Area=0.0026km^2)', 'Delta6 (Area=0.029km^2)', 'Delta9 (Area=1.8km^2)'])
	#title_list = {'Delta1': 'Area=0.0026km^2', 'Delta3':'Area=0.029km^2', 'Delta4':'Area=0.12km^2' , 'Delta6': 'Area=0.35km^2' , 'Delta9': 'Area=1.8km^2'}
	title_list = {'Delta1': 'Length=60m', 'Delta3':'Length=200m', 'Delta4':'Length=350m$' , 'Delta6': 'Length=700m$' , 'Delta9': 'Length=1600m'}
	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	#y_title_list=np.array(['Sea Surface Density', 'AB+Ross', 'East Antarctia', 'Weddell','Depth Integrated Density'])
	y_title_list=np.array(['Sea Surface','Sea Surface', 'Sea Surface'])
	field_type_list=np.array(['Temperature','Salinity','Density'])
	#y_title_list=np.array(['SSD','Column Density', 'AB+Ross', 'East Antarctia', 'Weddell'])
	colorbar_unit_list=np.array(['Iron Concentration ($mol$/$m^2$)', 'melt (m/year)','Sea Ice Conc anomaly ', 'density (kg/$m^{3}$)','density (kg/$m^{3}$)'])
	nrows=1
	ncols=1
	temperature_scale=0.2
	salinity_scale=0.1
	density_scale=0.04
	
	datapath='/home/Alon.Stern/Iceberg_Project/iceberg_scripts/python_scripts/size_matters_paper/processed_data/'
        
	#Setting up the figure
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	
	#This is used to plot the figure without the data for playing this anotations.
	if plot_all_panels_without_data==True:
		save_figure=False
	else:

		#Plotting Panels


		#Plotting Trajectory Fields
		if plot_traj_fields==True:
			boundinglat=-45
			projection = 'splaea'
			panel_num=1

			print 'Plotting Traj Fields for ' + Delta_list[0]
			#ax=subplot(nrows,ncols,(row_num-1)*ncols+(k+1))
			filename=datapath + Delta_list[0] + '_traj_1994_to_1995_all_south.mat'
			(lat ,lon,Total_berg_lon,Total_berg_lat,mass0)=load_traj_data_from_mat_file(filename)

			m = Basemap(projection=projection, boundinglat=boundinglat, lon_0=180)
			ax=subplot(nrows,ncols,panel_num)
			plot_polar_field(lat,lon,None,pole,difference_on=0.,title=None)
			#m.scatter(Total_berg_lon,Total_berg_lat,1,marker='o',color='black')

			cNorm = mpl.colors.LogNorm(vmin=8.8e7, vmax=7.4e11)
			datamap=m.scatter(Total_berg_lon,Total_berg_lat, c=mass0, marker='o',cmap='jet',norm=cNorm)

			cbar=fig.colorbar(datamap, orientation="horizontal")
			cbar.set_label('Calving Mass (kg)')


			#ax.set_title(title_list[Delta_list[k]], fontsize=18,y=1.13)
			#plt.ylabel(y_title_list[k], fontsize=18, labelpad=25)


		if plot_melt_fields==True:
			cscale=None
			panel_num=1
			difference_on=0
			axes_direction='xy'
			field_name='melt'
			file_extension='_melt_1959_to_2019_all_south_xy.mat'

			produce_figure_panel(fig,axes,nrows,ncols,Delta_list, title_list,colorbar_unit_list,y_title_list, datapath, file_extension,cscale, \
					panel_num,difference_on, pole,field_name,field_type_list,axes_direction)

		if plot_CN_fields==True:
			cscale=None
			panel_num=3
			difference_on=1
			axes_direction='xy'
			field_name='CN'
			file_extension='_CN_1959_to_2019_all_south_xy_anomaly.mat'  #Remember to change this!!!

			produce_figure_panel(fig,axes,nrows,ncols,Delta_list, title_list,colorbar_unit_list,y_title_list, datapath, file_extension,cscale, \
					panel_num,difference_on, pole,field_name,field_type_list,axes_direction)












	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
	fig.set_size_inches(10.5, 10.5,forward=True)
	if save_figure==True:
		plt.savefig('iron_melt_figure.png')
	plt.show()


if __name__ == '__main__':
	main()
	#sys.exit(main())

