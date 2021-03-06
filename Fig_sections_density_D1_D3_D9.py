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
	p_values=mat_contents['p_values'][:,:]
	return [lat, lon, data,p_values,field,layer_interface,lat_lon_bounds ]

def produce_figure_row(fig,axes,nrows,ncols,Delta_list, title_list,colorbar_unit_list,y_title_list,datapath, file_extension,cscale, row_num, \
		difference_on,pole,field_name,axes_direction='xy', significant_only='True'):
			title=None
			letter_labels=np.array(['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'])
			for k in range(ncols):
				print 'Plotting ' + field_name + ' fields for ' + Delta_list[k]
				ax=subplot(nrows,ncols,(row_num-1)*ncols+(k+1))
				filename=datapath + Delta_list[k] + file_extension
				#print filename
				(lat ,lon, data, p_values,field,layer_interface,lat_lon_bounds)=load_data_from_mat_file(filename)
				if significant_only==False:
					p_values=None
				if axes_direction=='xy':
					(cscale,datamap)=plot_polar_field(lat,lon,data,pole,difference_on,title,p_values,cscale,field,colorbar_on=False,return_data_map=True)
				else:
					file_type='ocean_month_z'
					datamap=plot_vertical_section(lon,lat,data,difference_on,title,p_values,cscale,field,layer_interface,axes_direction,lat_lon_bounds,file_type,\
							colorbar_on=False,return_data_map=True)
				
				if row_num==1:
					#ax.set_title(title_list[k], fontsize=18, y=1.13)
					ax.set_title(title_list[Delta_list[k]], fontsize=18, y=1.13)
				
				#if k==0:
				#	plt.ylabel(y_title_list[row_num-1], fontsize=18, labelpad=25)
				if k>0:
					ax.set_yticks([])
				
				if k==0 and row_num<4:
					ax.annotate(y_title_list[row_num-1], xy=(-0.6, 0.5), xycoords='axes fraction', fontsize=21,rotation=90,verticalalignment='center')
					#ax.annotate('Sector', xy=(-0.45, 0.5), xycoords='axes fraction', fontsize=21,rotation=90,verticalalignment='center')
				#	plt.ylabel(y_title_list[i], fontsize=18, labelpad=25)
				if k==0 and row_num>0:
					plt.ylabel('depth (m)', fontsize=14)
					ax.annotate(y_title_list[row_num-1], xy=(-0.6, 0.5), xycoords='axes fraction', fontsize=21,rotation=90,verticalalignment='center')
				
				if row_num!=nrows:
					ax.set_xticks([])
				else:
					plt.xlabel('latitude (deg)', fontsize=14)
					ax.set_xticks([-40 ,-60, -80])

				text(1,1,letter_labels[(row_num-1)*ncols+(k)], ha='right', va='bottom',transform=ax.transAxes)
			
			#Creating colorbar
			fig.subplots_adjust(right=0.8)
			cbar_ax = fig.add_axes([0.85, 0.12+((nrows-row_num)*0.28), 0.033, 0.2])
			cbar=fig.colorbar(datamap, cax=cbar_ax)
			cbar.set_label(colorbar_unit_list[row_num-1], rotation=90)
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
	significant_only=True
	plot_SSD_fields=False
	plot_depth_int_rho_fields=False
	plot_AB_plus_Ross_Density_Section=True
	plot_EastA_plus_Adele_Density_Section=True
	plot_Weddell_Density_Section=True

	plot_all_panels_without_data=False   #False makes the real figure
	

	#Parameters;
	pole='south'
	title=None
	cscale=None
	Delta_list=np.array(['Delta1', 'Delta3', 'Delta9'])
	#yitle_list=np.array(['Delta1 (Area=0.0026km^2)', 'Delta6 (Area=0.029km^2)', 'Delta9 (Area=1.8km^2)'])
	#title_list = {'Delta1': 'Area=0.0026km^2', 'Delta3':'Area=0.029km^2', 'Delta4':'Area=0.12km^2' , 'Delta6': 'Area=0.35km^2' , 'Delta9': 'Area=1.8km^2'}
	title_list = {'Delta1': 'Length=60m', 'Delta3':'Length=200m', 'Delta4':'Length=350m$' , 'Delta6': 'Length=700m$' , 'Delta9': 'Length=1600m'}
	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	
	#Change this to switch to one slice
	y_title_list=np.array(['AB+Ross', 'East Antarctia', 'Weddell Sea'])	;  Sector_list=np.array(['AB_plus_Ross', 'EastA_plus_AdeleLand','Weddell'])
	#y_title_list=np.array(['120West', '52West','60East']) ;	Sector_list=np.array(['120West', '52West','60East'])

	colorbar_unit_list=np.array(['density (kg/m$^{3}$)', 'density (kg/m$^{3}$)','density (kg/m$^{3}$)'])
	nrows=3
	ncols=len(Delta_list)#3
	density_scale=0.03
	
	datapath='/home/Alon.Stern/Iceberg_Project/iceberg_scripts/python_scripts/size_matters_paper/processed_data/'
        
	#Setting up the figure
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	
	#This is used to plot the figure without the data for playing this anotations.
	if plot_all_panels_without_data==True:
		save_figure=False
		count=0
		cscale=None
		difference_on=0
		boundinglat=-45
		projection = 'splaea'
		m = Basemap(projection=projection, boundinglat=boundinglat, lon_0=180)
		#filename=datapath + 'Delta9_melt_1959_to_2019_all_south_xy.mat'
		filename=datapath + 'Delta9_rho_1959_to_2019_all_south_yz_anomaly_EastA_plus_AdeleLand.mat'
		(lat ,lon, data, p_values,field ,layer_interface,lat_lon_bounds )=load_data_from_mat_file(filename)
		datamap=plot_vertical_section(lon,lat,data,difference_on,title,p_values,cscale,field,layer_interface,'yz',lat_lon_bounds,'ocean_month_z',\
							                                                        colorbar_on=False,return_data_map=True)
		for i in range(nrows):
			for j in range(ncols):
				count=count+1
				ax=subplot(nrows,ncols,count)
				print 'Working on subplot ' ,str(count)

				if i!=0 and i!=1:
					datamap=plot_vertical_section(lon,lat,data,difference_on,title,p_values,cscale,field,layer_interface,'yz',lat_lon_bounds,'ocean_month_z',\
							                                                        colorbar_on=False,return_data_map=True)
				else:
					plot_polar_field(lat,lon,None,pole,difference_on=0.,title=None)
				
				if i==0:
					#ax.set_title(title_list[j], fontsize=18, y=1.13)
					ax.set_title(title_list[Delta_list[j]], fontsize=18, y=1.13)
				if j==0 and i<2:
					ax.annotate(y_title_list[i], xy=(-0.6, 0.5), xycoords='axes fraction', fontsize=21,rotation=90,verticalalignment='center')
					ax.annotate('Density', xy=(-0.45, 0.5), xycoords='axes fraction', fontsize=21,rotation=90,verticalalignment='center')
				#	plt.ylabel(y_title_list[i], fontsize=18, labelpad=25)
				if j==0 and i>1:
					plt.ylabel('depth (m)', fontsize=14)
					ax.annotate(y_title_list[i], xy=(-0.6, 0.5), xycoords='axes fraction', fontsize=21,rotation=90,verticalalignment='center')

				if i!=nrows-1:
					ax.set_xticks([])
				else:
					plt.xlabel('latitude (deg)', fontsize=14)
					ax.set_xticks([-40 ,-60, -80])

				if j>0:
					ax.set_yticks([])

				if j==ncols-1 :
					fig.subplots_adjust(right=0.8)
					cbar_ax = fig.add_axes([0.85, 0.11+((nrows-i-1)*0.165), 0.03, 0.13])
					cbar=fig.colorbar(datamap, cax=cbar_ax)
					cbar.set_label(colorbar_unit_list[i], rotation=90)
	
	else:

		#Plotting Melt
		#Plotting AB_plus_Ross Density
		if plot_AB_plus_Ross_Density_Section==True:
			cscale=density_scale
			row_num=1
			difference_on=1
			axes_direction='yz'
			field_name='AB_plus_Ross Density Section'
			#file_extension='_rho_1959_to_2019_all_south_yz_anomaly_AB_plus_Ross.mat'
			file_extension='_rho_1959_to_2019_all_south_yz_anomaly_'+ Sector_list[row_num-1]  +'.mat'
			produce_figure_row(fig,axes,nrows,ncols,Delta_list, title_list,colorbar_unit_list,y_title_list, datapath, file_extension,cscale, \
					row_num,difference_on, pole,field_name,axes_direction,significant_only)

		#Plotting EastA_plus_Adele Density
		if plot_EastA_plus_Adele_Density_Section==True:
			cscale=density_scale
			row_num=2
			difference_on=1
			axes_direction='yz'
			field_name='EastA_plus_Adele Density Section'
			#file_extension='_rho_1959_to_2019_all_south_yz_anomaly_EastA_plus_AdeleLand.mat'
			file_extension='_rho_1959_to_2019_all_south_yz_anomaly_'+ Sector_list[row_num-1]  +'.mat'
			produce_figure_row(fig,axes,nrows,ncols,Delta_list, title_list,colorbar_unit_list,y_title_list, datapath, file_extension,cscale, \
					row_num,difference_on, pole,field_name,axes_direction,significant_only)
		#Plotting Weddell Density
		if plot_Weddell_Density_Section==True:
			cscale=density_scale
			row_num=3
			difference_on=1
			axes_direction='yz'
			field_name='Weddell Density Section'
			#file_extension='_rho_1959_to_2019_all_south_yz_anomaly_Weddell.mat'
			file_extension='_rho_1959_to_2019_all_south_yz_anomaly_'+ Sector_list[row_num-1]  +'.mat'
			produce_figure_row(fig,axes,nrows,ncols,Delta_list, title_list,colorbar_unit_list,y_title_list, datapath, file_extension,cscale, \
					row_num,difference_on, pole,field_name,axes_direction,significant_only)
			

			

	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
	fig.set_size_inches(13.5, 20.5,forward=True)
	if save_figure==True:
		plt.savefig('paper_figures/Fig_sections_density_D1_D3_D9.png')
	plt.show()


if __name__ == '__main__':
	main()
	#sys.exit(main())

