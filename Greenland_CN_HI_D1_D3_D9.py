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


def dunne_pm(N=256):
	"""
	Plus-minus  colormap from John Dunne.
	"""
	cdict = {'red':   [(0.00, 0.3, 0.3),
		(0.05, 0.5, 0.5),
		(0.20, 0.0, 0.0),
		(0.30, 0.4, 0.4),
		(0.40, 0.8, 0.8),
		(0.50, 1.0, 1.0),
		(0.95, 0.6, 0.6),
		(1.00, 0.4, 0.4)],

		'green': [(0.00, 0.0, 0.0),
			(0.30, 0.5, 0.5),
			(0.40, 1.0, 1.0),
			(0.70, 1.0, 1.0),
			(1.00, 0.0, 0.0)],

		'blue':  [(0.00, 0.3, 0.3),
			(0.05, 0.5, 0.5),
			(0.20, 1.0, 1.0),
			(0.50, 1.0, 1.0),
			(0.60, 0.7, 0.7),
			(0.70, 0.0, 0.0),
			(1.00, 0.0, 0.0)]}
	import matplotlib
	cmap = matplotlib.colors.LinearSegmentedColormap('dunnePM', cdict, N=N)
	cmap.set_under([.1,.0,.1]); cmap.set_over([.2,0.,0.])
	#cmap.set_bad('w')
	matplotlib.cm.register_cmap(cmap=cmap)
	return cmap


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
	field=mat_contents['field']
	title=mat_contents['title']
	p_values=mat_contents['p_values']
	return [lat, lon, data,p_values,field ]

def produce_figure_row(fig,axes,nrows,ncols,Delta_list, title_list,colorbar_unit_list,y_title_list,datapath, file_extension,cscale, row_num,\
		difference_on,pole,field_name,significant_only,y_title_list2):
			title=None
			letter_labels=np.array(['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'])
			for k in range(ncols):
				print 'Plotting ' + field_name + ' fields for ' + Delta_list[k]
				ax=subplot(nrows,ncols,(row_num-1)*ncols+(k+1))
				filename=datapath + Delta_list[k] + file_extension
				(lat ,lon, data, p_values,field)=load_data_from_mat_file(filename)
				#print min(p_values)
				if significant_only==False:
					p_values=None
				if field=='melt':
					data=data*(60*60*24*356.25) #Multiply to get units /year rather than /s
					data=data/850 # Multiply to get units of m/year
				(cscale,datamap)=plot_polar_field(lat,lon,data,pole,difference_on,title,p_values,cscale,field,colorbar_on=False,return_data_map=True)
				if k==0:
					plt.ylabel(y_title_list[row_num-1], fontsize=20, labelpad=25, multialignment='center')
					#ax.annotate(y_title_list[row_num-1], xy=(-0.4, 0.5), xycoords='axes fraction', fontsize=21,rotation=90,verticalalignment='center')
					#ax.annotate(y_title_list2[row_num-1], xy=(-0.25, 0.5), xycoords='axes fraction', fontsize=21,rotation=90,verticalalignment='center')

				if row_num==1:
					ax.set_title(title_list[Delta_list[k]], fontsize=20,y=1.13)
				text(1,1,letter_labels[(row_num-1)*ncols+(k)], ha='right', va='bottom',transform=ax.transAxes,fontsize=15)	
			#Creating colorbar
			fig.subplots_adjust(right=0.8)
			#cbar_ax = fig.add_axes([0.85, 0.103+((nrows-row_num)*0.165), 0.03, 0.13])
			cbar_ax = fig.add_axes([0.825, 0.115+((nrows-row_num)*0.435), 0.03, 0.33])
			if row_num==1:
	                	ticks=np.array([-0.04 ,-0.02 ,0 ,0.02,0.04])
			if row_num==2:
				ticks=np.array([-0.1 ,-0.05 ,0 ,0.05,0.1])
			cbar=fig.colorbar(datamap, cax=cbar_ax,ticks=ticks)
			cbar.set_label(colorbar_unit_list[row_num-1], rotation=90,fontsize=20)
			cbar.ax.tick_params(labelsize=20)
			#if field!='melt':
			#	tick_locator = ticker.MaxNLocator(nbins=5)
        	        #       cbar.locator = tick_locator
                	#        cbar.update_ticks()




##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################

def main():

	parser = argparse.ArgumentParser()
   	args = parser.parse_args()
	
	
	#Flags
	save_figure=True
	significant_only=True
	plot_traj_fields=False
	plot_melt_fields=False
	plot_sea_ice_conc_fields=True
	plot_sea_ice_thickness_fields=True
	plot_traj_AB_fields=False

	plot_all_panels_without_data=False   #False makes the real figure
	

	#Parameters;
	pole='north'
	title=None
	cscale=None
	Delta_list=np.array(['Delta1', 'Delta3', 'Delta9'])
	#Title_list=np.array(['Delta1 (Area=0.0026km^2)', 'Delta6 (Area=0.029km^2)', 'Delta9 (Area=1.8km^2)'])
	#Title_list=np.array(['Area=0.0026km^2', 'Area=0.35km^2', 'Area=1.8km^2'])
	#title_list = {'Delta1': 'Area=0.0026km^2', 'Delta3':'Area=0.029km^2', 'Delta4':'Area=0.12km^2' , 'Delta6': 'Area=0.35km^2' , 'Delta9': 'Area=1.8km^2'}
	title_list = {'Delta1': 'Length=60m', 'Delta3':'Length=200m', 'Delta4':'Length=350m$' , 'Delta6': 'Length=700m$' , 'Delta9': 'Length=1600m'}
	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	y_title_list=np.array(['Sea Ice Conc\nAnomaly', 'Sea Ice Thickenssn\nAnomaly'])
	#y_title_list=np.array(['Sea Ice Concentration', 'Sea Ice Thickenss'])
	y_title_list2=np.array(['Anomaly', 'Anomaly'])
	#colorbar_unit_list=np.array(['', 'melt (kg/m^2/s)', 'Conc (non-dim)', 'Thickness (m)',''])
	colorbar_unit_list=np.array(['Conc (non-dim)', 'Thickness (m)'])
	nrows=2
	ncols=len(Delta_list)#3
	
	datapath='/home/Alon.Stern/Iceberg_Project/iceberg_scripts/python_scripts/size_matters_paper/processed_data/'
        
	#Setting up the figure
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	
	#This is used to plot the figure without the data for playing this anotations.
	if plot_all_panels_without_data==True:
		save_figure=False
	
	else:
		#Plotting Trajectory Fields
		if plot_traj_fields==True:
			boundinglat=-45
			projection = 'splaea'
			row_num=1
			
			for k in range(ncols):
				print 'Plotting Traj Fields for ' + Delta_list[k]
				ax=subplot(nrows,ncols,(row_num-1)*ncols+(k+1))
				filename=datapath + Delta_list[k] + '_traj_1994_to_1995_all_south.mat'
				(lat ,lon,Total_berg_lon,Total_berg_lat)=load_traj_data_from_mat_file(filename)

				m = Basemap(projection=projection, boundinglat=boundinglat, lon_0=180)
				plot_polar_field(lat,lon,None,pole,difference_on=0.,title=None)
				m.scatter(Total_berg_lon,Total_berg_lat,1,marker='o',color=color_vec[k])
				ax.set_title(title_list[Delta_list[k]], fontsize=20,y=1.13)
				if k==0:
					plt.ylabel(y_title_list[k], fontsize=20, labelpad=25)
			
			fig.subplots_adjust(right=0.8)


		#Plotting Sea Ice Conc Fields
		if plot_sea_ice_conc_fields==True:
			cscale=0.04
			row_num=1
			difference_on=1
			field_name='Sea Ice Conc'
			file_extension='_CN_1959_to_2019_all_south_xy_anomaly.mat'
			produce_figure_row(fig,axes,nrows,ncols,Delta_list, title_list,colorbar_unit_list,y_title_list, datapath, file_extension,cscale, \
					row_num,difference_on, pole,field_name,significant_only,y_title_list2)
			

		#Plotting Sea Ice Thickness Fields
		if plot_sea_ice_thickness_fields==True:
			cscale=0.1
			row_num=2
			difference_on=1
			field_name='Sea Ice Thickness'
			file_extension='_HI_1959_to_2019_all_south_xy_anomaly.mat'
			produce_figure_row(fig,axes,nrows,ncols,Delta_list, title_list,colorbar_unit_list,y_title_list, datapath, file_extension,cscale, \
					row_num,difference_on, pole,field_name,significant_only,y_title_list2)
			
			
			fig.subplots_adjust(right=0.8)
	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
	fig.set_size_inches(18.5, 10.5,forward=True)
	if save_figure==True:
		plt.savefig('figures/Greenland_CN_HI_D1_D3_D9.png',dpi=300,bbox_inches='tight')
		#plt.savefig('paper_figures/Fig_CN_HI_D1_D3_D9.eps',bbox_inches='tight')
	plt.show()


if __name__ == '__main__':
	main()
	#sys.exit(main())

