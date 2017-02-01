#!/usr/bin/env python
from netCDF4 import Dataset
import numpy as np
from PIL import *
from pylab import *
import math
import os
import matplotlib.pyplot as plt
import pdb
import argparse
import sys
import scipy.io as sc
import os.path
from os import listdir
from os.path import isfile, join
from scipy.interpolate import interp1d
from sea_ice_concentrations import plot_polar_field
from sea_ice_concentrations import generate_data_file
from mpl_toolkits.basemap import Basemap


def plot_data_on_to_map(lat,lon,m,plot_four_maps,first_letter,ncols,starting_lon_lat,not_file_origin,fig,Three_colors_on_one=False):
	color_vec=np.array(['blue', 'red', 'green', 'purple', 'cyan', 'magenta', 'black','grey', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	grey_vec=np.array(['blue', 'red', 'green', 'purple', 'cyan', 'magenta', 'black','grey', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	grey_vec=np.array(['darkgrey','lightgrey', 'lightgrey',  'cyan', 'magenta', 'black','grey', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	grey_vec=np.array(['mistyrose','honeydew', 'lightgrey',  'cyan', 'magenta', 'black','grey', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	grey_vec=np.array(['red','green', 'lightgrey',  'cyan', 'magenta', 'black','grey', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	#longitudes vary between -180<180
	x, y = m(lon, lat)

	starting_lon=starting_lon_lat[0]
	starting_lat=starting_lon_lat[1]
	plot_method='sector'
	#plot_method='letter'
	#count=0
	if (plot_four_maps==True) or (Three_colors_on_one==True):
		if plot_method=='sector':
			lon1=-210
			if (starting_lon>lon1 and starting_lon<-55)  or (starting_lon>(lon1+360)):
				count=0
			if (starting_lon>0 and starting_lon<(lon1+360)) or starting_lon<lon1:
				count=1
			if starting_lon>-60 and starting_lon<0:
				count=2


		if plot_method=='letter':
			if first_letter=='b':
				count=0
			if first_letter=='c' or first_letter=='d':
			#if first_letter=='b' :
				count=1
			if first_letter=='a':
				count=2
		
		if plot_four_maps==True:
			subplot(1,ncols,count+1)
			m.scatter(x,y,1,marker='o',color=color_vec[count])
		else:
			#m.scatter(x,y,1,marker='o',color=grey_vec[count],zorder=(3-count),lw=0.0)
			m.scatter(x,y,1,marker='o',color=grey_vec[count],alpha=0.2,zorder=(3-count),lw=0.0)

	if (plot_four_maps==False) and (Three_colors_on_one==False):
		decimate_time_series=False
		if decimate_time_series==True:
			N=10
			js=np.arange(0, len(x)-1, N)
			print js
			m.scatter(x[js[:]],y[js[:]],1,marker='o',color='lightgrey')
		else:
			m.scatter(x,y,0.4,marker='o',color='lightgrey',lw=0.0)

def create_an_empty_map(plot_four_maps=False,ncols=None):
	start_year=1947
	input_folder='/ptmp/Alon.Stern/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta2/'
	letter_labels=np.array(['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'])
	field='CN'
	months_str='all'
	file_type='ice_month'
	pole='south'
	(All_data, lat, lon,z) =generate_data_file(start_year, start_year,input_folder,field, months_str,file_type,section_axes='xy',lat_lon_bounds=None)
	projection = 'splaea'
	boundinglat=-45
	#boundinglat=-57.5
			
	data_mean=None

	#Plot blank map.
	if plot_four_maps==True:
		for k in range(ncols):
			subplot(1,ncols,k+1)
			plot_polar_field(lat,lon,data_mean,pole,difference_on=0.,title='',\
				p_values=None,cscale=None,field=None,colorbar_on=True,return_data_map=False,plot_lat_lon_lines=True,boundinglat=boundinglat)
			ax=gca()
			text(1,1,letter_labels[k], ha='right', va='bottom',transform=ax.transAxes,fontsize=15)
	else:
		plot_polar_field(lat,lon,data_mean,pole,difference_on=0.,title='',\
				p_values=None,cscale=None,field=None,colorbar_on=True,return_data_map=False,plot_lat_lon_lines=True,boundinglat=boundinglat)

	m = Basemap(projection=projection, boundinglat=boundinglat, lon_0=180)
	#x, y = m(lon, lat)
	#x, y = m(0, -70)
	#m.scatter(x,y,3,marker='o',color='black')


	return [lat,lon,m]


def plot_data_from_NIC_data_base(lat,lon,m,plot_four_maps,ncols,fig,Three_colors_on_one=False):
	root_path='/home/aas/Iceberg_Project/Iceberg_Satelite_data/'
	#dirname='ascat/'
	#all_file_types=np.array(['qscat', 'ascat', 'oscat', 'nscat',  'ers',  'sass' , 'seawinds'])
	all_file_types=np.array([ 'ascat', 'oscat', 'nscat',  'ers',  'sass' , 'seawinds','qscat'])
	#all_file_types=np.array([  'ers'])
	#all_file_types=np.array([  'ascat'])

	count=0
	for file_type in all_file_types:
		print 'Downloading files of type: ' ,file_type
		#file_type='qscat'
		mypath=root_path+file_type+'/'
		#filename=root_path +dirname + 'a23a.ascat'

		#Getting all file names
		file_list=[]
		for file in os.listdir(mypath):
			if file.endswith('.'+file_type):
				file_list.append(file)

		N=len(file_list)
		for k in range(N):
			#filename=root_path +file_type+ '/' +file_list[k]
			#if file_list[k][0:3]=='a22'  and file_type=='qscat':
			#if file_list[k][0]=='a':
			filename=file_list[k]

			first_letter=filename[0]  #first letter
			first_bit=filename.split('.')[0]
			second_bit=filename.split('.')[1]
			
			if first_bit[-1].isdigit()==True:
				origin_filename=mypath+filename
				not_file_origin=False
			else:
				origin_filename= mypath + first_bit[:-1] + '.' +second_bit
				if os.path.isfile(origin_filename)==False:
					origin_filename=mypath+filename
				not_file_origin=True

			#if first_letter=='a' or first_letter=='b' or first_letter=='c' or first_letter=='d':
			if True:
				filename=mypath +file_list[k]

				f = open(filename, 'r')
				lat=[]; lon=[]

				for line in f.readlines():
					line=line.strip()
					columns = line.split()
					if len(columns)>3:    #Incase the file has lines missing.
						lat.append(float(columns[1]))
						lon.append(float(columns[3]))
				f.close()
				
				starting_lon_lat=range(2)
				starting_lon_lat[0]=lon[0]
				starting_lon_lat[1]=lat[0]

				if origin_filename!=filename:
					f = open(origin_filename, 'r')
					lat_new=[]; lon_new=[]
					for line in f.readlines():
						line=line.strip()
						columns = line.split()
						if len(columns)>3:    #Incase the file has lines missing.
							lat_new.append(float(columns[1]))
							lon_new.append(float(columns[3]))
					f.close()
					starting_lon_lat[0]=lon_new[0]
					starting_lon_lat[1]=lat_new[0]

				plot_data_on_to_map(lat,lon,m,plot_four_maps,first_letter,ncols,starting_lon_lat,not_file_origin,fig,Three_colors_on_one)


def plot_data_from_tournadre(lat,lon,m,plot_four_maps):
	data_base = sc.loadmat('/home/aas/Iceberg_Project/Satelite_data/Tournadre_data/Tournadre_data.mat')
	#field_list=np.array(['LA_I','LO_I'])
	#for k in range(len(field_list)):
	LA_I=data_base['LA_I']
	LO_I=data_base['LO_I']
	print LO_I
	print LA_I
	print len(LA_I)
	print len(LO_I)
	x, y = m(LO_I, LA_I)
        m.scatter(x,y,1,marker='o',color='black')


###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
def main():
	#Clear screen
	#os.system('clear')

	
	#Plotting empty map  (perhaps not the most effitient way of doing this)
	plot_four_maps=True
	plot_NIC_data=True
	plot_tournadre_data=False
	Three_colors_on_one=False
	ncols=3
	nrows=1
	save_figure=True
        
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	(lat,lon,m)=create_an_empty_map(plot_four_maps,ncols)

        #Setting up the figure

	if plot_NIC_data==True and plot_tournadre_data==True:
		print 'Can not plot both data sets at once'
		return

	if plot_NIC_data==True:
		plot_data_from_NIC_data_base(lat,lon,m,plot_four_maps,ncols,fig,Three_colors_on_one)
	
	if plot_tournadre_data==True:
		plot_data_from_tournadre(lat,lon,m,plot_four_maps)

	plt.tight_layout()
	fig.set_size_inches(19.5, 10.,forward=True)
        
	if save_figure==True:
		plt.savefig('paper_figures/Fig_observational_iceberg_plotting.png')
	plt.show()



	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

