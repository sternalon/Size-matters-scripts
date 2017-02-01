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


def plot_data_on_to_map(lat,lon,m,plot_four_maps,first_letter,ncols,starting_lon_lat,not_file_origin,fig):
	color_vec=np.array(['blue', 'red', 'green', 'purple', 'cyan', 'magenta', 'black','grey', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	#longitudes vary between -180<180

	starting_lon=starting_lon_lat[0]
	starting_lat=starting_lon_lat[1]
	plot_method='sector'
	#plot_method='letter'
	#count=0
	if plot_four_maps==True:
		if plot_method=='sector':
			lon1=-210
			if (starting_lon>lon1 and starting_lon<-55)  or (starting_lon>(lon1+360)):
				count=0
				subplot(1,ncols,1)
			if (starting_lon>0 and starting_lon<(lon1+360)) or starting_lon<lon1:
				count=1
				subplot(1,ncols,2)
			if starting_lon>-60 and starting_lon<0:
				count=2
				subplot(1,ncols,3)


		if plot_method=='letter':
			if first_letter=='b':
				count=0
				subplot(1,ncols,1)
			if first_letter=='c' or first_letter=='d':
			#if first_letter=='b' :
				count=1
				subplot(1,ncols,2)
			if first_letter=='a':
				count=2
				subplot(1,ncols,3)
			#if first_letter=='b':
			#	count=2
			#	subplot(1,ncols,4)

	if plot_four_maps==False:
		count=6
	#if lon[0]>lon_lower_bound and  lon[0]<lon_upper_bound:


	x, y = m(lon, lat)
	m.scatter(x,y,1,marker='o',color=color_vec[count])
	#if starting_lon==lon[0]:
	#if not_file_origin==False:
	m.scatter(x[0],y[0],3,marker='o',color='black')

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
			text(1,1,letter_labels[k], ha='right', va='bottom',transform=ax.transAxes)
	else:
		plot_polar_field(lat,lon,data_mean,pole,difference_on=0.,title='',\
				p_values=None,cscale=None,field=None,colorbar_on=True,return_data_map=False,plot_lat_lon_lines=True,boundinglat=boundinglat)

	m = Basemap(projection=projection, boundinglat=boundinglat, lon_0=180)
	#x, y = m(lon, lat)
	#x, y = m(0, -70)
	#m.scatter(x,y,3,marker='o',color='black')


	return [lat,lon,m]

def remove_repeated_values(W,L,F, lat, lon, SST, Time, codename):
	isgood=zeros([len(W),1])
	for k in range(0,len(W)-1):
			if W[k]==W[k-1]:
				isgood[k]=0.
			else:
				isgood[k]=1.
	print isgood[:,0]
	js=np.where((isgood[:,0]>0.5)*(W>0.))

	W=W[js]
	L=L[js]
	F=F[js]
	lat=lat[js]
	lon=lon[js]
	SST=SST[js]
	Time=Time[js]
	codename=codename[js]

	return [W,L,F, lat, lon, SST, Time, codename]

def load_data_from_tournadre(data_set_name,filter_repeats):
	print 'Loading data from matfile....'
	data_base = sc.loadmat('/home/aas/Iceberg_Project/Satelite_data/Tournadre_data/Tournadre_data.mat')
	
	lat=data_base['LA_I'] #Latitude 
	lon=data_base['LO_I'] #Longitude
	
	SST=data_base['SST_I'] #Longitude
	codename=data_base['INDn_I'] #Name of iceberg
	
	Time=data_base['T_I'] #Time in Julian Day
	#day=data_base['DAY_I'] #Day
	#month=data_base['MO_I'] #Month
	#year=data_base['AN_I'] #Year
	#Time=((year-2000)*365.25) + day

	if data_set_name=='NIC':
		W=data_base['W_I'] #Width from NIC analysis
		L=data_base['L_I'] #Length from NIC analysis
		F=data_base['H_I'] #Freeboard from NIC analysis

	if data_set_name=='Altimeter':
		W=data_base['WW_I'] #Width from NIC analysis
		L=data_base['LL_I'] #Length from NIC analysis
		F=data_base['H_I'] #Freeboard from NIC analysis

	if data_set_name=='Interpolated':
		W=data_base['wedf'] #Interpolated width
		L=data_base['ledf'] #Interpolated length
		F=data_base['hedf'] #Interpolated freeboard

	#Regridding
	j=np.where(np.ones([len(W),1])>0)
	W=W[j] ; L=L[j]  ; F=F[j]  ;SST=SST[j] ; Time=Time[j] ; lat=lat[j]  ; lon=lon[j] ;codename=codename[j]
	#W = [l[0] for l in W]
	#L = [l[0] for l in L]
	#F = [l[0] for l in F]
	#SST = [l[0] for l in SST]
	#Time = [l[0] for l in Time]
	#lat = [l[0] for l in lat]
	#lon = [l[0] for l in lon]
	#codename = [l[0] for l in codename]


	W=W*1000 ; L=L*1000   #Coverting to meters

	if filter_repeats==True:
		(W,L,F, lat, lon, SST, Time, codename)=remove_repeated_values(W,L,F, lat, lon, SST, Time, codename)


	print 'Finished Loading'

	return [W,L,F, lat, lon, SST, Time, codename]




###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
def main():
	#Clear screen
	#os.system('clear')

	
	#Plotting empty map  (perhaps not the most effitient way of doing this)
	download_tournadre_data=True
	ncols=1
	nrows=1
	save_figure=False
	#data_set_name='Interpolated'
	#data_set_name='Altimeter'
	data_set_name='NIC'

	#Parameters
	rho_i=800.
	rho_w=1020.


	#Flags
	filter_repeats_widths=True

        
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	#(lat,lon,m)=create_an_empty_map(plot_four_maps,ncols)

        #Setting up the figure
	
	(W,L,F, lat, lon, SST, Time, codename)=load_data_from_tournadre(data_set_name,filter_repeats_widths)
	H=F/(1-(rho_i/rho_w))
	D=F/((rho_w/rho_i)-1)
	Mass=rho_i*W*L*H
	N=len(Time)



	print 'There are ' ,N, 'many points left before filter'
	#Another filter
	time_diff=np.diff(Time)
	mass_diff=-np.diff(Mass)
	codename_diff=np.diff(codename)
	condition1=(abs(codename_diff)<0.5)
	condition2=(time_diff>1)
	condition3=(mass_diff>10.)
	js=np.where(condition1*condition2*condition3)
	Time=Time[js]
	Mass=Mass[js]
	codename=codename[js]
	N=len(Time)
	print 'There are ' ,N, 'many points left after filter'



	#N=5000
	#Mass=Mass[np.arange(0,N)]
	#Time=Time[np.arange(0,N)]
	#codename=codename[np.arange(0,N)]

	print 'There are ' ,N, 'many points left'
	time_diff_matrix=np.zeros([N,N])
	mass_diff_matrix=np.zeros([N,N])
	
	print 'Calculating the mass and time Matrix'
	for i in range(N):
		for j in range(i):
			if codename[i]==codename[j]:
				time_diff_matrix[i,j]=Time[i]-Time[j]
				#mass_diff_matrix[i,j]=-(Mass[i]-Mass[j])
				mass_diff_matrix[i,j]=(Mass[i]**(-1./3.)-Mass[j]**(-1./3))
	print 'done matrix'
	

	#One way of doing it	
	#time_diff=np.diff(Time)
	#mass_diff=-np.diff(Mass)
	#codename_diff=np.diff(codename)
	#condition1=(abs(codename_diff)<0.5)
	#condition2=(time_diff>1)
	#condition3=(mass_diff>10.)
	#js=np.where(condition1*condition2*condition3)
	
	condition2=(time_diff_matrix>1)
	condition3=(mass_diff_matrix>10.)
	#js=np.where(condition2*condition3)
	js=np.where(condition2)

	time_diff=time_diff_matrix[js]
	mass_diff=mass_diff_matrix[js]

	print 'There are ' ,len(time_diff), 'many points left'

	
	plot_mass_with_time=True
	if plot_mass_with_time==True:
		#plt.plot(Time,Mass,'o')
		#plt.plot(np.arange(len(Time)),Time,'o')
		plt.plot(time_diff,mass_diff,'o')
		plt.xlabel('Time (m)')
		plt.ylabel('Mass (kg)')
		#axes.set_yscale('log')
		#axes.set_xscale('log')



	plot_Width_vs_Height=False
	if plot_Width_vs_Height==True:
		plt.plot(W,H,'o')
		plt.xlabel('Width (m)')
		plt.ylabel('Thickness (m)')
	

	plot_rolling_equation=False
	if plot_rolling_equation==True:
		 D_vec=np.linspace(0.1,300,100)
		 L_vec=sqrt(0.92*(D_vec**2)+(58.32*D_vec))
		 L_MacAyeal=0.75*D_vec
		 L_MacAyeal=0.75*D_vec
		 L_Weeks=(0.49*D_vec) + np.sqrt(((0.49*D_vec)**2) +  (0.92*(D_vec**2)+(58.32*D_vec)) )
		 plt.plot(L_vec,D_vec,linestyle='--',color='blue')
		 plt.plot(L_MacAyeal,D_vec,linestyle=':',color='green')
		 plt.plot(L_Weeks,D_vec,linestyle='--',color='red')
	#axes.set_xscale('log')



	#plt.tight_layout()
	#fig.set_size_inches(19.5, 10.,forward=True)
        
	if save_figure==True:
		plt.savefig('observational_scatter.png')
	plt.show()



	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

