#!/usr/bin/env python
from netCDF4 import Dataset
import numpy as np
import netCDF4 as nc
from PIL import *
from pylab import *
import math
import os
import scipy.io as sc
import matplotlib.pyplot as plt
import pdb
import argparse
import sys
from scipy.interpolate import interp1d
from sea_ice_concentrations import MidpointNormalize
from sea_ice_concentrations import get_months_from_month_string
from mpl_toolkits.basemap import Basemap
from iceberg_mass_comparisons import find_max_min_years
from distributions_of_bergs import add_constraint_to_data
from sea_ice_concentrations import get_mean_abs_within_area


def fetch_field_from_file(filename, field, months_str, file_type, time_average,  section_axes):
	print 'Retriving data from file:  ' , filename
	months=get_months_from_month_string(file_type,months_str)
	#Loop through the years
	layer_interface=None

	with nc.Dataset(filename) as file:
		#Getting the lat and lon
		if file_type=='ocean_month_z' or  file_type=='ocean_month':
			lon = file.variables['xh'][:] ; lon[-1]= lon[0]+360 
			lat = file.variables['yh'][:]
		else:
			lon = file.variables['xT'][:] ; lon[-1]= lon[0]+360 
			lat = file.variables['yT'][:]

		if field=='rho':
			if file_type=='ocean_month'  or file_type=='ocean_month_z':
				T=file.variables['temp'][:]
				S=file.variables['salt'][:]
			else:
				T=file.variables['SST'][:]
				S=file.variables['SSS'][:]
			data=rho_Wright97(S, T, P=0)
		else:
			data=file.variables[field][:]
		
		if field == 'CN':
			data=np.squeeze(np.sum(data,axis=1))

	if time_average is True:
		data = np.squeeze( np.mean(data[months,:],axis=0) )
	
	

		print lat.shape
		print lon.shape
		print data.shape

	return [data, lat, lon]


def define_pojection(pole, lon, lat, bounding_lat=None):
	fig = plt.figure(1)
	lons, lats = np.meshgrid(lon, lat)
	
	if pole=='north':
		projection = 'nplaea'
		if bounding_lat is None:
			bounding_lat=45
		m = Basemap(projection=projection, boundinglat=bounding_lat, lon_0=180)
	if pole=='south':
		projection = 'splaea'
		if bounding_lat is None:
			bounding_lat=-45
		m = Basemap(projection=projection, boundinglat=bounding_lat, lon_0=180)

	if pole=='world':
		#m = Basemap(projection='gall', resolution='h', llcrnrlon=np.min(lons), llcrnrlat=np.min(lats), urcrnrlon=np.max(lons), urcrnrlat=np.max(lats), lon_0=0)
		m = Basemap(projection='gall',  llcrnrlon=np.min(lons), llcrnrlat=np.min(lats), urcrnrlon=np.max(lons), urcrnrlat=np.max(lats), lon_0=0)
		#m = Basemap(resolution='c',projection='ortho',lat_0=60.,lon_0=-60.)
		#m = Basemap(projection='gall',  lon_0=0)

	if pole=='Greenland':
		projection = 'lcc'
		lat_1=80
		lat_2=60
		lat_0=64.5
		lon_0=-50.
		m = Basemap(width=2700000,height=4500000, rsphere=(6378137.00,6356752.3142),resolution='l',area_thresh=1000.,\
				                                        projection=projection,lat_1=lat_1,lat_2=lat_2,lat_0=lat_0,lon_0=lon_0)
	if pole=='Greenland2':
		projection = 'lcc'
		lat_1=80
		lat_2=60
		lat_0=64.5
		lon_0=-50.
		m = Basemap(width=4500000,height=4500000, rsphere=(6378137.00,6356752.3142),resolution='l',area_thresh=1000.,\
				                                        projection=projection,lat_1=lat_1,lat_2=lat_2,lat_0=lat_0,lon_0=lon_0)

	if pole=='Peninsula':
		projection = 'lcc'
		lat_1=-80
		lat_2=-60
		lat_0=-64.5
		lon_0=-50.
		m = Basemap(width=3500000,height=4500000, rsphere=(6378137.00,6356752.3142),resolution='l',area_thresh=1000.,\
				                                        projection=projection,lat_1=lat_1,lat_2=lat_2,lat_0=lat_0,lon_0=lon_0)
	return m


def plot_data_on_map(data, lat, lon, m, cmap, cNorm, field):
        fig = plt.figurlons, lats = np.meshgrid(lon, lat)
        lons, lats = np.meshgrid(lon, lat)
	m.drawcoastlines(linewidth=.6, zorder=2)
        m.fillcontinents(color='0.8')
        m.drawparallels(np.arange(-80.,81.,20.), zorder=1)
        m.drawmeridians(np.arange(-180.,181.,20.), zorder=1)
	x, y = m(lons, lats)
	datamap = m.pcolor(x, y, data, zorder=0,cmap=cmap,norm=cNorm)
	plt.colorbar()
	plt.title(field)

def load_default_colorbar_values(field, cNorm, data, field_is_anomaly, lat, boundinglat, pole):
	if field_is_anomaly is True:
		cscale=None
		if field=='SSS':
			cscale=0.1
		elif field=='SST':
			cscale=0.5
		elif field=='UO':
			cscale=0.005
		elif field=='h_ML':
			cscale=15
		elif field=='HI':
			cscale=0.1
		elif field=='CN':
			cscale=0.04
		elif field=='melt'  or  field=='melt_buoy' or field=='melt_eros' or  field=='melt_conv':
			cscale=0.0000001
		elif field=='mass':
			cscale=0.1
		elif (lat is not None)	 and (boundinglat is not None)  and (pole is not None): 
        		lons, lats = np.meshgrid(lon, lat)
			cscale =4*round(get_mean_abs_within_area(data, lats, boundinglat, pole)*100)/100

		if cscale is not None:	
			cNorm = MidpointNormalize(vmin=-cscale, vmax=cscale,midpoint=0)

	if field_is_anomaly is False:
		if field=='ice' or field=='CN':
			cNorm = mpl.colors.Normalize(vmin=0, vmax=1)
		if field=='SSS':
			cNorm = mpl.colors.Normalize(vmin=33, vmax=35)
		if field=='SST':
			cNorm = mpl.colors.Normalize(vmin=-3, vmax=10)
		if field=='SSD':
			cNorm = mpl.colors.Normalize(vmin=1024, vmax=1028)
		if field=='HI':
			cNorm = mpl.colors.Normalize(vmin=0, vmax=1)
		if field=='UO':
			cNorm = MidpointNormalize(vmin=-0.05, vmax=0.05,midpoint=0)
		if field=='melt'  or  field=='melt_buoy' or field=='melt_eros' or  field=='melt_conv':
			cNorm = mpl.colors.LogNorm(vmin=10**-11*(60*60*24*365.5/850), vmax=10**-4*(60*60*24*365.5/850)) # In meters per year.
		if field=='mass':
			cNorm = mpl.colors.LogNorm(vmin=10**-3, vmax=10**3)

	return cNorm

def decide_on_color_scale(data, field, Centered_on_zero=False, Use_log_scale=False, vmax=None, vmin=None, use_set_default_values=True, field_is_anomaly=False, \
	lat=None, boundinglat=None, pole=None):
	if vmax==None:
		vmax=np.max(data)
		print 'Max value is : ' , vmax
		vmax=vmax/10.
	if vmin==None:
		vmin=np.min(data)
		print 'Min value is : ' , vmin
		vmin=vmin/10.

	if Centered_on_zero is True:
		vmax=max(abs(vmin),abs(vmax)) 
		vmin=-vmax
		cmap='bwr'
	else:
		cmap='jet'

	if Use_log_scale is True:
		cNorm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
	else:
		cNorm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)	
		#cNorm = MidpointNormalize(vmin=vmin, vmax=vmax,midpoint=0)
	
	if use_set_default_values is True:
		(cNorm)=load_default_colorbar_values(field, cNorm , data, field_is_anomaly, lat, boundinglat, pole)

	return [cNorm, cmap]

def decide_on_filetype(filename):
	if len(filename.split('ocean'))>1:
		file_type='ocean_month'
	elif len(filename.split('ice'))>1:
		file_type='ice_month'
	else:
		print 'Unknown filetype!!'
		halt
	return file_type

###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
def main():
	#Clear screen
	#os.system('clear')

	field='SSH'
	field='UO'
	field='v_accel_bt'
	field='TKE_mixing'
	#field='p_surf'
	#field='tauy'
	#field='depth_ocean'
	months_str='all'
	constraint_name=None
	pole='south'
	time_average=True
	#time_average=False
	section_axes='xy'
	bounding_lat=None
	cmap='jet'


	#Plotting flats
	Centered_on_zero=False
	Use_log_scale=False
	use_set_default_values=True
	field_is_anomaly=False
	vmax=None
	vmin=None
	#vmax=0.00000001
	#vmin=-vmax


	#Choose a filename
	root_path='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg/'
	#filename = root_path + 'time_averaged_from_1959_to_2019_ice_month.nc'
	#filename = root_path + 'time_averaged_from_1959_to_2019_ocean_month_z.nc'
	filename = root_path + 'time_averaged_from_1959_to_2019_ocean_month.nc'
	
	#Deciding on file_type
	file_type=decide_on_filetype(filename)
		
	(data, lat, lon) = fetch_field_from_file(filename, field, months_str, file_type, time_average , section_axes)
	m=define_pojection(pole, lon, lat,  bounding_lat)

	(cNorm, cmap) = decide_on_color_scale(data, field, Centered_on_zero, Use_log_scale, vmax, vmin, use_set_default_values, field_is_anomaly, lat, bounding_lat, pole)

	plot_data_on_map(data, lat, lon, m, cmap, cNorm, field)
		

	#plt.savefig(output_file, dpi=150, bbox_inches='tight', pad_inches=0.4)

	plt.show()



	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

