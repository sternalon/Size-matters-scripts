#!/usr/bin/env python

import sys
import argparse
import netCDF4 as nc
import numpy as np
import scipy.io as sc
from scipy import stats
import scipy.interpolate
import datetime
import matplotlib
import matplotlib as mpl
import numpy.ma as ma
from pylab import *
#matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from numpy.ma import masked_array
from matplotlib.colors import Normalize
from iceberg_mass_comparisons import find_max_min_years
from iceberg_mass_comparisons import define_mask
from m6toolbox  import section2quadmesh
from m6toolbox  import rho_Wright97

class MidpointNormalize(Normalize):
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		Normalize.__init__(self, vmin, vmax, clip)
	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y))


def get_max_within_area(field, lats, latitude, pole='north'):
	"""
	Get the maximum value of field within the area defined by latitude and 'pole'.
	"""
	assert(pole == 'south' or pole == 'north')
        idx = (np.abs(lats[:,0] - latitude)).argmin()
        if pole == 'north':
	        return np.max(np.abs(field[idx:,:]))
	else:
		return np.max(np.abs(field[:idx,:]))


def get_mean_abs_within_area(field, lats, latitude, pole):
	"""
	Get the maximum value of field within the area defined by latitude and 'pole'.
	"""
	assert(pole == 'south' or pole == 'north')
        idx = (np.abs(lats[:,0] - latitude)).argmin()
        if pole == 'north':
	        return np.mean(np.abs(field[idx:,:]))
	else:
		return np.mean(np.abs(field[:idx,:]))





def plot_polar_field(lat,lon,data,pole,difference_on=0.,title=None,p_values=None,cscale=None,field=None,colorbar_on=True,return_data_map=False,plot_lat_lon_lines=True,boundinglat=-57.5):

	if pole == 'north':
	        projection = 'nplaea'
		boundinglat = 45.
	else:
		projection = 'splaea'
                #boundinglat = -57.5




	fig = plt.figure(1)
	lons, lats = np.meshgrid(lon, lat)
	
	if pole!='world':
		m = Basemap(projection=projection, boundinglat=boundinglat, lon_0=180)
		if plot_lat_lon_lines==True:
			m.drawparallels(np.arange(-80.,81.,20.), zorder=1)
			m.drawmeridians(np.arange(-180.,181.,20.), zorder=1)

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


	m.drawcoastlines(linewidth=.6, zorder=2)
	m.fillcontinents(color='0.8')
	x, y = m(lons, lats)
	
	if data is not None:
		
		if difference_on==1:
			if cscale==None:
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
				elif pole=='north' or pole=='south' :   # For non global cases which have not yet been defined
					cscale =4*round(get_mean_abs_within_area(data, lats, boundinglat, pole)*100)/100
			
			cNorm = MidpointNormalize(vmin=-cscale, vmax=cscale,midpoint=0)
			cmap = 'bwr'

		if difference_on==0:
			#if cscale==None and pole!='world':   # For non global cases which have not yet been defined
			#	cscale =4*round(get_mean_abs_within_area(data, lats, boundinglat, pole)*100)/100

			cmap = 'jet'
			if field=='ice' or field=='CN':
				cNorm = mpl.colors.Normalize(vmin=0, vmax=1)
			if field!='ice':
				vmax=np.max(data)  ; vmin=np.min(data)
				#cNorm = MidpointNormalize(vmin=-cscale, vmax=cscale,midpoint=0)
				cNorm = MidpointNormalize(vmin=vmin, vmax=vmax,midpoint=0)
			if field=='SSS':
				cNorm = mpl.colors.Normalize(vmin=33, vmax=35)
			if field=='SST':
				cNorm = mpl.colors.Normalize(vmin=-3, vmax=10)
			if field=='SSD':
				cNorm = mpl.colors.Normalize(vmin=1024, vmax=1028)
			if field=='HI':
				cNorm = mpl.colors.Normalize(vmin=0, vmax=1)
			if field=='UO':
				cscale=0.05
				cNorm = MidpointNormalize(vmin=-cscale, vmax=cscale,midpoint=0)
				#cmap = plt.get_cmap('bwr')
				cmap = 'bwr'
		#Special cases using log colors
		if field=='melt'  or  field=='melt_buoy' or field=='melt_eros' or  field=='melt_conv':
			#cNorm = mpl.colors.LogNorm(vmin=10**-8, vmax=10**-4)
			#cNorm = mpl.colors.LogNorm(vmin=10**-9, vmax=10**-4)
			if difference_on==0:
				cNorm = mpl.colors.LogNorm(vmin=10**-11*(60*60*24*365.5/850), vmax=10**-4*(60*60*24*365.5/850)) # In meters per year.
			if difference_on==1:
				cNorm = MidpointNormalize(vmin=-0.1, vmax=0.1,midpoint=0)
			#cmap = 'Blues'
		if field=='iron':
			cNorm = mpl.colors.LogNorm(vmin=10**-11*((0.00069)*60*60*24*365.5/850), vmax=10**-4*((0.00069)*60*60*24*365.5/850)) # In meters per year.
			cmap = 'hot_r'
		if field=='mass':
			#cNorm = mpl.colors.LogNorm(vmin=10**-3, vmax=10**3)
			cNorm = MidpointNormalize(vmin=-5, vmax=5,midpoint=0)

		if p_values is None:	
			datamap = m.pcolor(x, y, data, zorder=0,cmap=cmap,norm=cNorm)
		
		else:
			#data_s=masked_array(data,p_values>0.01)
			#data_u=masked_array(data,p_values<0.01)
			#cmap='seismic'
			#datamap = m.pcolor(x, y, data_s, cmap='seismic',norm=cNorm)
			#datamap2 = m.pcolor(x, y, data_u ,cmap='bwr',norm=cNorm)
			#data[np.where(p_values>0.05)]=np.NaN
			datamap = m.pcolor(x, y, data, zorder=0,cmap=cmap,norm=cNorm)
			levels1=np.array([0.05,999.])
			#cmap.set_bad('grey',1.)
			cNorm2=mpl.colors.Normalize(vmin=400, vmax=1000)
			CS = m.contourf(x, y, p_values,levels=levels1, hatches=[' '], fill=False,cmap='Greys',norm=cNorm2 ) # hatches=['+'] also an option,  hatches=['//']
			#CS = m.contourf(x, y, p_values,levels=levels1, hatches=[' '], fill=False,cmap='Greys_r' ) # hatches=['+'] also an option,  hatches=['//']
			#CS = m.contourf(x, y, p_values,levels=levels1, fill=True,cmap='bwr',color=np.array(['grey','grey']) ) # hatches=['+'] also an option,  hatches=['//']
			#CS = m.contour(x, y, p_values,levels=levels1,  fill=False,cmap='bwr' ) # hatches=['+'] also an option
		
		if colorbar_on==True:
			plt.colorbar(datamap, cmap=cmap, norm=cNorm, shrink=0.5)
			#plt.colorbar(datamap)
		if title!=None:
			plt.title(title)

	if return_data_map==False:
		return cscale
	if return_data_map==True:
		return [cscale,datamap]


def plot_vertical_section(lon,lat,data_mean,difference_on,title,p_values,cscale,field,layer_interface,axes_direction,lat_lon_bounds,file_type,colorbar_on=True,return_data_map=False):
	if axes_direction=='yz':
		x=lat
		upper_bound=lat_lon_bounds[0];
		lower_bound=lat_lon_bounds[1];

	if axes_direction=='xz':
		x=lon
		upper_bound=lat_lon_bounds[2];
		lower_bound=lat_lon_bounds[3];

	if file_type=='ocean_month':
		representation='linear'  #'pcm'
		#representation='pcm'
		M=layer_interface.shape
		layers=range(0,M[0]-1)
		q=data_mean[range(0,M[0]-1),:] ;q=q[:,range(0,len(x)-1) ]
		layer_interface=layer_interface[:,range(0,len(x)-1)]

		(X, Z, Q)=section2quadmesh(x, layer_interface, q, representation)
	
	if file_type=='ocean_month_z':
		X=x 
		Z=layer_interface

		Q=data_mean

	#Plotting
 	if difference_on==1:
		 #cmap = plt.set_cmap('bwr')
		 cmap = 'bwr'
		 if cscale==None:
			 if field=='temp':
				 cscale=0.5
			 if field=='salt':
				 cscale=0.1
			 if field=='u':
				 cscale=0.005
			 if field=='rho':
				 cscale=0.05
		 #cNorm = MidpointNormalize(vmin=-cscale, vmax=cscale,midpoint=0)
		 cNorm = mpl.colors.Normalize(vmin=-cscale, vmax=cscale)
	
 	if difference_on==0:
		#cmap = plt.set_cmap('jet')
		cmap = 'jet'
		if field=='temp':
			cscale_min=-3 ; cscale_max=10
	        if field=='rho':
			cscale_min=1026 ; cscale_max=1029
		if field=='salt':
			cscale_min=33 ; cscale_max=35
		#cmap = plt.set_cmap('jet')
		if field=='u':
			cscale_min=-.05 ; cscale_max=.05
		 	#cmap = plt.set_cmap('bwr')
		 	cmap = 'bwr'

		#cNorm = MidpointNormalize(vmin=cscale_min, vmax=cscale_max,midpoint=0)
		cNorm = mpl.colors.Normalize(vmin=cscale_min, vmax=cscale_max)
		#cmap = plt.set_cmap('jet')

	print X.shape
	print Z.shape
	print Q.shape
	datamap=pcolormesh(X,Z,Q,cmap=cmap,norm=cNorm,facecolor='white')
	
	#datamap=pcolormesh(X,-Z[::-1],Q[:,:],cmap=cmap,norm=cNorm,facecolor='white')
	#datamap=pcolormesh(X,-Z[::-1],Q[:,:],cmap=cmap,norm=cNorm,facecolor='white')

	
	#Overlaying significance
	
	Z2=Z[:-1]

	#if p_values!=None:
	if p_values is not None:
		#print 'plotting p_values'
		p_values[np.where(Q.mask)]=0.001
		levels1=np.array([0.05,999.])
		#levels1=np.array([0.01,0.02,0.03,0.04,0.05,999.])
		cNorm2=mpl.colors.Normalize(vmin=400, vmax=1000)
		CS = contourf(X, Z2, p_values,levels=levels1, hatches=[' '], fill=False,cmap='Greys',norm=cNorm2 ) # hatches=['+'] also an option,  hatches=['//']

		#CS = contourf(X, Z2, p_values,levels=levels1, hatches=['//'], fill=False,cmap='bwr' ) # hatches=['+'] ,hatches=['///] also an option
		#CS = contour(X, Z2, p_values,levels=levels1,linewidth=5.,color='black' ) # hatches=['+'] also an option
		#datamap=pcolormesh(X,Z,p_values,cmap=cmap,norm=cNorm,facecolor='white')
		#CS = m.contourf(x, y, p_values,levels=levels1,  fill=False,cmap='bwr' ) # hatches=['+'] also an option

	axis = plt.gca()
	landcolor=[.5,.5,.5]
	axis.set_axis_bgcolor(landcolor)
	plt.xlim( (lower_bound,upper_bound))
        #plt.ylim((-1500,0)  )
        plt.ylim((0 ,1500))
        #plt.ylim((-1500 ,1500))
	axis.invert_yaxis()
	#plt.gca.invert_yaxis()
        #plt.ylim((-5000,0)  )
	if title!=None:
		plt.title(title)
	if colorbar_on==True:
		plt.colorbar()
	
	if return_data_map==True:
		return datamap



def depth_integrate_and_time_average_ocean_data(file,data,months,field,file_type,h=None):
	if file_type=='ocean_month':
		if field!='h':
			h = np.squeeze( np.mean(file.variables['h'][months,:,:,:],axis=0) )
			data=np.squeeze(np.sum(data*h,axis=0))
		if field=='h':
			data=np.squeeze(np.sum(h,axis=0))
	
	if file_type=='ocean_month_z':
			data=np.squeeze(np.sum(data*h,axis=0))

	#data=np.squeeze(np.sum(h,axis=0))
	return data

	#layers=range(41,60)
	#data=np.squeeze(np.sum(data[layers,:,:]*h[layers,:,:],axis=0))
	#data=np.squeeze(np.sum(data[layers,:,:],axis=0))

def horizonal_averaging(file,data, field, months,lon,lat,section_axes,lat_lon_bounds,file_type):
	lon_inds=range(len(lon))
	lat_inds=range(len(lat))
	layer_interface=None

	if section_axes=='yz':
		lon_min=lat_lon_bounds[2]  ;  lon_max=lat_lon_bounds[3]
		if lon_min<lon_max:
			lon_inds=np.where((lon>lon_min) * (lon<lon_max))[:][0]
		else:
			lon_inds=np.where((lon>lon_min) + (lon<lon_max))[:][0]
		ax_num=2
	if section_axes=='xz':
		lat_min=lat_lon_bounds[0]  ;  lat_max=lat_lon_bounds[1]
		lat_inds=np.where((lat>lat_min)*(lat<lat_max))[:][0]
		ax_num=1
	data = data[months,:,:,:]
	data = data[:,:,lat_inds,:]
	data = data[:,:,:,lon_inds]

	if file_type=='ocean_month_z':
		data=np.squeeze( np.mean(data,axis=0))
		#data[data.mask]=np.nan
		#data = np.squeeze(   data[:,:,10])
		#data = np.squeeze( np.mean(  data[:,:,range(10,20)] , axis=ax_num))
		data = np.squeeze( np.mean(  data , axis=ax_num))

	if file_type=='ocean_month':
		data = np.squeeze( np.mean(  np.squeeze( np.mean(data,axis=0) ) , axis=ax_num))
		#print data.shape
		layer_interface = np.squeeze( np.mean(  np.squeeze( np.mean(file.variables['e'][months,:,lat_inds,lon_inds],axis=0) ) , axis=ax_num))
	return [data ,layer_interface]


def generate_data_file(start_year, end_year,input_folder,field,months_str,file_type,section_axes='xy',lat_lon_bounds=None):
	months=get_months_from_month_string(file_type,months_str)
	#Loop through the years
	layer_interface=None
	num_files=end_year-start_year+1
	count=-1
	for year in range(start_year, end_year+1):
		count=count+1
		if year >=1000:
			filename ='/' + str(year) + '0101.' +file_type + '.nc'
		elif year >= 100:
			filename ='/' + '0'+ str(year) + '0101.' +file_type + '.nc'
		elif year >= 10:
			filename ='/' + '00'+ str(year) + '0101.' +file_type + '.nc'
		elif year >= 0:
			filename ='/' + '000'+ str(year) + '0101.' +file_type + '.nc'

		input_file=input_folder + filename
		print input_file
		#Downloading the data
		with nc.Dataset(input_file) as file:
			if field=='rho':
				if file_type=='ocean_month'  or file_type=='ocean_month_z':
					T=file.variables['temp'][:,:,:,:]
					S=file.variables['salt'][:,:,:,:]
				else:
					T=file.variables['SST'][:,:,:]
					S=file.variables['SSS'][:,:,:]
				raw_data=rho_Wright97(S, T, P=0)
			else:
				raw_data=file.variables[field]

			if file_type!='ocean_month'  and file_type!='ocean_month_z':
				if field == 'CN':
					data = np.squeeze( np.mean(raw_data[months,:,:,:],axis=0) )
					data=np.squeeze(np.sum(data,axis=0))

				else:
				#if field != 'CN':
					data = np.squeeze( np.mean(raw_data[months,:,:],axis=0) )

				if count==0:
					dimen=data.shape
					All_data=np.zeros((num_files, dimen[0], dimen[1]))
					lon = file.variables['xT'][:] ; lon[-1]= lon[0]+360 
					lat = file.variables['yT'][:]
			else:
				if count==0:
					lon = file.variables['xh'][:] ; lon[-1]= lon[0]+360 
					lat = file.variables['yh'][:]
					h=None
					if file_type=='ocean_month_z':
						z = file.variables['zw'][:]

						if section_axes=='xy':   #setting up thickness grid of integrating in z
							dz = np.diff(file.variables['zw'][:])
							h=np.zeros((len(dz),len(lat),len(lon)))
							for i in range(len(lon)):
								for j in range(len(lat)):
									h[:,j,i]=dz

				if section_axes=='xy':
					if field=='h_ML':
						data = np.squeeze( np.mean(raw_data[months,:,:],axis=0) )
					else:
						data = np.squeeze( np.mean(raw_data[months,:,:,:],axis=0) )
						data=depth_integrate_and_time_average_ocean_data(file,data,months,field,file_type,h)
				if section_axes!='xy':
					(data, layer_interface) = horizonal_averaging(file,raw_data,field, months,lon,lat,section_axes,lat_lon_bounds,file_type)
					
					#data=np.squeeze(raw_data[1,:,:,0])
					#print data.mask
					#datamap=pcolormesh(lat,-z,data)
					#plt.show()

				if count==0:
					dimen=data.shape
					All_data=np.ma.zeros((num_files, dimen[0], dimen[1]))
					if section_axes!='xy' and file_type=='ocean_month':
						All_layer_interface=np.zeros((num_files, dimen[0]+1, dimen[1]))

			
			All_data[count,:,:]=data
			#Save layer thickness for vertical sections in isopycnal coordinates.
			if section_axes!='xy' and file_type=='ocean_month':
				All_layer_interface[count,:,:]=layer_interface
	
	
	if section_axes!='xy':
		if file_type=='ocean_month':
			layer_interface=np.squeeze(np.mean(All_layer_interface,axis=0))
		if file_type=='ocean_month_z':
			layer_interface=z


	return [All_data, lat, lon, layer_interface]


def plot_latitude_sum(lat,lon,data_mean,input_folder, start_year, pole, file_type):
	filename ='/' + str(start_year) + '0101.' +file_type + '.nc'
	input_file=input_folder + filename
	with nc.Dataset(input_file) as file:
		if file_type=='icebe/rgs_month' or file_type=='icebergs_day':
			area = file.variables['area'][:, :]
		if file_type=='ice_month' or file_type=='ice_day':
			area = file.variables['CELL_AREA'][:, :]
			#Defining the mask
			mask=define_mask(area,lat,pole)
			#mask=define_mask(area,lat,pole,lower_lat_bound,upper_lat_bound)
	data_total=data_mean*area#*mask
	data_sum=np.squeeze(np.sum(data_total,axis=0))
	
	M=data_mean.shape
	for k in range(M[0]):
		data_total[k,:]=data_sum

	plot_polar_field(lat,lon,data_total,pole,difference_on=1,title=None,p_values=None,cscale=None,field=None)
	#plt.plot(lon,data_sum)



def get_months_from_month_string(file_type,months_str):
	#Deciding which months to average over
	month_or_day=file_type.split('_')[1]
	if month_or_day=='month':
		if months_str=='all':
			months=range(12)
		if months_str=='jfm':
			months=np.array([0,1,2])
		if months_str=='jas':
			months=np.array([6,7,8])
		if months_str=='march':
			months=np.array([2])
		if months_str=='sept':
			months=np.array([8])

	if month_or_day=='day':
		if months_str=='all':
			months=range(365)
	
	return months


def produce_map_displaying_data(start_year, end_year,input_folder,field, second_folder, months_str, pole,file_type, axes_direction='xy',lat_lon_bounds=None,cscale=None,control_mean=None, \
		significant_vals_only=False,plot_lon_sum=False,difference_color_off=False,title=None,second_file_type=None,save_mat_file=False,lon_sector='all'):

	difference_on=1  # This flag allows you to plot the difference btw the first and second files
	if second_folder==None:
		difference_on=0
	if second_folder=='blank':
		difference_on=1
		second_folder=None

	#Collecting data first folder:
        (All_data, lat, lon, layer_interface) =generate_data_file(start_year, end_year, input_folder, field,months_str,file_type,axes_direction,lat_lon_bounds)
	data_mean=np.squeeze(np.mean(All_data,axis=0))

	#Collecting data second folder:
	p_values=None
	t_values=None
	if difference_on==1:
		if control_mean==None:
			if second_file_type==None:
				second_file_type=file_type
			(All_data2, lat, lon,layer_interface) =generate_data_file(start_year, end_year, second_folder, field, months_str,second_file_type,axes_direction,lat_lon_bounds)
			control_mean=np.squeeze(np.mean(All_data2,axis=0))
			if significant_vals_only=='True':
				#Calculating statistics:
				stats_vals=stats.ttest_ind(All_data,All_data2,axis=0, equal_var=False)
				t_values= stats_vals[0]
				p_values= stats_vals[1]
				p_values[np.where(np.isnan(p_values))]=999.
				t_values[np.where(np.isnan(t_values))]=999.	
		#data_mean=p_values
		#data_mean=t_values
		data_mean=np.squeeze(data_mean-control_mean)

	if control_mean==None:
		control_mean=data_mean
	#Plotting data:
	if title==None:
		print 'year',start_year, end_year
		title= str.split(input_folder,'_bergs_')[-1] + ' (' + str(start_year) + ' to ' + str(end_year) + ')'

	if difference_color_off==True:
		difference_on=0
	if plot_lon_sum==True:
		plot_latitude_sum(lat,lon,data_mean,input_folder,start_year,pole,file_type)
	else:
		if save_mat_file==True:
			save_the_mat_file(lat,lon,data_mean,p_values, field, difference_on,pole,title,layer_interface,file_type,lat_lon_bounds,cscale,axes_direction, \
			input_folder,start_year,end_year,months_str,lon_sector)

		if axes_direction=='xy':
			cscale=plot_polar_field(lat,lon,data_mean,pole,difference_on,title,p_values,cscale,field)
		if axes_direction!='xy':
			plot_vertical_section(lon,lat,data_mean,difference_on,title,p_values,cscale,field,layer_interface,axes_direction,lat_lon_bounds,file_type)

	return [control_mean, cscale]

def save_the_mat_file(lat,lon,data_mean,p_values, field, difference_on,pole,title,layer_interface,file_type,lat_lon_bounds,cscale,axes_direction,input_folder,\
		start_year,end_year,months_str,lon_sector ):
	#Setting absent variables to zero for saving purposes (can't save None)
	#if p_values==None:
	if p_values is None:
		p_values=0
	if layer_interface is None:
		#if layer_interface==None:
		layer_interface=0
	if cscale==None:
		cscale=0
	if field=='rho' and file_type=='ice_month':
		field='SSD'
	
	#Creating filename
	exp_name=str.split(input_folder,'_bergs_')[-1]
	if len(exp_name)>15:
		exp_name='Control'
	mat_filename='processed_data/'+exp_name+'_'+field+'_'+str(start_year)+'_to_'+ \
	str(end_year)+'_'+months_str+'_'+pole + '_' + axes_direction
	if difference_on==1:
		mat_filename=mat_filename+'_anomaly'
	if axes_direction!='xy' and lon_sector!='all':
		mat_filename=mat_filename+'_'+lon_sector
	mat_filename=mat_filename +'.mat'
	
	#Setting the mask equal to 1.e20 for saving reasons, to avoid the format being screwed up.
	#if axes_direction=='yz':
	fill_extra_values=False  #This is a hack.
	if fill_extra_values==True:
		ma.set_fill_value(data_mean,1.e20)
		data_mean=data_mean.filled()

	print p_values

	sc.savemat(mat_filename, {'lat':lat, 'lon':lon, 'data':data_mean,'p_values':p_values,'field':field, \
			'difference_on':difference_on, 'cscale':cscale,'pole': pole, 'title':title, 'layer_interface':layer_interface, \
			'file_type': file_type, 'axes_direction':axes_direction, 'lat_lon_bounds' : lat_lon_bounds})

	print 'File: ' + mat_filename + ' saved.'


def create_output_file_and_save(field,start_year,end_year,months_str,pole,plot_multiple_panels=True,input_folder=None,second_folder=None,naming_flag=None,save_fig=False):

	if plot_multiple_panels==1:
		input_folder=None
		second_folder=None
	#Creating name for output file:
	#label= str.split(input_folder,'_')[-1]
	output_file = 'figures/' + field + '_' + str(start_year) + '_to_' + str(end_year) + '_' + months_str + '_'+ pole

	if plot_multiple_panels==True:
		output_file=output_file+'_anomaly_mult_panels'
	else:
		if input_folder!=None:
			label= str.split(input_folder,'_')[-1]
			output_file=output_file+'_'+ label

		if second_folder!=None:
			output_file=output_file+ '_anomaly'
	
	if naming_flag!=None:
		output_file=output_file + '_' + naming_flag

	output_file= output_file + '.png'
	
	#Saving figure
	if save_fig==True:
		print 'Saving data in file: ' + output_file
		plt.savefig(output_file, dpi=150, bbox_inches='tight', pad_inches=0.4)


def define_longitudes(lon_sector):
	#Format for lat_lon_bounds:  (lat_min ,lat_max ,lon_min ,lon_max).
	#Also note, -279.5 < lon < 80.5

	#Definition of each sector:
	#lon_sector='all'  #Complete zonal average
	#lon_sector='AB'    #Amundsen and Bellinghausen
	#lon_sector='Ross'    #Ross only
	#lon_sector='Weddell'    #Weddell
	#lon_sector='EastA'    #East Anaractica
	#lon_sector='AdeleLand'    #AdeleLand

	lat_lon_bounds=np.array([-90 ,90 , -280 ,90])  #Format:  (lat_min ,lat_max  , lon_min  , lon_max)     Note, -279.5 < lon < 80.5
	if lon_sector=='AB':
		lat_lon_bounds[2]=-150;   lat_lon_bounds[3]=-55
	if lon_sector=='Ross':
		lat_lon_bounds[2]=-210 ;  lat_lon_bounds[3]=-150	
	if lon_sector=='Weddell':
		lat_lon_bounds[2]=-60 ;  lat_lon_bounds[3]=0
	if lon_sector=='EastA':
		lat_lon_bounds[2]=0 ;  lat_lon_bounds[3]=80
	if lon_sector=='AdeleLand':
		lat_lon_bounds[2]=-280 ;  lat_lon_bounds[3]=-210

	if lon_sector=='AB_plus_Ross':
		lat_lon_bounds[2]=-210;   lat_lon_bounds[3]=-55
	if lon_sector=='EastA_plus_AdeleLand':
		lat_lon_bounds[2]=0;   lat_lon_bounds[3]=-210

	if lon_sector=='120West':
		lat_lon_bounds[2]=-120;   lat_lon_bounds[3]=-118
	if lon_sector=='52West':
		lat_lon_bounds[2]=-53;   lat_lon_bounds[3]=-51
	if lon_sector=='60East':
		lat_lon_bounds[2]=59;   lat_lon_bounds[3]=61


	return lat_lon_bounds

##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument('--input_folder', default=None, help='The input data file in NetCDF format.')
	parser.add_argument('--second_folder', help='The input data file in NetCDF format.')
	parser.add_argument('--field',default='CN', help='The field to plot.')
	#parser.add_argument('--output_file',default = 'output_seaice_conc.png', help='The output figure name .png')
	parser.add_argument('--output_file', help='The output figure name .png')
    	parser.add_argument('--pole', default='south', help="""Plot South, North or both poles.""")
	parser.add_argument('--months', default='all', help="""Which models is the averaging done over: all, djf, jja.""")
	parser.add_argument('--file_type', default='ice_month',help='The fields included in this plot.')
	parser.add_argument('--cscale', default=None,help='Scale for the axes in the plot.')
	parser.add_argument('--lon_sector', default='all',help='Which longitude sector you are in for zonal averages.')
	parser.add_argument('--exp_name', default='delta',help='delta or groups are the experiment names.')
	parser.add_argument('--significance', default=True,help='Only plot significant values') 
   	args = parser.parse_args()
	
	#Cscale
	cscale=None
	if args.cscale!=None:
		cscale=float(args.cscale)
	
	#Flags
	plot_multiple_panels=0
	save_fig=False
	save_mat_file=True
	plot_lon_sum=False
	#significant_vals_only=False
	significant_vals_only=args.significance



	#If no imput file is used, then plot multiple panels.
	if args.input_folder==None:
		plot_multiple_panels=1

	#Parameters;
	start_year=1990
	end_year=2000
	#end_year=150

	#start_year=1950#1925
	#end_year=2015
	Number_of_years=10
	axes_direction='xy'
	if args.file_type=='ice_month' or args.file_type=='icebergs_month':
		axes_direction='xy'

	lat_lon_bounds=define_longitudes(args.lon_sector)

	if plot_multiple_panels==0:
		produce_map_displaying_data(start_year, end_year,args.input_folder,args.field,  args.second_folder, args.months, \
				args.pole,args.file_type,axes_direction,lat_lon_bounds,cscale ,control_mean=None,significant_vals_only=significant_vals_only,\
				plot_lon_sum=plot_lon_sum,difference_color_off=False,title='')
	 	naming_flag=''
	

	if plot_multiple_panels==1:
		root_path='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg'
		second_folder='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg'
		if args.file_type=='icebergs_month':
			second_folder=None
		second_file_type=args.file_type
		if args.exp_name=='Gold1' or args.exp_name=='Gold2' or args.exp_name=='Gold3':
			root_path='/ptmp/aas/model_output/ulm_201510_mom6_2015.10.23/MOM6_GOLD_SIS2'
			second_folder='/ptmp/aas/model_output/ulm_201510_mom6_2015.10.23/MOM6_GOLD_SIS2_bergs_none'
			#second_file_type='ice_month'


		#second_folder=None
		run_names=np.array(['']) ; naming_flag='Control'
		run_names=np.array(['Delta1', 'Delta3', 'Delta6', 'Delta9']) ; naming_flag='Delta'
		#run_names=np.array(['Delta1', 'Delta2', 'Delta3', 'Delta6', 'Delta9','Delta10']) ; naming_flag='Delta'

		if args.exp_name=='groups':
			run_names=np.array(['freq','all_big', 'mass', 'all_small']) ; naming_flag='groups'
		if args.exp_name=='delta':
			run_names=np.array(['Delta1', 'Delta3', 'Delta9']) ; naming_flag='Delta' ;
		if args.exp_name=='Delta1':
			run_names=np.array(['Delta1']) ; naming_flag='Delta' ;
		if args.exp_name=='Delta3':
			run_names=np.array(['Delta3']) ; naming_flag='Delta' ;
		if args.exp_name=='Delta9':
			run_names=np.array(['Delta9']) ; naming_flag='Delta' ;
		if args.exp_name=='tournadre':
			run_names=np.array(['tournadre']) ; naming_flag='Tournadre' ;
		if args.exp_name=='rolling':
			#run_names=np.array(['rolling_tournadre_Burton_new','rolling_tournadre_WM']) ; naming_flag='Rolling' ;
			run_names=np.array(['rolling_tournadre_Burton_new','rolling_tournadre_Burton_new']) ; naming_flag='Rolling' ;
			run_names=np.array(['rolling_tournadre_Burton_new']) ; naming_flag='Rolling' ;
			#run_names=np.array(['tournadre_C3']) ; naming_flag='Rolling' ;
		if args.exp_name=='Gold1':
			run_names=np.array(['none_AB','freq_AB', 'Delta1_AB', 'Delta10_AB']) ; naming_flag=args.exp_name #+ 'anom_from_AB';
		if args.exp_name=='Gold2':
			run_names=np.array(['none', 'freq','Delta1', 'Delta10']) ; naming_flag=args.exp_name #+ 'anom_from_AB';
		if args.exp_name=='Gold3':
			run_names=np.array(['freq','Delta1', 'Delta10']) ; naming_flag=args.exp_name #Subtracting control from each one.
			
		#run_names=np.array([ 'Delta9', 'freq','tournadre']) ; naming_flag='groups'
		#run_names=np.array(['Delta1' ,'Delta10','Delta1_thick', 'Delta10_thick', 'Delta1b_thick', 'Delta10b_thick']) ; naming_flag='Thick'
		if significant_vals_only==True:
			naming_flag=naming_flag+'_sig'
		if axes_direction!='xy':
			naming_flag=naming_flag+'_'+args.lon_sector
		control_mean=None
		cscale=None
		for k in range(1):
			print root_path
			input_folder=root_path + '_bergs_' + run_names[k]
			if args.exp_name=='Gold2':
				second_folder=root_path + '_bergs_' + run_names[k]
				input_folder=second_folder +'_AB'
			if args.exp_name=='Control':
				input_folder=second_folder
				second_folder=None
			if args.exp_name=='rolling':
				input_folder=root_path +'_bergs_'+ run_names[k]
				if k==0:
					#second_folder=root_path + '_bergs_'+ 'tournadre_C3'
					second_folder=root_path + '_bergs_'+ 'rolling_tournadre_Rolling_off'
					#second_folder=None
				#if k==1:
				#	second_folder=root_path + '_bergs_'+ 'tournadre_C3'
			print input_folder
			(min_year, max_year)=find_max_min_years(input_folder,'.' + args.file_type + '.nc')
			if second_folder!=None:
				(min_year2, max_year2)=find_max_min_years(second_folder,'.' + second_file_type + '.nc')
				min_year=max(min_year,min_year2) ; max_year=min(max_year,max_year2)  ;
			end_year=min(max_year,end_year)
			start_year=max(end_year-Number_of_years,min_year)
			print 'Start year: ', start_year,', End year: ', end_year
			print 'File: ' + input_folder +' from '+ str(start_year) + ' to ' + str(end_year)

			#subplot(2,2,k+1)
			subplot(1,1,k+1)
			control_mean=None # This line makes the contrl reload each time
			#(control_mean, cscale)=produce_map_displaying_data(start_year, end_year,input_folder,args.field, second_folder, args.months,args.pole,args.file_type,cscale,control_mean)
			(control_mean, cscale)=produce_map_displaying_data(start_year, end_year,input_folder,args.field, second_folder, args.months,args.pole,\
					args.file_type,axes_direction,lat_lon_bounds,cscale,control_mean,significant_vals_only,\
					plot_lon_sum=False,difference_color_off=False,title=None,second_file_type=second_file_type,save_mat_file=save_mat_file,lon_sector=args.lon_sector)
		

	if plot_multiple_panels==2: #Subtracting from the same run
		half_way_year=start_year+((end_year-start_year)/2)
		subplot(1,2,1)
		(control_mean, cscale)=produce_map_displaying_data(start_year,half_way_year ,args.input_folder,args.field, None , args.months, \
				args.pole,args.file_type,axes_direction,lat_lon_bounds,cscale ,control_mean=None,significant_vals_only=significant_vals_only,\
				plot_lon_sum=plot_lon_sum,difference_color_off=False,title='')
		cscale=None
		subplot(1,2,2)
		produce_map_displaying_data(half_way_year,end_year ,args.input_folder,args.field, 'blank' , args.months, \
				args.pole,args.file_type,axes_direction,lat_lon_bounds,cscale ,control_mean=control_mean,significant_vals_only=significant_vals_only,\
				plot_lon_sum=plot_lon_sum,difference_color_off=False,title='')
	 	naming_flag=''
	create_output_file_and_save(args.field,start_year,end_year,args.months,args.pole,plot_multiple_panels,args.input_folder,args.second_folder,naming_flag,save_fig)

	plt.show()


if __name__ == '__main__':
	    sys.exit(main())

