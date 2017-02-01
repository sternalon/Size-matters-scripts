#!/usr/bin/env python

import sys
import os
import argparse
import netCDF4 as nc
import numpy as np
import pandas as pd
import scipy.io as sc
import re
import datetime
import matplotlib
from pylab import *
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join

"""
What this script does:

Plot a full timeseries fields from mom ocean_scalar.nc

Example:

"""

def define_mask(area,lat,lon,pole,lat_lon_bounds):
	lat_lower_bound=lat_lon_bounds[0]
	lat_upper_bound=lat_lon_bounds[1]
	lon_lower_bound=lat_lon_bounds[2]
	lon_upper_bound=lat_lon_bounds[3]

	M=area.shape
	north_mask=np.zeros((M[0], M[1]))
	south_mask=np.zeros((M[0], M[1]))
	both_mask=np.ones((M[0], M[1]))
	south_mask[0:(M[0]/2),: ]=1.
	north_mask[((M[0]/2)+1): M[0],: ]=1.
	if pole=='south':
		mask=south_mask
	if pole=='north':
		mask=north_mask
	if pole=='both':
		mask=both_mask

	#Applying the latitude bound
	for k in range(len(lat)):
		if lat[k]<lat_lower_bound:
			mask[k,:]=0.
		if lat[k]>lat_upper_bound:
			mask[k,:]=0.


	#Applying the longitude bound
	if lon_lower_bound<lon_upper_bound:
		for k in range(len(lon)):
			if lon[k]<lon_lower_bound:
				mask[:,k]=0.
			if lon[k]>lon_upper_bound:
				mask[:,k]=0.
	else:
		for k in range(len(lon)):
			if lon[k]<lon_lower_bound and lon[k]>lon_upper_bound:
				mask[:,k]=0.
	
	return mask

def define_paths_array():
        all_paths=[]
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta1')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta2')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta3')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta4')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta5')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta6')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta7')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta8')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta9')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta10')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta12')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta13')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta11')
        
	all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_all_big')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_all_small')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_freq')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_mass')
        all_paths.append('/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg')
	return all_paths



def times_to_periods(time_var):
    """
    Return the time dimension as a list of pandas Periods.

    We use periods here because pandas doesn't support datetime objects as low
    as 0001-01-01.
    """

    m = re.search('\d{4}-\d{2}-\d{2}', time_var.units)
    start_date = datetime.datetime.strptime(m.group(0), '%Y-%m-%d')

    periods = None
    if 'hours since' in time_var.units:
        # Pandas likes the periods to be the same here, so make them both
        # daily.
        periods = [pd.Period(start_date + datetime.timedelta(days=int(h / 24.0)), 'D') for h in time_var[:]]
    elif 'days since' in time_var.units:
        periods = [pd.Period(start_date + datetime.timedelta(days=int(d)), 'D') for d in time_var[:]]
    else:
        assert(False)

    return periods


def calc_area_integral(data, area,mask,field,file_type):
	"""Calculate total ice area. """
	#if field!='calving':
	#	#if field=='mass':
	#	#	Total = np.sum(data *mask)
	#	#else:
	#	#	Total = np.sum(area * data *mask)
	Total = np.sum(area * data *mask)
	
	if file_type=='ice':
		#Note that the cell area in the ice_month field is measured in units of sphere.
		Area_of_earth=510.1*(10**12)  #in m^2
		Total=Total*Area_of_earth

	#print data.shape, area.shape
	if field=='calving':
		Total = np.sum( data *mask)


        return Total



	
def find_max_min_years(input_folder,file_type):
	year_list=[]
	for file in os.listdir(input_folder):
	    	if file.endswith(file_type):
			year=str.split(file,'0101.')[0]
			year_list.append(int(year))
	max_year=max(year_list)
	min_year=min(year_list)
	return [min_year, max_year]


def define_longitudes(lon_sector):  #Same as the one defined in sea_ice_concentration.py
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
		lat_lon_bounds[2]=-150;   lat_lon_bounds[3]=-60
	if lon_sector=='Ross':
		lat_lon_bounds[2]=-210 ;  lat_lon_bounds[3]=-150
	if lon_sector=='Weddell':
		lat_lon_bounds[2]=-60 ;  lat_lon_bounds[3]=0
	if lon_sector=='EastA':
		lat_lon_bounds[2]=0 ;  lat_lon_bounds[3]=80
	if lon_sector=='AdeleLand':
		lat_lon_bounds[2]=-280 ;  lat_lon_bounds[3]=-210

	if lon_sector=='AB_plus_Ross':
		lat_lon_bounds[2]=-210;   lat_lon_bounds[3]=-60
	if lon_sector=='EastA_plus_AdeleLand':
		lat_lon_bounds[2]=0;   lat_lon_bounds[3]=-210

	return lat_lon_bounds




def generate_time_series_data(start_year, end_year,input_folder,field,pole,file_type,av_months,lat_lon_bounds):
	data1 = []
	time1 = []
	
	#Make sure that data exists of the years allocated, othersize change it.
	(min_year, max_year)=find_max_min_years(input_folder,file_type + '_month.nc')
	if max_year<end_year:
		end_year=max_year

	# Go through input files one at a time building the timeseries as we go
	count=0
	#for file in input_files:
	for year in range(start_year, end_year+1):
		if year>=1000:
			filename ='/' + str(year) + '0101.' + file_type + '_month.nc'
		elif year>=100:
			filename ='/' + '0'+ str(year) + '0101.' + file_type + '_month.nc'
		elif year>=10:
			filename ='/' + '00' + str(year) + '0101.' + file_type + '_month.nc'
		elif year>=0:
			filename ='/' + '000'+ str(year) + '0101.' + file_type + '_month.nc'

		input_file=input_folder + filename
		print input_file
		count=count+1
		with nc.Dataset(input_file) as file:
			time_var = file.variables['time']
			if field=='Volume':
				data_var = file.variables['CN']
			else:
				data_var = file.variables[field]


			if count==1:
				title = data_var.long_name
				ylabel = data_var.units
				if file_type=='icebergs':
					area = file.variables['area'][:, :]
					lat =  file.variables['yT'][:]
					lon =  file.variables['xT'][:]
				if file_type=='ice':
					#print start_year
					area = file.variables['CELL_AREA'][:, :]
					#if start_year>500: #non-gold runs
					#	lat =  file.variables['yT'][:, :]
					#else:
					lat =  file.variables['yT'][:]
					lon =  file.variables['xT'][:]
				#Defining the mask
				mask=define_mask(area,lat,lon,pole,lat_lon_bounds)

			


			# Calculate the times/dates, these will be our indices.
			periods = times_to_periods(file.variables['time'])
			time = file.variables['time']

			if field=='CN':
				data=np.squeeze(np.sum(file.variables[field][:,:,:,:],axis=1))

			elif field=='Volume':
				data_CN=np.squeeze(np.sum(file.variables['CN'][:,:,:,:],axis=1))
				data_HI = file.variables['HI'][:,:,:]
				data=data_CN*data_HI

			else:
				data = file.variables[field][:,:,:]

			
	   		
	    		if av_months=='':
				for month in range(12):
					temp = calc_area_integral(np.squeeze(data[month,:,:]), area,mask,field,file_type)
					data1.append(temp)
	 				time1.append(time[month]/365.25)
			
			if av_months!='':
				if av_months!='max':
					temp = calc_area_integral(np.squeeze(np.mean(data[av_months,:,:],axis=0)), area,mask,field,file_type)
					data1.append(temp)
					time1.append(np.squeeze(np.mean(time[av_months]/365.25)))
				if av_months=='max':
					temp = calc_area_integral(np.squeeze(np.max(data[:,:,:],axis=0)), area,mask,field,file_type)
					data1.append(temp)
					time1.append(np.squeeze(np.mean(time[:]/365.25)))


	return [data1, time1]

def save_timeseries_matfile(time, data,Delta_name,field,lon_sector=''):
	mat_filename='processed_data/'+Delta_name+'_timeseries_'+ field
	if lon_sector!='':
		mat_filename=mat_filename+'_'+lon_sector
	mat_filename=mat_filename +'.mat'
	sc.savemat(mat_filename, {'time':time, 'data':data, 'Delta_name':Delta_name })
	print 'File: ' + mat_filename + ' saved.'



##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################


def main():
	parser = argparse.ArgumentParser()
        parser.add_argument('--field', default='mass',help='The fields included in this plot.')
        parser.add_argument('--file_type', default='icebergs',help='The fields included in this plot.')
        parser.add_argument('--av_months', default='', help='Directory where plots will be written.')
        parser.add_argument('--lon_sector', default='', help='Directory where plots will be written.')
        parser.add_argument('--exp_name', default='', help='Directory where plots will be written.')
        args = parser.parse_args()
        print 'Plotting the %s distribution' %args.field

	#Flags
	plot_time_series=1
	plot_averaged_values=0
        pole='south'
        #pole='north'
	start_year=1900
	end_year=1905
	#start_year=1
	#end_year=150
	save_figure=False
	save_matfile=False

	#Setting lower and upper latitude bounds
	lat_lon_bounds=define_longitudes(args.lon_sector)
	lat_lon_bounds[0]=-90  #lower bound latitude
	lat_lon_bounds[1]=0    #upper bound latitude


	#Parameters
	Q=3.1e6 #Approximate calving flux
	rho_i=860 #density of iceberg

	#List of directories:
	all_paths=define_paths_array()

	initial_mass_vec=np.array([ 8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11, 1.05e13 ])
	initial_thickness_vec=np.array([ 40., 67., 133., 175., 250., 250., 250., 250., 250., 250, 250 ])
	initial_area_vec=initial_mass_vec/(initial_thickness_vec*rho_i*np.pi)
	initial_radius_vec=np.sqrt(initial_mass_vec/(initial_thickness_vec*rho_i*np.pi))

	#Deciding which months to average over
	av_months=''
	if args.av_months=='jfm':
		av_months=np.array([0,1,2])
	if args.av_months=='jas':
		av_months=np.array([6,7,8])
	if args.av_months=='all':
		av_months=range(12)
	if args.av_months=='max':
		av_months='max'
	if args.av_months=='max':
		av_months='max'

        #Creating name for output file:
        output_file = 'figures/' + 'time_series_' + args.field + '_' + np.str(start_year) + '_to_' + np.str(end_year) + '_av_months_' + args.av_months + '.png'
	print output_file

	#Defining a list of colors
	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	#color_vec=np.array(['red', 'red', 'red', 'black','black', 'black', 'black', 'black', 'black', 'black', 'black', 'black',  'black', 'orange', 'coral', 'yellow', 'orchid' ])

	#root_path='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg'	
	root_path='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg'		
	if args.exp_name=='Gold1' or args.exp_name=='Gold2':
		                        root_path='/ptmp/aas/model_output/ulm_201510_mom6_2015.10.23/MOM6_GOLD_SIS2'
	#run_names=np.array(['Delta1', 'Delta2', 'Delta3','Delta4','Delta5','Delta6','Delta7','Delta8', 'Delta9','Delta10','Delta11']) ; naming_flag='Delta'
	run_names=np.array(['tournadre_C3']) ; naming_flag='Tournadre'
	#run_names=np.array(['tournadre', 'rolling_tournadre_Burton_new','rolling_tournadre_WM' ]) ; naming_flag='Tournadre'
	#run_names=np.array(['Control','freq', 'Delta1', 'Delta10','Control_AB','freq_AB', 'Delta1_AB', 'Delta10_AB'])
	
	label_list=run_names
	#label_list=np.array(['Martin and Adcroft', 'Burton','Weeks and Mellor' ]) 

	#root_path='/ptmp/aas/model_output/ulm_201510_mom6_2015.10.23/MOM6_GOLD_SIS2'
	#run_names=np.array(['none','freq','Delta1','Delta_10','none_AB','freq_AB','Delta1_AB', 'Delta10_AB']) ; naming_flag='Gold'
	#run_names=np.array(['Delta10','none_AB','freq_AB','Delta1_AB', 'Delta10_AB']) ; naming_flag='Gold'

	####################################################################################
	########################### Starting the loops #####################################
	####################################################################################
	
	mean_data=[]
	xdata=[]
	xdata2=[]
	for k in range(1):	
	#for k in np.array([5]):	
		#input_folder=all_paths[k]
		input_folder=root_path + '_bergs_' + run_names[k]
		if args.exp_name=='Control':
			input_folder=root_path
		print 'Input folder: ', input_folder
		
		#Set max and min seperately each time
		(min_year, max_year)=find_max_min_years(input_folder,args.file_type+'_month.nc')
		max_year=2049
		start_year=min(start_year, min_year)
		end_year=max_year
		print 'Start year: ' , start_year, ' , End year: ' , end_year

		#Fetching the time series data:
		(data1, time1)=generate_time_series_data(start_year, end_year , input_folder, args.field,pole, args.file_type,av_months,lat_lon_bounds)
		label= label_list[k]
		
		#Old labeling system
		#label= str.split(input_folder,'_')[-1]
		#if len(label)>10 or label[0]=='1':
		#	label='Control'
		print label
	
		#Finding averaged values:
		if plot_averaged_values==1:
			mean_data.append(np.mean((data1)))
			xdata.append(initial_radius_vec[k])
			xdata2.append(initial_mass_vec[k])
			#xdata2.append(initial_radius_vec[k]**(0.001))
			##xdata=xdata.append((initial_mass_vec[k]))
			#mass_estimate=mass_estimate.append(Q*sqrt((initial_mass_vec[k]/(rho_i*np.pi)))

		if plot_time_series==1:
			plt.plot(time1, data1 , color=color_vec[k], linewidth=2.5, linestyle="-", label=label)

			if save_matfile==True:
				save_timeseries_matfile(time1, data1,label,args.field,args.lon_sector)

	#This is when you are plotting an averaged quanitiy from each time series
	if plot_averaged_values==1:
		#subplot(2,1,2)
		plt.plot(xdata,mean_data, linewidth=3.5, linestyle="-")
		plt.plot(xdata,mean_data, 'ro')
		plt.xlabel('Initial Radius')
		plt.ylabel('Total iceberg mass (Gt)')

		#subplot(2,1,1)
		#plt.plot(xdata2,mean_data, linewidth=3.5, linestyle="-")
		#plt.plot(xdata2,mean_data, 'ro')
		#plt.xlabel('Initial Mass')
		#plt.ylabel('Total iceberg mass (Gt)')
		#plt.title('Total iceberg mass')


	if plot_time_series==1:
		plt.xlabel('Time (years)')
		plt.ylabel(args.field)
		#plt.ylabel('Total iceberg mass (Gt)')
		#plt.title('Total iceberg mass')
		plt.legend(loc='upper left', frameon=True)	

	fig = matplotlib.pyplot.gcf()
	fig.set_size_inches(9,4.5)
	fig.set_size_inches(9,4.5)
	if save_figure==True:
		plt.savefig(output_file, dpi=150, bbox_inches='tight', pad_inches=0.4)
	plt.show()

	print 'Script complete'

if __name__ == '__main__':
    sys.exit(main())


