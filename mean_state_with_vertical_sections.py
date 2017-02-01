#!/usr/bin/env python

import sys
import argparse
import netCDF4 as nc
import numpy as np
import scipy as sc
from scipy import stats
import scipy.interpolate
import datetime
import matplotlib
import matplotlib as mpl
from pylab import *
#matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from iceberg_mass_comparisons import find_max_min_years
from iceberg_mass_comparisons import define_mask
from sea_ice_concentrations import produce_map_displaying_data
from sea_ice_concentrations import define_longitudes
from sea_ice_concentrations import create_output_file_and_save

##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument('--input_folder', default=None, help='The input data file in NetCDF format.')
	parser.add_argument('--second_folder', help='The input data file in NetCDF format.')
	parser.add_argument('--field',default='salt', help='The field to plot.')
	#parser.add_argument('--output_file',default = 'output_seaice_conc.png', help='The output figure name .png')
	parser.add_argument('--output_file', help='The output figure name .png')
    	parser.add_argument('--pole', default='south', help="""Plot South, North or both poles.""")
	parser.add_argument('--months', default='all', help="""Which models is the averaging done over: all, djf, jja.""")
	parser.add_argument('--file_type', default='ocean_month_z',help='The fields included in this plot.')
	parser.add_argument('--cscale', default=None,help='Scale for the axes in the plot.')
	parser.add_argument('--berg_type', default=None,help='which berg experiment are you using (eg: Delta1.')
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
	naming_flag='Mean_state' ;
	
	#If no imput file is used, then plot multiple panels.
	if args.input_folder==None:
		plot_multiple_panels=1

	#Parameters;
	#start_year=30#1925
	#end_year=39

	start_year=1950#1925
	end_year=2019
	Number_of_years=60
	axes_direction='yz'
	#axes_direction='xy'
	lat_min=-85 ;lat_max=-40

	if args.file_type=='ice_month':
		axes_direction='xy'

	if plot_multiple_panels==1:
		input_folder='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg'
		#input_folder='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta9'
		second_folder=None
		if args.berg_type!=None: 
			second_folder=input_folder
			input_folder='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg' + '_bergs_' + args.berg_type
		print 'input_folder: ', input_folder
		print 'second_folder: ', second_folder

		#run_type=np.array(['surface', 'AB', 'Ross','AdeleLand', 'EastA','Weddell'])  
		#run_type=np.array(['AB_plus_Ross', 'Weddell', 'AB', 'Ross','AdeleLand', 'EastA'])  
		#run_type=np.array(['Ross','AdeleLand', 'EastA'])  
		#run_type=np.array(['120West', '52West','60East'])  
		#run_type=np.array(['surface','AB_plus_Ross', 'Weddell','EastA_plus_AdeleLand'])  
		run_type=np.array(['AB_plus_Ross', 'Weddell','EastA_plus_AdeleLand'])  
		#run_type=np.array(['AB_plus_Ross'])  
		#run_type=np.array(['Weddell','EastA_plus_AdeleLand'])  
		#run_type=np.array(['AB_plus_Ross', 'EastA_plus_AdeleLand'])  
		#run_type=np.array(['EastA_plus_AdeleLand'])  
		if args.berg_type!=None:
			naming_flag=naming_flag +'_'+args.berg_type
		if significant_vals_only==True:
			naming_flag=naming_flag+'_sig'

		(min_year, max_year)=find_max_min_years(input_folder,'.' + args.file_type + '.nc')
		if second_folder!=None:
			(min_year2, max_year2)=find_max_min_years(second_folder,'.' + args.file_type + '.nc')
			min_year=max(min_year,min_year2) ; max_year=min(max_year,max_year2)  ;
		end_year=max_year
		start_year=max(end_year-Number_of_years,min_year)
		print 'File: ' + input_folder +' from '+ str(start_year) + ' to ' + str(end_year)
		cscale=None
		for k in range(3):
			print 'Next in loop: ' , k , run_type[k]
			field=args.field
			if run_type[k]=='surface':
				axes_direction='xy'
				file_type='ice_month'
				if args.field=='temp':
					field='SST'
				if args.field=='salt':
					field='SSS'
				if args.field=='u':
					field='UO'
				if args.field=='rho':
					field='rho'
			else:
				axes_direction='yz'
				#file_type='ocean_month_z'
				file_type=args.file_type
				if k==1:
					naming_flag=naming_flag+'_'+file_type
				field=args.field
			lon_sector=run_type[k]
			lat_lon_bounds=define_longitudes(lon_sector)
			lat_lon_bounds[0]=lat_min ;lat_lon_bounds[1]=lat_max ;
			subplot(3,1,k+1)
			control_mean=None # This line makes the contrl reload each time
			#(control_mean, cscale)=produce_map_displaying_data(start_year, end_year,input_folder,args.field, second_folder, args.months,args.pole,args.file_type,cscale,control_mean)
			(control_mean, cscale)=produce_map_displaying_data(start_year, end_year,input_folder,field, second_folder, args.months,args.pole,\
					file_type,axes_direction,lat_lon_bounds,cscale,control_mean,significant_vals_only,plot_lon_sum=False,difference_color_off=False,title=run_type[k]\
					,second_file_type=None,save_mat_file=save_mat_file,lon_sector=lon_sector)
		

	create_output_file_and_save(args.field,start_year,end_year,args.months,args.pole,plot_multiple_panels,args.input_folder,args.second_folder,naming_flag,save_fig)

	plt.show()


if __name__ == '__main__':
	    sys.exit(main())

