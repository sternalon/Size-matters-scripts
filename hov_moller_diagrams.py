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
from pylab import *
#matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from iceberg_mass_comparisons import find_max_min_years
from iceberg_mass_comparisons import define_mask
from m6toolbox  import section2quadmesh
from m6toolbox  import rho_Wright97

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


	(control_mean, cscale)=produce_hov_moller_diagram(start_year, end_year,input_folder,args.field, second_folder, args.months,args.pole,\

			args.file_type,axes_direction,lat_lon_bounds,cscale,control_mean,significant_vals_only,\
			
			plot_lon_sum=False,difference_color_off=False,title=None,second_file_type=second_file_type,save_mat_file=save_mat_file,lon_sector=args.lon_sector)



	plt.show()


if __name__ == '__main__':
	    sys.exit(main())

