#!/usr/bin/env python
from netCDF4 import Dataset
import numpy as np
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
from sea_ice_concentrations import plot_polar_field
from sea_ice_concentrations import generate_data_file
from distributions_of_bergs import get_valid_data_from_traj_file
from mpl_toolkits.basemap import Basemap
from iceberg_mass_comparisons import find_max_min_years
from distributions_of_bergs import add_constraint_to_data



###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
def main():
	#Clear screen
	#os.system('clear')


	input_folder='/ptmp/Alon.Stern/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta9/'
	field='CN'
	pole='north'
	pole=None
	save_traj_data=False

	input_folder='/home/Alon.Stern/Iceberg_Project/iceberg_scripts/python_scripts/size_matters_paper/Delta9_CN_1959_to_2019_all_south_xy_anomaly.mat'
	(All_data, lat, lon,z) =generate_data_file(start_year, start_year, input_folder, field, months_str, file_type)
	data_mean=np.squeeze(np.mean(All_data,axis=0))
	data_mean=None

	#Deciding on the projection.
	if pole=='south':
		projection = 'splaea'
		boundinglat=-45
	if pole=='north':
		projection = 'nplaea'
		boundinglat=45
	if pole==None:
		projection = 'lcc'
		lat_1=80
		lat_2=60
		lat_0=63.5
		lon_0=-59.


	run_names=np.array(['tournadre']) ; naming_flag='Delta'
	#run_names=np.array(['Delta1', 'Delta3', 'Delta9']) ; naming_flag='Delta'
	#run_names=np.array(['Delta1', 'Delta2','Delta3','Delta6', 'Delta9','Delta10']) ; naming_flag='Delta'
	#run_names=np.array(['Delta1', 'Delta2','Delta3', 'Delta9']) ; naming_flag='Delta'

		datamap=m.scatter(Total_berg_lon,Total_berg_lat, c=Total_berg_mass0, marker='o',cmap='jet',norm=cNorm)


	#output_file= 'figures/Leigh_Rink_figure2.png'
	#plt.savefig(output_file, dpi=150, bbox_inches='tight', pad_inches=0.4)

	plt.show()



	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

