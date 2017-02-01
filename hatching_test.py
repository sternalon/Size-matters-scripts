#!/usr/bin/env python
from netCDF4 import Dataset
import numpy as np
import netCDF4 as nc
from PIL import *
from pylab import *
import math
import os
import matplotlib.pyplot as plt
import pdb
import argparse
import sys
from scipy.interpolate import interp1d
from iceberg_mass_comparisons import define_paths_array
from iceberg_mass_comparisons import find_max_min_years
from distributions_of_bergs import get_valid_data_from_traj_file
from distributions_of_bergs import add_constraint_to_data
from sea_ice_concentrations import plot_polar_field
from sea_ice_concentrations import generate_data_file


###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
	
	
	
def main():
	start_year=1947
        input_folder='/ptmp/Alon.Stern/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta2/'
	field='CN'
	months=range(0,11)
	file_type='ice_month'
	pole='south'
	(All_data, lat, lon) =generate_data_file(start_year, start_year, input_folder, field,months,file_type)
	data=np.squeeze(All_data)
	print data.shape
	#plt.show()

	#data(volcano)
	volcano=data
	m <- volcano
	dimx <- nrow(m)
	dimy <- ncol(m)

	d1 <- list(x = seq(0, 1, length = dimx), y = seq(0, 1, length = dimy), z = m)



	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

