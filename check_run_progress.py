#!/usr/bin/env python

import sys
import os
import argparse
import netCDF4 as nc
import numpy as np
import pandas as pd
import re
import datetime
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
from iceberg_mass_comparisons import find_max_min_years

##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################


parser = argparse.ArgumentParser()
parser.add_argument('--dirtype', default='history',help='Histroy vs Restart files.')
parser.add_argument('--source', default='archive',help='archive vs ptmp.')
args = parser.parse_args()

print	'**************************************************************************************'
print	'    ********************* Listing ' + args.dirtype + ' files *********************************'
print	'**************************************************************************************'
#Defining possible paths
all_paths=[]
if args.source=='archive':
	root_path='/archive/aas/ulm_mom6_2015.07.20_again_myIcebergTag/'
if args.source=='ptmp':
	root_path='/ptmp/Alon.Stern/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/'
extend_path='/gfdl.ncrc2-intel-prod-openmp/'
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_mass' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_freq' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_all_small' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_all_big' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta1' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta2' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta3' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta4' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta5' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta6' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta7' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta8' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta9' )
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta10')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta11')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta12')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta13')

all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta1b')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta9b')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta10b')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta11b')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta12b')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta13b')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta14b')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta1_passive')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta9_passive')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_tournadre')

all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta1_thick')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta10_thick')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta1b_thick')
all_paths.append('AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta10b_thick')

#	all_paths=define_paths_array()
for k in range(32):
	if args.source=='archive':
		input_folder=root_path+all_paths[k] + extend_path + args.dirtype +'/'
		(min_year, max_year)=find_max_min_years(input_folder,'.tar')
	if args.source=='ptmp':
		input_folder=root_path+all_paths[k] +'/'
		(min_year, max_year)=find_max_min_years(input_folder,'.nc')
	print 'Max year in ' +all_paths[k] + ' is:  ' + str(max_year)
print 'Script complete'

################################################################################################
################################# Ocean Only Runs ##############################################
################################################################################################

print	'**************************************************************************************'
print	'    ********************* Listing ' + args.dirtype + ' files *********************************'
print	'**************************************************************************************'
#Defining possible paths
all_paths=[]
if args.source=='archive':
	root_path='/archive/aas/ulm_201510_mom6_2015.10.23/'
if args.source=='ptmp':
	root_path='/ptmp/Alon.Stern/model_output/ulm_201510_mom6_2015.10.23/'
extend_path='/gfdl.ncrc2-intel-prod/'
all_paths.append('MOM6_GOLD_SIS2_bergs_freq')
all_paths.append('MOM6_GOLD_SIS2_bergs_none')
all_paths.append('MOM6_GOLD_SIS2_bergs_Delta1')
all_paths.append('MOM6_GOLD_SIS2_bergs_Delta10')
all_paths.append('MOM6_GOLD_SIS2_bergs_freq_AB')
all_paths.append('MOM6_GOLD_SIS2_bergs_none_AB')
all_paths.append('MOM6_GOLD_SIS2_bergs_Delta1_AB')
all_paths.append('MOM6_GOLD_SIS2_bergs_Delta10_AB')
all_paths.append('MOM6_GOLD_SIS2_bergs_Delta9')

#	all_paths=define_paths_array()
for k in range(8):
	if args.source=='archive':
		input_folder=root_path+all_paths[k] + extend_path + args.dirtype +'/'
		(min_year, max_year)=find_max_min_years(input_folder,'.tar')
	if args.source=='ptmp':
		input_folder=root_path+all_paths[k] +'/'
		(min_year, max_year)=find_max_min_years(input_folder,'.nc')
	print 'Max year in ' +all_paths[k] + ' is:  ' + str(max_year)
print 'Script complete'
