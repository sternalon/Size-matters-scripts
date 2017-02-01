#!/usr/bin/env python

import sys
import argparse
import netCDF4 as nc
import numpy as np
import scipy.io as sc
from scipy import stats
import matplotlib
import matplotlib as mpl
from pylab import *
from matplotlib import ticker
#matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from sea_ice_concentrations import plot_polar_field
#from sea_ice_concentrations import define_longitudes
from sea_ice_concentrations import create_output_file_and_save


def dunne_pm(N=256):
	"""
	Plus-minus  colormap from John Dunne.
	"""
	cdict = {'red':   [(0.00, 0.3, 0.3),
		(0.05, 0.5, 0.5),
		(0.20, 0.0, 0.0),
		(0.30, 0.4, 0.4),
		(0.40, 0.8, 0.8),
		(0.50, 1.0, 1.0),
		(0.95, 0.6, 0.6),
		(1.00, 0.4, 0.4)],

		'green': [(0.00, 0.0, 0.0),
			(0.30, 0.5, 0.5),
			(0.40, 1.0, 1.0),
			(0.70, 1.0, 1.0),
			(1.00, 0.0, 0.0)],

		'blue':  [(0.00, 0.3, 0.3),
			(0.05, 0.5, 0.5),
			(0.20, 1.0, 1.0),
			(0.50, 1.0, 1.0),
			(0.60, 0.7, 0.7),
			(0.70, 0.0, 0.0),
			(1.00, 0.0, 0.0)]}
	import matplotlib
	cmap = matplotlib.colors.LinearSegmentedColormap('dunnePM', cdict, N=N)
	cmap.set_under([.1,.0,.1]); cmap.set_over([.2,0.,0.])
	#cmap.set_bad('w')
	matplotlib.cm.register_cmap(cmap=cmap)
	return cmap



def load_data_from_mat_file(filename,significant_vals_only):
        mat_contents=sc.loadmat(filename)
	lat=mat_contents['lat'][0,:]
	lon=mat_contents['lon'][0,:]
	data=mat_contents['data'][:,:]
	field=mat_contents['field']
	title=mat_contents['title']
	p_values=mat_contents['p_values']
	#if p_values==0:
	if significant_vals_only is False:
		p_values=None
	return [lat, lon, data,p_values,field ]


def load_mass_vs_age_from_mat_file(filename,load_type='lognormal'):

        mat_contents=sc.loadmat(filename)
        Total_x=mat_contents['Total_x'][:]

        if load_type=='occurance':
                occurance_x=np.squeeze(mat_contents['occurance_x'][:])
                occurance_y=np.squeeze(mat_contents['occurance_y'][:])
                occurance_matrix=mat_contents['occurance_matrix'][:]
                print occurance_x.shape
                print occurance_y.shape
                print occurance_matrix.shape
                return [ occurance_x, occurance_y, occurance_matrix]
        if load_type=='median':
                Total_median=mat_contents['Total_median'][:]
                Total_mean=mat_contents['Total_mean'][:]
                median_matrix=mat_contents['median_matrix'][:]
                mean_matrix=mat_contents['mean_matrix'][:]
                x_matrix=mat_contents['median_matrix'][:]
                return [Total_median, median_matrix, Total_mean, mean_matrix, Total_x, x_matrix]

def smooth_the_line(field):
	N=len(field)
	print N
	R= 8
	new_field=field.copy()

	S=30  #Starting point
	for k in range(R-1,N-R):
		if k <N-S:
			new_field[k]= np.sum(field[np.arange(k-R,k+R+1)]) /(2*float(R)+1 )
		else:
			new_field[k]=field[k]
	return new_field


		



##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################

def main():

	parser = argparse.ArgumentParser()
   	args = parser.parse_args()
	
	
	#Flags
	save_figure=True
	significant_vals_only=False
	plot_melt_fields=False
	plot_occurance_matrix=True
	

	#Parameters;
	pole='south'
	title=None
	cscale=None
	Delta_list=np.array(['Rolling','Rolling'])
	Pole_list=np.array(['south','Greenland2'])
	#Title_list=np.array(['Delta1 (Area=0.0026km^2)', 'Delta6 (Area=0.029km^2)', 'Delta9 (Area=1.8km^2)'])
	#Title_list=np.array(['Area=0.0026km^2', 'Area=0.35km^2', 'Area=1.8km^2'])
	#title_list = {'Delta1': 'Area=0.0026km^2', 'Delta3':'Area=0.029km^2', 'Delta4':'Area=0.12km^2' , 'Delta6': 'Area=0.35km^2' , 'Delta9': 'Area=1.8km^2'}
	#title_list = {'Delta1': 'Length=60m', 'Delta3':'Length=200m', 'Delta4':'Length=350m$' , 'Delta6': 'Length=700m$' , 'Delta9': 'Length=1600m'}
	title_list = {'Delta1': 'Length=62m', 'Delta3':'Length=209m', 'Delta9': 'Length=1659m'}
	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	letter_labels=np.array(['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'])
	y_title_list=np.array(['','','','','',''])
	#colorbar_unit_list=np.array(['', 'melt (kg/m^2/s)', 'Conc (non-dim)', 'Thickness (m)',''])
	colorbar_unit_list=np.array([ 'melt (m/yr)', 'melt (m/yr)', 'Thickness (m)',''])
	nrows=1
	ncols=2
	boundinglat=-45
	projection = 'splaea'
	Name_list=np.array(['Rolling', 'No Rolling'])


        #Rolling_file='processed_data/rolling_tournadre_Burton_new_age_vs_mass_1946_to_2046_expected_values_logaxis__fixed_upperboundNumbins100.mat'
	#No_Rolling_file='processed_data/rolling_tournadre_Rolling_off_age_vs_mass_1943_to_2043_expected_values_logaxis__fixed_upperboundNumbins100.mat'
        #Southern_Hemisphere
        #No_Rolling_file='processed_data/rolling_tournadre_Rolling_off_age_vs_mass_1980_to_2030_expected_values_constraints_lat_50.0_to_75.0_logaxis__fixed_upperboundNumbins100.mat'
        #Rolling_file='processed_data/rolling_tournadre_Burton_new_age_vs_mass_1980_to_2030_expected_values_constraints_lat_50.0_to_75.0_logaxis__fixed_upperboundNumbins100.mat'
	Rolling_file='processed_data/rolling_tournadre_Burton_new_age_vs_mass_1950_to_2000_expected_values_logaxis__fixed_upperboundNumbins100.mat'
	No_Rolling_file='processed_data/rolling_tournadre_Rolling_off_age_vs_mass_1950_to_2000_expected_values_logaxis__fixed_upperboundNumbins100.mat'


	load_filename_list=np.array([Rolling_file,No_Rolling_file])

	
	datapath='/home/Alon.Stern/Iceberg_Project/iceberg_scripts/python_scripts/size_matters_paper/processed_data/'
        
	#Setting up the figure
	#fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	fig = plt.figure(1)
        ax = fig.add_subplot(1,1,1)

	
	if plot_occurance_matrix is True:
                for k in np.array([0,1]):
			mass_num=5
                        filename=load_filename_list[k]
                        (occurance_x, occurance_y, occurance_matrix) =load_mass_vs_age_from_mat_file(filename,load_type='occurance')
                        dp=occurance_x[1]-occurance_x[0]
                        max_occurance=np.zeros((len(occurance_y)))
                        std_occurance=np.zeros((len(occurance_y)))
                        for j in range(len(occurance_y)):
                                #max_occurance[j]=np.sum(  occurance_matrix[mass_num,:,j]*occurance_x*dp /np.sum(occurance_matrix[mass_num,:,j])  )
                                max_occurance[j]=np.sum(  occurance_matrix[mass_num,:,j]*occurance_x /np.sum(occurance_matrix[mass_num,:,j])  )
				std_occurance[j]=sqrt(np.sum(  occurance_matrix[mass_num,:,j]*((occurance_x-max_occurance[j])**2) /np.sum(occurance_matrix[mass_num,:,j])  ))
                        vmax=np.max(occurance_matrix[mass_num,:,:])
                        #cNorm = MidpointNormalize(vmin=0, vmax=5,midpoint=0)
                        #cNorm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)
                        cNorm = mpl.colors.LogNorm(vmin=10**-5, vmax=0.1)
                        #if k==0:
                        #       max_occurance0=max_occurance
			max_occurance=smooth_the_line(max_occurance)
			std_occurance=smooth_the_line(std_occurance)
                        if k==0:
                                plt.pcolor(occurance_x, occurance_y, (np.squeeze(occurance_matrix[mass_num,:,:])).transpose(),cmap='Blues', norm=cNorm)
                                cbar=plt.colorbar()
				cbar.set_label('Probability', rotation=90,fontsize=20)
                                #plt.title(Name_list[k])
                        #plt.plot(max_occurance,occurance_y,color=color_vec[k],linewidth=2)
                        plt.plot(max_occurance,occurance_y,color=color_vec[k],linewidth=4, label=Name_list[k])
                       	if k==0: 
				plt.plot(max_occurance+std_occurance,occurance_y,  color=color_vec[k],linewidth=2,linestyle=':')
				plt.plot(max_occurance-std_occurance,occurance_y,  color=color_vec[k],linewidth=2,linestyle=':')
                        ax.set_yscale('log')
		plt.xlim([0,25])
		#plt.ylim([10**6.7,7.5e10])
		plt.ylim([10**6.0,7.5e10])
		plt.xlabel('Age (years)')
		plt.ylabel('Mass (kg)')
		plt.legend(loc='upper right', frameon=True)

	

	subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
	fig.set_size_inches(12.0, 10.5,forward=True)
	if save_figure==True:
		plt.savefig('Till_paper_figures/Till_Fig_lifetime.png',dpi=300,bbox_inches='tight')
	plt.show()


if __name__ == '__main__':
	main()
	#sys.exit(main())

