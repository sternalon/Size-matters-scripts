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



def load_time_series_data_from_mat_file(filename):
        mat_contents=sc.loadmat(filename)
	time=mat_contents['time'][0]
	data=mat_contents['data'][0,]
	Delta_name=mat_contents['Delta_name']

	return [time, data, Delta_name]

def produce_figure_panel(datapath,file_extension,field_name, Delta_list,nrows,ncols,plot_num,color_vec,y_label_list,Area_list):
	print 'Plotting ' + field_name + ' time_series'
	subplot(nrows,ncols,plot_num-3)
	#for k in range(10):
	for k in np.array([9,8,7,6,5,4,3,2,1,0]):
		filename=datapath + Delta_list[k] + file_extension
		(time, data, Delta_name)=load_time_series_data_from_mat_file(filename)
		plt.plot(time, data , color=color_vec[k], linewidth=2.5, linestyle="-", label=Area_list[k])
		#plt.ylabel(y_title_list[row_num-1], fontsize=18, labelpad=25)
		
	plt.xlabel('Time (years)')
	plt.ylabel(y_label_list[field_name])
	plt.xlim(0, 100)


def take_mean_over_end_of_time_series(data,time,number_of_years):
	max_time=np.max(time)
	min_year=max_time-number_of_years
	mean_value=(np.sum((time>min_year)*data))/np.sum(time>min_year)
	

	return mean_value


def produce_bar_figure_panel(datapath,file_extension,field_name, Delta_list,nrows,ncols,plot_num,color_vec,y_label_list,title_list,Gladstone_A,scale_by_control,ax):
	print 'Plotting ' + field_name + ' time_series'
	scale_by_control=False
	Num_runs=len(Delta_list)
	mean_list=np.zeros([Num_runs,4])
	ind=np.arange(Num_runs)
	x_tick_list=[]
	
	fetch_list=ind

	color_vec=np.array(['blue','red','green','black'])
	#ind2=ind[:]+width/2.
	width=.35
	number_of_years=60
	sector_list=np.array(['', '_AB_plus_Ross','_EastA_plus_AdeleLand','_Weddell'])
	sector_title_list=np.array(['All', 'AB + Ross Seas','East Antarctica','Weddell Sea'])

	#N=Num_runs
 	if field_name=='mass':
		scale_by_control=False
		#Num_runs=
		fetch_list=np.arange(1,Num_runs)
		x_tick_list.append('')
		#print fetch_list

	if field_name=='melt':
		control_extension='_timeseries_' + 'CALVING'
	else:
		control_extension=file_extension

	for j in range(4):
	#for j in np.array([1,2,3]):
		if scale_by_control==True:
				filename=datapath + 'Control' + control_extension + sector_list[0]+'.mat'
				(time, data, Delta_name)=load_time_series_data_from_mat_file(filename)
				Control_value=np.mean(data)
				Control_value=take_mean_over_end_of_time_series(data,time,number_of_years)
		#for k in range(Num_runs):
		for k in fetch_list:
			filename=datapath + Delta_list[k] + file_extension + sector_list[j]+'.mat'
			if Delta_list[k]=='Control' and field_name=='melt':
				filename=datapath + 'Control' + control_extension + sector_list[j]+'.mat'
			(time, data, Delta_name)=load_time_series_data_from_mat_file(filename)
			#mean_list[k,j]=np.mean(data)
			mean_list[k,j]=take_mean_over_end_of_time_series(data,time,number_of_years)
			#if Delta_list[k]=='Delta8' or Delta_list[k]=='Delta7' or Delta_list[k]=='Delta10':
			#	mean_list[k,j]=0
			x_tick_list.append(title_list[Delta_list[k]])
		
		if scale_by_control==True:
			mean_list[:,j]=mean_list[:,j]/Control_value
			#mean_list[:,j]=100*(mean_list[:,j]-1) #Make it a percent change
	

		#Unit conversion
		if field_name=='mass':
			mean_list=mean_list/(10**12)  #Convert to Gton
		if field_name=='melt' or field_name=='CALVING':
			mean_list=(60*60*24*365.25)*mean_list/(10**12) #Convert to Gton per year

		scale_by_total_anatarctic=False
		if scale_by_total_anatarctic==True:
			if j==0:
				Full_Ant=mean_list[:,0]
			if j>0:
				mean_list[:,j]=np.squeeze(mean_list[:,j])/Full_Ant
				print 'Second time:', mean_list[:,j]



		if field_name=='mass':
			ax.bar(ind+((width/2)*j), mean_list[:,j],width=width/2,color=color_vec[j],label=sector_title_list[j])
			ax.set_xticks(ind + width/2)
			plt.xlim([1, Num_runs]) #For area
			plt.xlabel('Calving Length (m)', fontsize=18)
		else:
			ax.bar(ind+((width/2)*j)+(3*width*(ind>0)), mean_list[:,j],width=width/2,color=color_vec[j])
			ax.set_xticks(ind + width/2  +(3*width*(ind>0)))
			plt.xlim([0.0, Num_runs+1]) #For area
			plt.xlabel('       Calving Length (m)', fontsize=18)
			ax.text(2, -800, '                    Iceberg Melt ( F$_{IM}$ )                            ', style='normal',
					        bbox={'facecolor':'grey', 'alpha':0.5, 'pad':20},fontsize=16)
			ax.text(-0.5, -800, ' Antarctic \n Freshwater \n Runoff (F$_{R}$ ) \n (all simulations)', style='normal',
					        bbox={'facecolor':'grey', 'alpha':0.5, 'pad':20},fontsize=13)

		if j<10:
			print 'Working out budgets for the ' , sector_list[j]
			print mean_list[:,j]
			print ' '
			mass_diff=(mean_list[:,j] - mean_list[0,j])*900
			L_ice=334*1000 #Joule/kg
			Year_sec=365*24*60*60
			SH_sea_ice_area = 10**13
			Watt_diff=mass_diff*L_ice/SH_sea_ice_area/Year_sec
			percent_diff=((mean_list[:,j] / mean_list[0,j]) -1)*100

			print 'Mass in control', mean_list[:,0]*900
			print ' '
			print 'Mass differernce', mass_diff
			print ' '
			print 'Percent differernce', percent_diff
			print ' '
			print 'Watt differernce', Watt_diff
			print ' '



		#plt.ylim([0.1, 1.3]) #For area
		#ax.set_yscale('log')
		
		
		ax.set_xticklabels(x_tick_list, rotation=70)
		plt.ylabel(y_label_list[field_name], fontsize=18)
		#plt.xlim(0, 100)

		#ax.bar(ind+((width/2)*4), np.sum(mean_list[:,1:4],axis=1),width=width/2,color='orange')
		if field_name=='mass':
			#plt.legend(loc=(0.8,0.93), frameon=True,prop={'size':14})
			plt.legend(loc=2, frameon=True,prop={'size':14})


##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################

def main():

	parser = argparse.ArgumentParser()
   	args = parser.parse_args()
	
	#Flags
	save_figure=True
	plot_time_series=False
	scale_by_control=True

	#Parameters;
	pole='south'
	title=None
	cscale=None
	#Delta_list=np.array(['Delta1', 'Delta2', 'Delta3','Delta4', 'Delta5', 'Delta6','Delta7', 'Delta8', 'Delta9','Delta10' ])
	#Delta_list=np.array(['Delta1', 'Delta2', 'Delta3','Delta4', 'Delta5', 'Delta6','Delta7', 'Delta8', 'Delta9','Delta10','Control' ])
	#Delta_list=np.array(['Control','Delta1', 'Delta2', 'Delta3','Delta4', 'Delta5', 'Delta6','Delta7', 'Delta8', 'Delta9','Delta10' ])
	Delta_list=np.array(['Control','Delta1', 'Delta2', 'Delta3','Delta4', 'Delta5', 'Delta6', 'Delta9' ])
	#Delta_list=np.array(['Control','Delta1', 'Delta2', 'Delta3','Delta4', 'Delta5', 'Delta6', 'tournadre' ])
	Area_list=np.array(['Area$_{1}$=0.0026km$^2$', 'Area$_{2}$=0.0072km$^2$','Area$_{3}$=0.029km$^2$',\
			'Area$_{4}$=0.12km$^2$','Area$_{5}$=0.18km$^2$', 'Area$_{6}$=0.35km$^2$','Area$_{7}$=0.56km$^2$', 'Area$_{8}$=1.0km$^2$','Area$_{9}$=1.8km$^2$', 'Area$_{10}$=3.5km$^2$'])
	Gladstone_A=np.array([0.0026, 0.0072, 0.0292, 0.1210, 0.1788, 0.3529, 0.5647, 1.0353, 1.8353, 3.4824])  ; Gladstone_A=Gladstone_A*(10**(6))  #Convert to m^2
	letter_labels=np.array(['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'])
	
	Area_list=np.array(['F_{R}'])
	nrows=1
	ncols=2#len(Delta_list)#3
	
	number = 10
	cmap = plt.get_cmap('hot')
	colors = [cmap(i) for i in np.linspace(0, 0.5, number)]
	#color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	#color_vec=np.array(['red', 'red', 'red', 'black','black', 'black', 'black', 'black', 'black', 'black', 'black', 'black',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	color_vec=colors
	
	datapath='/home/Alon.Stern/Iceberg_Project/iceberg_scripts/python_scripts/size_matters_paper/processed_data/'
        
	#field_list=np.array(['mass','CN','HI','Volume'])
	#field_list=np.array(['mass','melt','CALVING','Volume'])
	field_list=np.array(['mass','melt'])
	field_list=np.array(['Volume','Volume'])
	#field_list=np.array(['mass','calving','CN_AB','HI_AB'])
	#nrows=len(field_list)
	
	#extension_list = {'mass': '_timeseries_mass.mat', 'calving':'_timeseries_calving.mat', 'CN_AB':'_timeseries_CN_AB.mat' 'melt':'_timeseries_melt.mat'}
	y_label_list = {'mass': 'Iceberg Mass (Gton)', 'CALVING':'Calving Flux (Gton/year)', 'CN':'Sea Ice Concentration (m$^2$)' ,'HI': 'Sea Ice Thickness (m$^3$)',\
			'melt':'Freshwater Flux (Gton/year)', 'Volume':'Volume (m$^{3}$'}
	#y_label_list = {'mass': 'Iceberg Mass (Gton)', 'CALVING':'Calving Flux Change (%)', 'CN':'Sea Ice Concentration Change (%)' ,'HI': 'Sea Ice Thickness Change (%)',\
	#		'melt':'Melt Change(%)', 'Volume':'Ice Volume Change (%)'}


	title_list = {'Control':'F$_{R}$', 'Delta1': 'L=60m', 'Delta2':'L=100m', 'Delta3':'L=200m' , 'Delta4': 'L=350m' , 'Delta5': 'L=500m','Delta6': 'L=700m',\
			'Delta7':'L=900m', 'Delta8':'L=1200m$' , 'Delta9': 'L=1600m' , 'Delta10': 'L=2200m', 'tournadre': 'L=10000m'}

	#Setting up the figure
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
	#Plotting Melt Fields
	for plot_num in range(2):
	#for plot_num in range(nrows*ncols):
		field_name=field_list[plot_num]
		#file_extension=extension_list[field_name]
		file_extension='_timeseries_' + field_name
		if plot_time_series==True:
			produce_figure_panel(datapath,file_extension,field_name, Delta_list,nrows,ncols,plot_num,color_vec,y_label_list,Area_list)
		else:
			ax=plt.subplot(nrows,ncols,plot_num+1)
			produce_bar_figure_panel(datapath,file_extension,field_name, Delta_list,nrows,ncols,plot_num,color_vec,y_label_list,title_list,Gladstone_A,scale_by_control,ax)
		
		text(1,1,letter_labels[plot_num], ha='right', va='bottom',transform=ax.transAxes)
	fig.subplots_adjust(bottom=0.3)

	#plt.legend(loc='lower_right', frameon=True)
	#plt.legend(loc=(0.7,-0.25 ), frameon=True,prop={'size':12})
	#plt.tight_layout()



	fig.set_size_inches(20.5, 10.5,forward=True)
	if save_figure==True:
		plt.savefig('Ice_volume_calculation.png')
	plt.show()


if __name__ == '__main__':
	main()
	#sys.exit(main())

