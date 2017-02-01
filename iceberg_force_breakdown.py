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



def get_Total_berg_field(input_file,Total_berg_field,field,constraints):
	print 'Loading ' ,field
	get_original=False
	if field[-1]=='0':
		get_original=True
		field=field[:-1]
	all_bergs_field =  get_valid_data_from_traj_file(input_file,field,subtract_orig=False,get_original_values=get_original) # Getting iceberg mass too.
	all_bergs_field=add_constraint_to_data(input_file,all_bergs_field,constraints)
	Total_berg_field=np.concatenate((Total_berg_field,all_bergs_field),axis=0)

	return Total_berg_field

def generate_Total_berg_lon_lat(Total_berg_lon,Total_berg_lat, input_file,constraints,m):
	all_bergs_lat = get_valid_data_from_traj_file(input_file,'lat')
	all_bergs_lon = get_valid_data_from_traj_file(input_file,'lon')
			
	#Adding constraints to the data:
	print 'Lendth of original field: ' , len(all_bergs_lat)
	all_bergs_lat=add_constraint_to_data(input_file,all_bergs_lat,constraints)
	all_bergs_lon=add_constraint_to_data(input_file,all_bergs_lon,constraints)
	
	print 'Lendth of field after constraints: ' , len(all_bergs_lat)
	x, y = m(all_bergs_lon, all_bergs_lat)
	#m.scatter(x,y,3,marker='o',color=color_vec[k])
	
	Total_berg_lon=np.concatenate((Total_berg_lon,x),axis=0)
	Total_berg_lat=np.concatenate((Total_berg_lat,y),axis=0)

	return (Total_berg_lon, Total_berg_lat)



def calculate_forces(Total_berg_list,Total_berg_lat):
	#Parameters
	rho_i=850.
	rho_o=1025.
	rho_a=1.225
	c_av=1.3
	c_ah=0.0055
	c_ov=0.9
	c_oh=0.0012
	c_iv=0.9
	g=9.81
	omega=2*np.pi/(24*60*60)
	
	#Shortening variable names
	W=Total_berg_list['width']
	L=Total_berg_list['length']
	T=Total_berg_list['thickness']
	mass=Total_berg_list['mass']
	D=(rho_i/rho_o)*T
	F=T-D

	u=Total_berg_list['uvel']
	v=Total_berg_list['vvel']
	ua=Total_berg_list['ua']
	va=Total_berg_list['va']
	uo=Total_berg_list['uo']
	vo=Total_berg_list['vo']
	ui=Total_berg_list['ui']
	vi=Total_berg_list['vi']
	hi=Total_berg_list['hi']
	ssh_x=Total_berg_list['ssh_x']
	ssh_y=Total_berg_list['ssh_y']

	f=2*omega*sin(pi*Total_berg_lat/180)


	#Reconstruction the forces:
	#test=sqrt( ((v-va)**2) + ((u-ua)**2))
	#test=(u[:]-ua[:])
	#print u
	#print ua
	#print test
	
	#########

	#Atmospheric Drag
	Fa_x=rho_a*((0.5*c_av*W*F)+(c_ah*L*W))*sqrt( ((v-va)**2) + ((u-ua)**2))*(ua-u)
	Fa_y=rho_a*((0.5*c_av*W*F)+(c_ah*L*W))*sqrt( ((v-va)**2) + ((u-ua)**2))*(va-v)
	Fa_mod=sqrt((Fa_x**2)+(Fa_y**2))
	
	#Oceanic Drag
	Fo_x=rho_o*((0.5*c_ov*W*(D-hi))+(c_oh*L*W))*sqrt( ((v-vo)**2) + ((u-uo)**2))*(uo-u)
	Fo_y=rho_o*((0.5*c_ov*W*(D-hi))+(c_oh*L*W))*sqrt( ((v-vo)**2) + ((u-uo)**2))*(vo-v)
	Fo_mod=sqrt((Fo_x**2)+(Fo_y**2))

	#Sea Ice Drag
	Fi_x=rho_i*(0.5*c_iv*W*hi)*sqrt( ((v-vi)**2) + ((u-ui)**2))*(ui-u)
	Fi_y=rho_i*(0.5*c_iv*W*hi)*sqrt( ((v-vi)**2) + ((u-ui)**2))*(vi-v)
	Fi_mod=sqrt((Fi_x**2)+Fi_y**2)

	#Pressure Gradient Force
	Fp_x=-mass*g*ssh_x
	Fp_y=-mass*g*ssh_y
	Fp_mod=sqrt((Fp_x**2)+(Fp_y**2))
	
	#Coriolis Force
	Fc_x=-mass*f*v
	Fc_y=-mass*f*(-u)
	Fc_mod=sqrt((Fc_x**2)+(Fc_y**2))


	#Wave Radiation Force
	a=0.010125*( ((va-vo)**2) + ((ua-uo)**2))
	L_w=0.32*(((va-vo)**2) + ((ua-uo)**2))
	L_c=0.125*L_w
	L_t=0.25*L_w
	cr=L*0  #Initializing cr
	min_a_F=L*0  #Initializing min_a_F
	for k in range(len(L)):
		cr[k]=0.06*min(max(0,(L[k]-L_c[k])/(L_t[k]-L_c[k])),1)
		min_a_F[k]=min(a[k],F[k])

	Fr_x=0.5*rho_o*cr*g*a*min_a_F*2*(L*W/(L+W+0.0000001))*ua/sqrt(((va)**2) + ((ua)**2))
	Fr_y=0.5*rho_o*cr*g*a*min_a_F*2*L*W/(L+W)*va/sqrt(((va)**2) + ((ua)**2))
	for k in range(len(Fr_x)):
		if sqrt(((va[k])**2) + ((ua[k])**2))==0:
			Fr_x[k]=0
			Fr_y[k]=0

	Fr_mod=sqrt((Fr_x**2)+(Fr_y**2))
	
	F_total_x=Fa_x + Fo_x + Fi_x + Fr_x + Fc_x + Fp_x
	F_total_y=Fa_y + Fo_y + Fi_y + Fr_y + Fc_y + Fp_y
	F_total_mod=Fa_mod + Fo_mod + Fi_mod + Fr_mod + Fc_mod + Fp_mod

	###########



	#Creating dictionaries
	Force_list=np.array(['Fa','Fo','Fi','Fr','Fc','Fp','F_total'])
	Force_x={'Fa':Fa_x,'Fo':Fo_x,'Fi':Fi_x,'Fr':Fr_x,'Fc':Fc_x,'Fp':Fp_x, 'F_total':F_total_x}
	Force_y={'Fa':Fa_y,'Fo':Fo_y,'Fi':Fi_y,'Fr':Fr_y,'Fc':Fc_y,'Fp':Fp_y,'F_total':F_total_y}
	Force_mod={'Fa':Fa_mod,'Fo':Fo_mod,'Fi':Fi_mod,'Fr':Fr_mod,'Fc':Fc_mod,'Fp':Fp_mod,'F_total':F_total_mod}


	return (Force_x, Force_y, Force_mod,Force_list,mass)

###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
def main():
	#Clear screen
	#os.system('clear')


	start_year=1995
	end_year0=1995
	Number_of_years=0
	input_folder='/ptmp/Alon.Stern/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg_bergs_Delta9/'
	second_folder=None
	#field0='CN'
	months_str='all'
	constraint_name=None
	months=range(0,11)
	file_type='ice_month'
	pole='south'
	save_Force_data=False

	load_data_from_mat=False
	#load_data_from_mat=True
	file_path='processed_data/'
	#mat_filename='processed_data/open_ocean_iceberg_forces.mat'
	mat_filename='Forces_on_berg_Tournadre_1995_1995_constraints_depth_8000_to_1000.mat'
	mat_filename='Forces_on_berg_Delta9_1995_1995_constraints_depth_8000_to_1000.mat'
	mat_filename='Forces_on_berg_Delta9_1995_1995_constraints_lon_55_to_65_lat_82_to_78_mass_1000000000000000_to_10000000000.mat'
	mat_filename='Forces_on_berg_Delta9_1995_1995_constraints_depth_8000_to_1000_mass_1000000000000000_to_10000000000.mat'
	mat_filename='Forces_on_berg_Delta9_1990_1995_constraints_lat_0_to_90_depth_8000_to_1000_mass_1000000000000000_to_10000000000.mat'
	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	Force_list=np.array(['Fa','Fo','Fi','Fr','Fc','Fp','F_total'])



 	# The constraint has the form [field_name, lower_bonud, upper_bound]
 	constraint_lon_petermann = {'Constraint_field_name': 'lon', 'lower_bound': -65 , 'upper_bound': -55, 'original_values': True} #Longitude of AB_plus_Ross
 	constraint_lat_petermann = {'Constraint_field_name': 'lat', 'lower_bound': 78 , 'upper_bound': 82, 'original_values': True} #Longitude of AB_plus_Ross
	
	constraint_depth = {'Constraint_field_name': 'depth', 'lower_bound': 1000 , 'upper_bound': 8000, 'original_values': False}
	
	constraints=[]

	constraint_mass = {'Constraint_field_name': 'mass', 'lower_bound': 10**10 , 'upper_bound': 10**15, 'original_values': False}
 	constraint_lat_baffin_coast = {'Constraint_field_name': 'lat', 'lower_bound': 30 , 'upper_bound': 80, 'original_values': False} 
 	constraint_lat_southern_hemisphere = {'Constraint_field_name': 'lat', 'lower_bound': -90 , 'upper_bound': 0, 'original_values': False} 
	#constraint_m = {'Constraint_field_name': 'mass', 'lower_bound': mass-100 , 'upper_bound': mass+100, 'original_values': True}
	#constraints.append(constraint_m)
	#constraints.append(constraint_lon_petermann)
	#constraints.append(constraint_lat_petermann)
	#constraints.append(constraint_lat_baffin_coast)
	constraints.append(constraint_lat_southern_hemisphere)
	constraints.append(constraint_depth)
	#constraints.append(constraint_mass)

	#Definging fields to be used.
	field_list=np.array(['width','length','mass','thickness','uvel','vvel','uo','vo','ui','vi','ua','va','ssh_x','ssh_y','hi'])
		

	root_path='/ptmp/aas/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg'

	(All_data, lat, lon,z) =generate_data_file(start_year, start_year, input_folder, 'CN', months_str, file_type)
	data_mean=np.squeeze(np.mean(All_data,axis=0))
	data_mean=None
	
	#Plot blank map.
	#plot_polar_field(lat,lon,data_mean,pole,difference_on=0.,title=title1)
	boundinglat=-45
	projection = 'splaea'
	m = Basemap(projection=projection, boundinglat=boundinglat, lon_0=180)

	#run_names=np.array(['tournadre']) ; naming_flag='Tournadre'
	run_names=np.array(['Delta1']) ; naming_flag='Delta1'

	if load_data_from_mat==False:
		for k in range(1):
			input_folder=root_path + '_bergs_' + run_names[k]
			#subplot(2,4,k+1)
			print input_folder	
			#Making sure the years work
			(min_year, max_year)=find_max_min_years(input_folder,'.iceberg_trajectories.nc')
			end_year=min(max_year,end_year0)
			start_year=max(end_year-Number_of_years,min_year)
			title1= run_names[k] + ' (' + str(start_year) + ' to ' + str(end_year) + ')'

			
			Total_berg_lat=np.array([])
			Total_berg_lon=np.array([])
			#Total_berg_mass0=np.array([])
			
			#initializing the Total_berg_list
			Total_berg_list={}
			for k in range(len(field_list)):
				Total_berg_list['Total_berg_'+field_list[k]]=np.array([])
				Total_berg_list[field_list[k]]=np.array([])
				#Total_berg_list['Total_berg_'+field_list[k]]=[]





			count=0
			for year in range(start_year, end_year+1):
				count=count+1
				filename ='/' + str(year) + '0101.iceberg_trajectories.nc'
				input_file=input_folder + filename
				print input_file
				
				#Handling lon and lat first
				[Total_berg_lon, Total_berg_lat]=generate_Total_berg_lon_lat(Total_berg_lon,Total_berg_lat, input_file,constraints,m)

				for k in range(len(field_list)):
					Total_berg_list[field_list[k]]=get_Total_berg_field(input_file,Total_berg_list[field_list[k]],field_list[k],constraints)


			[Force_x, Force_y, Force_mod,Force_list, mass]=calculate_forces(Total_berg_list,Total_berg_lat)
		
		if save_Force_data==True:
			mat_filename='Forces_on_berg_' + naming_flag + '_' + str(start_year)+'_'+ str(end_year)
		
			if constraints!=None:
				if len(constraints)>0:
					mat_filename=mat_filename+'_constraints'
				for k in range(len(constraints)):
					current_constraint=constraints[k]
					constraint_name=current_constraint['Constraint_field_name']
					lower_bound=current_constraint['lower_bound']
					upper_bound=current_constraint['upper_bound']
					original_values=current_constraint['original_values']
					mat_filename=mat_filename+ '_'+constraint_name + '_' + str(abs(upper_bound)) +'_to_'+ str(abs(lower_bound))
			mat_filename=mat_filename+'.mat'
			
			print 'Saving file: ' ,(file_path+mat_filename)
			sc.savemat(file_path+mat_filename, {'Force_x':Force_x , 'Force_y':Force_y , 'Force_mod':Force_mod, 'mass':mass, 'Force_list':Force_list,'Total_berg_lat':Total_berg_lat,\
				'Total_berg_lon':Total_berg_lon})

	if load_data_from_mat==True:
		print 'Loading file: ' , mat_filename
		mat_contents=sc.loadmat(file_path+mat_filename)
		Force_x=mat_contents['Force_x'][0][0]
	        Force_y=mat_contents['Force_y'][0][0]
	        Force_mod=mat_contents['Force_mod'][0][0]
	        #Force_list=mat_contents['Force_list'][:]
	        Total_berg_lon=mat_contents['Total_berg_lon'][:]
	        Total_berg_lat=mat_contents['Total_berg_lat'][:]
	        mass=mat_contents['mass'][:]


	#Converting force to acceleration
	Accel_x={}
	Accel_y={}
	Accel_mod={}
	for k in range(7):
		Accel_x[Force_list[k]]=Force_x[Force_list[k]]/mass
		Accel_y[Force_list[k]]=Force_y[Force_list[k]]/mass
		Accel_mod[Force_list[k]]=Force_mod[Force_list[k]]/mass


	#Calculating the mean and std of each force
	ind_list=np.arange(7)
	
	Force_x_mean=ind_list*0.
	Force_y_mean=ind_list*0.
	Force_mod_mean=ind_list*0.
	
	Accel_x_mean=ind_list*0.
	Accel_y_mean=ind_list*0.
	Accel_mod_mean=ind_list*0.

	Force_x_std=ind_list*0.
	Force_y_std=ind_list*0.
	Force_mod_std=ind_list*0.
	
	Accel_x_std=ind_list*0.
	Accel_y_std=ind_list*0.
	Accel_mod_std=ind_list*0.
	
	for k in range(7):
		Force_x_mean[k]=np.mean(Force_x[Force_list[k]])
		Force_y_mean[k]=np.mean(Force_y[Force_list[k]])
		Force_mod_mean[k]=np.mean(Force_mod[Force_list[k]])
		Accel_x_mean[k]=np.mean(Accel_x[Force_list[k]])
		Accel_y_mean[k]=np.mean(Accel_y[Force_list[k]])
		Accel_mod_mean[k]=np.mean(Accel_mod[Force_list[k]])
		Force_x_std[k]=np.std(Force_x[Force_list[k]])
		Force_y_std[k]=np.std(Force_y[Force_list[k]])
		Force_mod_std[k]=np.std(Force_mod[Force_list[k]])
		Accel_x_std[k]=np.std(Accel_x[Force_list[k]])
		Accel_y_std[k]=np.std(Accel_y[Force_list[k]])
		Accel_mod_std[k]=np.std(Accel_mod[Force_list[k]])

	plot_map_traj=True
	if plot_map_traj==True:
		 #plt.figure(3)
		 ax1=subplot(1,1,1)
		 plot_polar_field(lat,lon,data_mean,pole,difference_on=0.,title=''\
				 ,p_values=None,cscale=None,field=None,colorbar_on=True,return_data_map=False,plot_lat_lon_lines=True,boundinglat=boundinglat )

		 m.scatter(Total_berg_lon,Total_berg_lat,1,marker='o',color=color_vec[k])



	#F_max_lim=10**12
	F_max_lim=10**-4
	F_min_lim=-10**-4
	#F_min_lim=10**(-10)
	do_scatter_plot=False
	if do_scatter_plot==True:
		plt.figure(2)	
		for k in range(7):
		#for k in np.array([4,5]):

			#A=Force_y[Force_list[k]]
			#print Force_list[k],'max', np.max(A,axis=1) ,'min',np.min(A,axis=1)
			ax=subplot(3,3,k+1)
			ax.scatter(mass[:],Accel_y[Force_list[k]][:],color=color_vec[k],label=Force_list[k])
			#plt.plot(mass[:],Accel_y[Force_list[k]][:])
			plt.ylim([F_min_lim,F_max_lim])
			#ax.set_xscale('log')
			#ax.set_yscale('log')
			#plt.legend(loc='upper left', frameon=True)
			plt.xlabel('Mass (kg)')
			plt.ylabel('Force (N)')

	
	plot_bar_graph=False
	if plot_bar_graph==True:
		plt.figure(3)
		label_list=np.array(['Force_x','Force_y','Accel_x','Accel_y'])
		Mean_Force_and_Accel={'Force_x':Force_x_mean, 'Force_y':Force_y_mean, 'Accel_x':Accel_x_mean, 'Accel_y':Accel_y_mean}
		for k in range(4):
			ax=subplot(2,2,k+1)
			width=3
			#ax.bar(ind+((bar_width/2)*j), mean_list[:,j],width=width/2,color=color_vec[j],label=sector_title_list[j])
			#ax.bar(ind_list,Force_x_mean)
			ax.bar(ind_list,Mean_Force_and_Accel[label_list[k]],width=width/3,label=label_list[k])
			plt.title(label_list[k])
			ax.set_xticks(ind_list + 0.5)
			ax.set_xticklabels(Force_list)		

			#output_file= 'figures/traj_plot_' + str(start_year)+ '_to_' +  str(end_year) + '_with_' + run_names[k] + '.png'
			#plt.savefig(output_file, dpi=150, bbox_inches='tight', pad_inches=0.4)




	plt.show()



	print 'Script complete'

if __name__ == '__main__':
	sys.exit(main())

