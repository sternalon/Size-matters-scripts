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
from sea_ice_concentrations import generate_data_file
from Fig_observational_iceberg_plotting import plot_data_from_NIC_data_base

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

def load_data_from_mat_file(filename):
	mat_contents=sc.loadmat(filename)
	lat=mat_contents['lat'][0,:]
	lon=mat_contents['lon'][0,:]
	data=mat_contents['data'][:,:]
	data=ma.masked_values(data,1.e20)  #Recreating the mask over the data.
	field=mat_contents['field']
	#title=mat_contents['title']
	#lat_lon_bounds=mat_contents['lat_lon_bounds'][0]
	#layer_interface=mat_contents['layer_interface'][0]
	#p_values=mat_contents['p_values']
	#if p_values==0:
	#	p_values=None
	return [lat, lon, data,field ]



def create_an_empty_map():
	start_year=1947
	field='UO'
	months_str='all'
	file_type='ice_month'
	pole='south'
	input_folder='/ptmp/Alon.Stern/model_output/ulm_mom6_2015.07.20_again_myIcebergTag/AM2_LM3_SIS2_MOM6i_1deg/'
	#(All_data, lat, lon,z) =generate_data_file(start_year, start_year,input_folder,field, months_str,file_type,section_axes='xy',lat_lon_bounds=None)
	
	filename='processed_data/Control_UO_1959_to_2019_all_south_xy.mat'
	(lat ,lon, data_mean, field )=load_data_from_mat_file(filename)

	#data_mean=None
	#data_mean=np.squeeze(np.mean(All_data,axis=0))
	projection = 'splaea'
	#boundinglat=-45
	boundinglat=-55
	#Plot blank map.
	#data_mean[np.where(data_mean>0)]=1000
	#plot_polar_field(lat,lon,data_mean,pole,difference_on=1.,title='',p_values=None,cscale=None,field=None,colorbar_on=True,\
	#		return_data_map=False,plot_lat_lon_lines=False,boundinglat=boundinglat)  
	
	#m = Basemap(projection=projection, boundinglat=boundinglat, lon_0=180)
	m = Basemap(width=8000000,height=8000000,resolution='l',projection='stere', lat_ts=-70,lat_0=-90,lon_0=0)

	m.drawcoastlines()
        m.fillcontinents(color='grey',lake_color='white',zorder=0)
	#m.bluemarble()
	#m.etopo()
	m.shadedrelief()
        #m.fillcontinents(color='grey')
        #m.fillcontinents(lake_color='white')
	#m.drawmapboundary(fill_color='white')
	#m.drawlsmask(ocean_color='white',lakes=True)

	lons, lats = np.meshgrid(lon, lat)
	x, y = m(lons, lats)

	#data_mean[np.where((lons>-67)*(lons<-60))]=-1000
	data_mean[np.where((lons>-67)*(lons<-60))]=np.NaN
	data_mean[np.where((lat>-55))]=np.NaN
	data_mean[np.where((lat<-73))]=np.NaN
	#data_mean[np.where((lat<-68))]=-1000
	data_mean[np.where((lons>0)*(lons<150)*lats<-68)]=np.NaN  #Killing East Antarctica eddies
	
	levels1=np.array([0.00])

	datamap = m.contour(x, y, data_mean, levels=levels1,colors='k',linewidths=3)
	cNorm = mpl.colors.Normalize(vmin=-10., vmax=10.)
	cmap='Greys_r'
	#datamap = m.pcolor(x, y, data_mean, cmap=cmap,norm=cNorm)
	#plt.colorbar(datamap, cmap=cmap, norm=cNorm, shrink=0.5)
	return [lat,lon,m]


def plot_sectors_on_map(m,lat,sector_boundaries):
	lat_min=-75
	for k in range(len(sector_boundaries)):
		#lat=lat_min+lat		
		lat_vec=lat*0-60
		#m.plot(lon_vec,lat)
		
		if k==0:
			lat1=lat+13
		if k==1:
			lat1=lat+8
		if k==2:
			lat1=lat+7

		lon=lat1*0+sector_boundaries[k]
		x, y = m(lon, lat1)



		#m.scatter(x,y,1,marker='o',color='black')
		m.plot(x,y,linewidth=5,color='black',linestyle='--')
		#m.drawparallels(np.arange(-80.,81.,90.), zorder=1)
		#m.drawmeridians(np.arange(-180.,181.,90.), zorder=1)

		#Plotting longitude labels
		if k==0:
			x2, y2 =m(lon[0]+4,-51)
			ann1 = plt.annotate("$60^{o}$W", xy=(0,0), xycoords='data', xytext=(x2,y2), textcoords='data', size=18, va="center", ha="center")
		if k==1:
			x2, y2 =m(lon[0]+5,-51.5)
			ann1 = plt.annotate("$150^{o}$E", xy=(0,0), xycoords='data', xytext=(x2,y2), textcoords='data', size=18, va="center", ha="center")
		if k==2:  #
			x2, y2 =m(lon[0]+5,-56)
			ann1 = plt.annotate("$0^{o}$E", xy=(0,0), xycoords='data', xytext=(x2,y2), textcoords='data', size=18, va="center", ha="center")


def draw_flux_schematic(ax):
	height=10
	width=20
	xpos1=10
	xpos2=40
	xpos3=70

	rectangle1 = plt.Rectangle((xpos1, 10), width, height, fc='grey')
	rectangle2 = plt.Rectangle((xpos2, 10), width, height, fc='grey')
	rectangle3 = plt.Rectangle((xpos3, 10), width, height, fc='grey')
	plt.gca().add_patch(rectangle1)
	plt.gca().add_patch(rectangle2)
	plt.gca().add_patch(rectangle3)

	x0=xpos1+width/2
	y0=10+height
	x1=xpos1+width
	y1=10+height*2
	x2=xpos2+width/2
	y2=10+height
	x3=xpos3+width/2
	y3=10+height

	label_height=10


	ann1 = ax.annotate("", xy=(x2, y2), xycoords='data', xytext=(x1, y1), textcoords='data', size=20, va="center", ha="center", \
		arrowprops=dict(arrowstyle="simple", connectionstyle="arc3,rad=-0.4") )
	ann2 = ax.annotate("", xy=(x3, y3), xycoords='data', xytext=(x2, y2), textcoords='data', size=20, va="center", ha="center", \
		arrowprops=dict(arrowstyle="simple", connectionstyle="arc3,rad=-0.5") )
	ann3 = ax.annotate("", xy=(x3, y3), xycoords='data', xytext=(x1, y1), textcoords='data', size=20, va="center", ha="center", \
		arrowprops=dict(arrowstyle="simple", connectionstyle="arc3,rad=-0.7") )
	
	ann4 = ax.annotate("", xy=(x1, y1), xycoords='data', xytext=(x0, y0), textcoords='data', size=20, va="center", ha="center", \
		arrowprops=dict(arrowstyle="simple", connectionstyle="arc3,rad=-0.0") )

	ann4 = ax.annotate("", xy=(x2+width/4, y2-height), xycoords='data', xytext=(x2-width/4, y2-height), textcoords='data', size=20, va="center", ha="center", \
		arrowprops=dict(arrowstyle="simple", color='cyan',connectionstyle="arc3,rad=0.9") )


	ann1 = ax.annotate("F$_{IM}$", xy=(x2,y2), xycoords='data', xytext=(0.5*(x2+x3),y2+1.1* label_height), textcoords='data', size=20, va="center", ha="center")
	ann1 = ax.annotate("F$_{control}$", xy=(x3,y3), xycoords='data', xytext=(x2+width*0.3,y3+ 2.7*label_height), textcoords='data', size=20, va="center", ha="center")
	ann1 = ax.annotate("F$_{FR}$", xy=(x3,y3), xycoords='data', xytext=(x0+0.3*width*0.2,y0+1*label_height), textcoords='data', size=20, va="center", ha="center")
	ann1 = ax.annotate("F$_{C}$", xy=(x3,y3), xycoords='data', xytext=(x1+width*0.8,y0+label_height), textcoords='data', size=20, va="center", ha="center")
	ann1 = ax.annotate("T$_{I}$", xy=(x3,y3), xycoords='data', xytext=(x2,y2-1.9*label_height), textcoords='data', size=20, va="center", ha="center")
	
	ann1 = ax.annotate("Antarctica", xy=(x3,y3), xycoords='data', xytext=(x0,y3-height/2), textcoords='data', size=13, va="center", ha="center")
	ann1 = ax.annotate("Icebergs", xy=(x3,y3), xycoords='data', xytext=(x2,y3-height/2), textcoords='data', size=13, va="center", ha="center")
	ann1 = ax.annotate("Ocean", xy=(x3,y3), xycoords='data', xytext=(x3,y3-height/2), textcoords='data', size=13, va="center", ha="center")

	#set_connectionstyle("arc,angleA=0,armA=30,rad=10")
	#set_arrowstyle("Fancy,head_length=0.2")


	box_pos1=12
	box_pos2=-33
	#ax.text(-0., 55,\
	text1='F$_{FR}$     = Total Antarctic frozen runoff \nF$_{control}$ = Flux of instantaniously\n             melted runoff (control '
	text2='simulation)\nF$_{C}$      = Calving flux (iceberg simulations) \nF$_{IM}$    = Iceberg melt\nT$_{I}$      = Iceberg transport'
	long_text=text1+text2
	#ax.text(box_pos1, box_pos2,\
	#'F$_{FR}$     = Total Antarctic Frozen Runoff \nF$_{control}$ = Flux of instantaniously\n             melted runoff (control 
	#simulation)\nF$_{C}$      = Calving Flux (iceberg simulations) \nF$_{IM}$    = Iceberg Melt\nT$_{I}$      = Iceberg \
	#Transport', style='normal', bbox={'facecolor':'grey', 'alpha':0.5, 'pad':20},fontsize=13)
	
	ax.text(box_pos1, box_pos2,long_text, style='normal', bbox={'facecolor':'grey', 'alpha':0.5, 'pad':20},fontsize=13)

	#plt.gca().add_patch(circle1)
	plt.axis('scaled')
	plt.ylim([0 ,40])
	plt.axis('off')
	
	ax=gca()
	text(1,1,'(b)', ha='right', va='bottom',transform=ax.transAxes,fontsize=15)
	#plt.show()


def draw_Antarctic_maps_schematic(fig):
	(lat,lon,m)=create_an_empty_map()
	sector_boundaries=np.array([-60, -210,0])
	plot_sectors_on_map(m,lat,sector_boundaries)
	
	#Iceberg Transport Arrows	
	x, y = m(-20,-68)
	x2, y2 =m(20, -66)
	plt.annotate('', xy=(x, y),  xycoords='data',  xytext=(x2, y2), textcoords='data', size=15, arrowprops=dict(arrowstyle="simple", color='cyan',connectionstyle="arc3,rad=0.3"))
	
	x, y = m(-230,-65)
	x2, y2 =m(-195, -68)
	plt.annotate('', xy=(x, y),  xycoords='data',  xytext=(x2, y2), textcoords='data', size=15, arrowprops=dict(arrowstyle="simple", color='cyan',connectionstyle="arc3,rad=0.3"))

	#Calving Flux Arrows
	arrow_size=40

	#Weddell Sector
	x, y = m(-40,-72)
	x2, y2 =m(-34, -80)
	plt.annotate('', xy=(x, y),  xycoords='data',  xytext=(x2, y2), textcoords='data', size=arrow_size, arrowprops=dict(arrowstyle="simple",color='grey', connectionstyle="arc3,rad=0.0"))

	#East Antarctic Sector
	x, y = m(60,-63)
	x2, y2 =m(60, -71)
	plt.annotate('', xy=(x, y),  xycoords='data',  xytext=(x2, y2), textcoords='data', size=arrow_size, arrowprops=dict(arrowstyle="simple",color='g', connectionstyle="arc3,rad=0.0"))

	#AB plus Ross Sector
	x, y = m(-145,-72)
	x2, y2 =m(-145, -80)
	plt.annotate('', xy=(x, y),  xycoords='data',  xytext=(x2, y2), textcoords='data', size=arrow_size, arrowprops=dict(arrowstyle="simple",color='r', connectionstyle="arc3,rad=0.0"))


	#Labeling Sectors
	#x2, y2 =m(-40, -55)
	x2, y2 =m(-35, -52)
	ann1 = plt.annotate("Weddell Sea\nSector (WSS)", xy=(0,0), xycoords='data', xytext=(x2,y2), textcoords='data', size=18, va="center", ha="center")

	x2, y2 =m(-152, -54)
	ann1 = plt.annotate("Amundsen, Bellingshausen \nand Ross Seas Sector (ABRSS)", xy=(0,0), xycoords='data', xytext=(x2,y2), textcoords='data', size=18, va="center", ha="center")
	
	#x2, y2 =m(60, -53)
	x2, y2 =m(40, -53)
	ann1 = plt.annotate("East Antarctic \nSector (EAS)", xy=(0,0), xycoords='data', xytext=(x2,y2), textcoords='data', size=18, va="center", ha="center")


	#Labeling Arrows
	#Weddell
	x2, y2 =m(-49, -75)
	ann1 = plt.annotate("F$_{FR}$", xy=(0,0), xycoords='data', xytext=(x2,y2), textcoords='data', size=18, va="center", ha="center")
	#East Antarctica
	x2, y2 =m(67.9, -65.9)
	ann1 = plt.annotate("F$_{FR}$", xy=(0,0), xycoords='data', xytext=(x2,y2), textcoords='data', size=18, va="center", ha="center")
	#Amundsen
	x2, y2 =m(-158, -73)
	ann1 = plt.annotate("F$_{FR}$", xy=(0,0), xycoords='data', xytext=(x2,y2), textcoords='data', size=18, va="center", ha="center")
	
	#Transport Weddell
	x2, y2 =m(5, -64)
	ann1 = plt.annotate("T$_{I}$", xy=(0,0), xycoords='data', xytext=(x2,y2), textcoords='data', size=18, va="center", ha="center")
	
	#Transport Amundsen
	#x2, y2 =m(-205, -66)
	x2, y2 =m(-220, -61)
	ann1 = plt.annotate("T$_{I}$", xy=(0,0), xycoords='data', xytext=(x2,y2), textcoords='data', size=18, va="center", ha="center")

	
	#Plotting NIC database
	plot_data_from_NIC=True
	if plot_data_from_NIC==True:
		plot_data_from_NIC_data_base(lat,lon,m,plot_four_maps=False,ncols=1,fig=fig,Three_colors_on_one=True)

	ax=gca()
	text(1,1,'(a)', ha='right', va='bottom',transform=ax.transAxes,fontsize=15)

	(lat,lon,m)=create_an_empty_map()
	

##############################################################################################################
#################################### Main body of code #######################################################
##############################################################################################################

def main():

	parser = argparse.ArgumentParser()
   	args = parser.parse_args()
	
	#Flags
	save_figure=True
	plot_time_series=False
	plot_bar_figure=False
	scale_by_control=True

	#Parameters;
	
	#Delta_list=np.array(['Delta1',  'Delta9' ])
	Delta_list=np.array(['Delta1', 'Delta2', 'Delta3','Delta4', 'Delta5', 'Delta6','Delta7', 'Delta8', 'Delta9','Delta10' ])
	#Delta_list=np.array(['Delta1', 'Delta2', 'Delta3','Delta4', 'Delta5', 'Delta6','Delta7', 'Delta8', 'Delta9','Delta10','Control' ])

	Area_list=np.array(['Area$_{1}$=0.0026km$^2$', 'Area$_{2}$=0.0072km$^2$','Area$_{3}$=0.029km$^2$',\
		'Area$_{4}$=0.12km$^2$','Area$_{5}$=0.18km$^2$', 'Area$_{6}$=0.35km$^2$','Area$_{7}$=0.56km$^2$', 'Area$_{8}$=1.0km$^2$','Area$_{9}$=1.8km$^2$', 'Area$_{10}$=3.5km$^2$'])
	
	Gladstone_A=np.array([0.0026, 0.0072, 0.0292, 0.1210, 0.1788, 0.3529, 0.5647, 1.0353, 1.8353, 3.4824])  ; Gladstone_A=Gladstone_A*(10**(6))  #Convert to m^2
	color_vec=np.array(['blue', 'red', 'green', 'grey','purple', 'cyan', 'magenta', 'black', 'orange', 'coral', 'yellow', 'orchid',  'black', 'orange', 'coral', 'yellow', 'orchid' ])
	y_label_list = {'mass': 'Iceberg Mass (Gton)', 'CALVING':'Calving Flux Change (%)', 'CN':'Sea Ice Concentration Change (%)' ,'HI': 'Sea Ice Thickness Change (%)',\
			'melt':'Melt Change(%)', 'Volume':'Ice Volume Change (%)'}
	nrows=1
	ncols=1#len(Delta_list)#3
	plot_num=1
	datapath='/home/Alon.Stern/Iceberg_Project/iceberg_scripts/python_scripts/size_matters_paper/processed_data/'
	sector_list=np.array(['', '_AB_plus_Ross','_EastA_plus_AdeleLand','_Weddell'])

	#Setting up the figure
	fig, axes = plt.subplots(nrows=nrows, ncols=ncols)


	subplot(1,1,1)
	draw_Antarctic_maps_schematic(fig)

	#ax =subplot(1,2,2)

	fig.subplots_adjust(right=0.65)
	ax = fig.add_axes([0.68, 0.31 , 0.3, 0.2])
	draw_flux_schematic(ax)

	
	#Add (a) and (b) markers
	#ax.text(-85, 60.5,'(a)',fontsize=15)
	#ax.text(10, 60.5,'(b)',fontsize=15)


	fig.set_size_inches(15., 10.5,forward=True)
	if save_figure==True:
		plt.savefig('paper_figures/Fig_Antarctic_summary_map.png',dpi=300,bbox_inches='tight')
	plt.show()
	


if __name__ == '__main__':
	main()
	#sys.exit(main())

