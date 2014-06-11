/*
	|||||   ||||||||||||||||| SIMEPAR
	|||||   ||||||||||||||||| Sistema Meteorológico do Paraná
	|||||               |||||
	|||||||||||||||||   ||||| Technology and Environmental Information


	(C) Copyleft Sistema Meteorológico Simepar (SIMEPAR)
        http://www.simepar.br
 
	Author: Anderson Luis Gama
	Manager: Cesar Beneti

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 3 of the License (29 June 
        2007), or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program (LICENCE.txt); if not, write to the Free 
	Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#define USE_RSL_VARS
#include <rsl.h>
#include <hdf5.h>

//#include "../include/radarreader.h"
typedef struct {
	float azimuth_start;
	float azimuth_stop;
	float elevation_start;
	float elevation_stop;
	unsigned long long timestamp;
	}RayHeader;



/**************** PROTOTYPE ****************/
void radar_to_hdf5gamic(Radar* radar, char *outfile);// This is the only function intended to be call externally, it receives a RSL/TRMM radar structure (trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/) and write a HDF5 file in the Gamic Convention named outfile. In respect to the convention this function expect all volumes of the radar strucure to have the same number of sweeps. Furthermore although compatible with RSL V1.44, only volumes with index smaller than 20 are converted, this is dow to the fact that we are not sure how to name the other volumes in the HDF5 file. To extend this functionaly update the variable rsl_to_name_hdf5

//This functions automate the tedious process of writing to a HDF5 file
void write_attr_text(hid_t loc_id,char *name,char * value);//write atribute text to HDF5
void write_attr_float(hid_t loc_id,char *name,float value);//write atribute float to HDF5
void write_attr_uint(hid_t loc_id,char *name,int value);//write atribute unsignet int to HDF5
void write_attr_double(hid_t loc_id,char *name,double value);//Write atribute double to HDF5


RayHeader *make_ray_header(Sweep * sweep);//alloc and fill with vector of rayheader in a buffer that can be writen to HDF5, make sure to free it after
unsigned char *make_buffer(int moment,Sweep * sweep,double max, double min);//alloc and fill rays*bins buffer with moments data in scaled char to write to HDF5 /scanXX/moment_YY dataset,  make sure to free it after

int get_unfolding(Ray * ray);//Use prf and prf2 from ray header to calculate unfolding code according with HDF5 convertion
unsigned long long int unix_time(Ray_header header);//convert time in a ray to unix_time (in microseconds)
double get_max(Sweep *sweep);// get the biggers value in sweep;
double get_min(Sweep *sweep);// get the smallest value in sweep;


/**************** CODE *********************/

void radar_to_hdf5gamic(Radar* radar, char *outfile){
	extern int radar_verbose_flag; /* Defined by RSL */  

	char string[100]; //all uses string
	char moment[20],scan[20];//moment and scan name string
	hid_t file_id,what_id,how_id,where_id,scan_id;//specific id's, although scan_id may change (roll) during the program 
	hid_t group_id, attr_id, dataspace_id,dataset_id; //all uses id
	hsize_t dims[2];//dimensions
	int i,j,k,len;// len: length of strings
	int nsweeps, nvolume=0;
	int moments[200];//moment that will be writen, in rsl numeration
	#define MAX_NVOLUME 20
	char *rsl_to_name_hdf5[MAX_NVOLUME]={"Z","V","W","AZh","UZ","ZDR","LDR","ZDR","SIGPOW","RHOHV","PHIDP","XZ","UZDR","MZ","MR","AZh","UVh","KDP","TI","SQI"};
	herr_t status;
	RayHeader * rays_header;
	unsigned char * buffer;
	double min,max;
	j=0;
	for(i=0;i<MAX_NVOLUME;i++)if(radar->v[i]!=NULL){moments[j]=i;j++;}//get list of not NULL volumes in Radar
	nvolume=j;//number of volumes to be writen
	nsweeps=radar->v[moments[0]]->h.nsweeps;//number of scans

		
	//Overwrite file if already exist
	file_id = H5Fcreate(outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if(file_id<0){ printf("Falling Creating a HDF5 file\n"); exit(-1);}
	if(radar_verbose_flag) printf("Start writing HDF5\n");
	
	//creat groups
	what_id=H5Gcreate2(file_id,"/what", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//create, white and close atributes of /what
		write_attr_text(what_id,"object","PVOL");
		write_attr_text(what_id,"version","1");
		time_t now=time(NULL);
		struct tm *now_tm=gmtime (&now);
		strftime (string, 100, "%FT%TZ",now_tm);
		write_attr_text(what_id,"date",string);
		write_attr_uint(what_id,"sets",nsweeps);
		H5Gclose(what_id);
		
	how_id=H5Gcreate2(file_id,"/how", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//create, white and close atributes of /how
		write_attr_double(how_id,"elevation_beam",radar->v[moments[0]]->sweep[0]->h.vert_half_bw);
		write_attr_double(how_id,"azimuth_beam",radar->v[moments[0]]->sweep[0]->h.horz_half_bw);
		write_attr_text(how_id,"site_name",radar->h.radar_name);
		H5Gclose(how_id);
		
	where_id=H5Gcreate2(file_id,"/where", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//create, white and close atributes of /where
		write_attr_double(where_id,"lat",radar->h.latd+(radar->h.latm/60)+(radar->h.lats/3600));
		write_attr_double(where_id,"lon",radar->h.lond+(radar->h.lonm/60)+(radar->h.lons/3600));
		write_attr_double(where_id,"height",radar->h.height);
		H5Gclose(where_id);
	
	for(i=0;i<nsweeps;i++){
		if(radar_verbose_flag) printf("Writing /scan%i...\n",i);
		sprintf(scan,"/scan%d",i); 
		scan_id=H5Gcreate2(file_id,scan, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		group_id=H5Gcreate2(scan_id,"what", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		//create, white and close atributes of /scanXX/what
			write_attr_text(group_id,"product","SCAN");
			write_attr_text(group_id,"scan_type","PPI");
			write_attr_uint(group_id,"descriptor_count",nvolume);
			H5Gclose(group_id);
			
		group_id=H5Gcreate2(scan_id,"how", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		//create, white and close atributes of /scanXX/how
			sprintf(string, "%04i-%02i-%02iT%02i:%02i:%02iZ", radar->v[moments[0]]->sweep[i]->ray[0]->h.year,
							radar->v[moments[0]]->sweep[i]->ray[0]->h.month,radar->v[moments[0]]->sweep[i]->ray[0]->h.day,
							radar->v[moments[0]]->sweep[i]->ray[0]->h.hour,radar->v[moments[0]]->sweep[i]->ray[0]->h.minute,
							(int)radar->v[moments[0]]->sweep[i]->ray[0]->h.sec);
			write_attr_text(group_id,"timestamp",string);
			write_attr_uint(group_id,"bin_count",radar->v[moments[0]]->sweep[i]->ray[0]->h.nbins);
			write_attr_double(group_id,"range_start",radar->v[moments[0]]->sweep[i]->ray[0]->h.range_bin1);
			write_attr_double(group_id,"scan_speed",radar->v[moments[0]]->sweep[i]->ray[0]->h.azim_rate);
			write_attr_double(group_id,"range",radar->v[moments[0]]->sweep[i]->ray[0]->h.range_bin1+radar->v[moments[0]]->sweep[i]->ray[0]->h.gate_size*radar->v[moments[0]]->sweep[i]->ray[0]->h.nbins);
			
			write_attr_double(group_id,"range_step",radar->v[moments[0]]->sweep[i]->ray[0]->h.gate_size);
			write_attr_uint(group_id,"range_samples",-1);
			write_attr_uint(group_id,"PRF",radar->v[moments[0]]->sweep[i]->ray[0]->h.prf);
			write_attr_uint(group_id,"angle_sync",1);
			write_attr_uint(group_id,"time_samples",-1);
			if(radar->v[moments[0]]->sweep[i]->h.nrays>1);
			write_attr_double(group_id,"angle_step",radar->v[moments[0]]->sweep[i]->ray[1]->h.azimuth-radar->v[moments[0]]->sweep[i]->ray[0]->h.azimuth);
			write_attr_uint(group_id,"unfolding ",get_unfolding(radar->v[moments[0]]->sweep[i]->ray[0]));
			write_attr_uint(group_id,"filter",-1);
			write_attr_double(group_id,"radar_wave_length ",radar->v[moments[0]]->sweep[i]->ray[0]->h.wavelength);
			write_attr_uint(group_id,"ray_count",radar->v[moments[0]]->sweep[i]->h.nrays);
			write_attr_double(group_id,"azi_start",radar->v[moments[0]]->sweep[i]->ray[0]->h.azimuth);
			write_attr_double(group_id,"azi_stop",radar->v[moments[0]]->sweep[i]->ray[radar->v[moments[0]]->sweep[i]->h.nrays-1]->h.azimuth);
			write_attr_double(group_id,"elevation",radar->v[moments[0]]->sweep[i]->ray[0]->h.elev);
			H5Gclose(group_id);
			
			//Now Ray_Header
			//create compond type					
			hid_t  ray_header_mem,ray_header_file;
				//create memory type
				ray_header_mem = H5Tcreate(H5T_COMPOUND, sizeof(RayHeader));
   				status = H5Tinsert(ray_header_mem, "azimuth_start", HOFFSET(RayHeader, azimuth_start), H5T_NATIVE_FLOAT);
				status = H5Tinsert(ray_header_mem, "azimuth_stop", HOFFSET(RayHeader, azimuth_stop), H5T_NATIVE_FLOAT);
				status = H5Tinsert(ray_header_mem, "elevation_start", HOFFSET(RayHeader, elevation_start), H5T_NATIVE_FLOAT);
				status = H5Tinsert(ray_header_mem, "elevation_stop", HOFFSET(RayHeader, elevation_stop), H5T_NATIVE_FLOAT);		 
				status = H5Tinsert(ray_header_mem, "timestamp", HOFFSET(RayHeader, timestamp), H5T_NATIVE_ULLONG); 
				//create file type
				ray_header_file = H5Tcreate (H5T_COMPOUND, 4+4+4+4+8);
				status = H5Tinsert (ray_header_file, "azimuth_start", 0, H5T_IEEE_F32LE);
				status = H5Tinsert (ray_header_file, "azimuth_stop", 4, H5T_IEEE_F32LE);
				status = H5Tinsert (ray_header_file, "elevation_start", 4+4, H5T_IEEE_F32LE);
				status = H5Tinsert (ray_header_file, "elevation_stop", 4+4+4, H5T_IEEE_F32LE);
				status = H5Tinsert (ray_header_file, "timestamp", 4+4+4+4,  H5T_STD_U64LE);
				//create Space, memory, write and close
					rays_header=make_ray_header(radar->v[moments[0]]->sweep[i]);
					dataspace_id = H5Screate_simple(1,(const hsize_t*) &radar->v[moments[0]]->sweep[i]->h.nrays, NULL);
					attr_id = H5Dcreate (scan_id, "ray_header", ray_header_file, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
					status = H5Dwrite (attr_id, ray_header_mem, H5S_ALL, H5S_ALL, H5P_DEFAULT, rays_header);
				   	status = H5Dclose(attr_id);
					status = H5Sclose(dataspace_id); 
					free(rays_header);
			
		for(j=0;j<nvolume;j++){
			if(radar_verbose_flag)printf("        /scan%i/moment_%i:  %s\n",i,j,rsl_to_name_hdf5[moments[j]]);
			sprintf(moment, "moment_%d", j);
			//printf("moment %i",moments[j]);
			//create dataset /scanXX/moment_YY
			dims[0]=radar->v[moments[j]]->sweep[i]->h.nrays;
			dims[1]=radar->v[moments[j]]->sweep[i]->ray[0]->h.nbins;
			dataspace_id = H5Screate_simple(2, dims, NULL);
			dataset_id = H5Dcreate2(scan_id, moment, H5T_STD_U8LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			//write attribute to /scanXX/moment_YY
			write_attr_text(dataset_id,"moment",rsl_to_name_hdf5[moments[j]]);
			write_attr_text(dataset_id,"format","UV8");
			max=get_max(radar->v[moments[j]]->sweep[i]);
			write_attr_double(dataset_id,"dyn_range_max",max);
			min=get_min(radar->v[moments[j]]->sweep[i]);
			write_attr_double(dataset_id,"dyn_range_min",min);
			//get date from sweep, write and close
			buffer=make_buffer(moments[j],radar->v[moments[j]]->sweep[i],max,min);
			status = H5Dwrite(dataset_id, H5T_STD_U8LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
			status = H5Dclose(dataset_id);
			status = H5Sclose(dataspace_id); 
			free(buffer);
			
		}
		status = H5Gclose(scan_id);
	}
	if(radar_verbose_flag)printf("Closing file...");
	status = H5Fclose(file_id);
	if(radar_verbose_flag)printf("OK\n");
}	
	
	
/********* ATTRIBUTE WRITING FUNCTIONS **********/
void write_attr_text(hid_t loc_id,char *name,char * value){//write atribute text
	herr_t status; 
	const hsize_t len=strlen(value); 
	hid_t dataspace_id = H5Screate_simple(1, &len, NULL); 
	hid_t attr_id=H5Acreate2(loc_id,name,H5T_C_S1,dataspace_id,H5P_DEFAULT, H5P_DEFAULT); 
	status = H5Awrite(attr_id,H5T_C_S1, value);
   	status = H5Aclose(attr_id);
	status = H5Sclose(dataspace_id); 
}
void write_attr_float(hid_t loc_id,char *name,float value){//write atribute float
	herr_t status;
	const hsize_t len=1;
	hid_t dataspace_id = H5Screate_simple(1, &len, NULL);
	hid_t attr_id=H5Acreate2(loc_id,name,H5T_IEEE_F32LE,dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr_id,H5T_NATIVE_FLOAT, &value);
   	status = H5Aclose(attr_id);
	status = H5Sclose(dataspace_id); 
}
void write_attr_uint(hid_t loc_id,char *name,int value){//write atribute unsignet int
	herr_t status; 
	const hsize_t len=1;
	hid_t dataspace_id = H5Screate_simple(1, &len, NULL);
	hid_t attr_id=H5Acreate2(loc_id,name,H5T_STD_U32LE,dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr_id,H5T_NATIVE_UINT, &value);
   	status = H5Aclose(attr_id);
	status = H5Sclose(dataspace_id); 
}
void write_attr_double(hid_t loc_id,char *name,double value){//write atribute double
	herr_t status; 
	const hsize_t len=1;
	hid_t dataspace_id = H5Screate_simple(1, &len, NULL);
	hid_t attr_id=H5Acreate2(loc_id,name,H5T_IEEE_F64LE,dataspace_id,H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr_id,H5T_NATIVE_DOUBLE, &value);
   	status = H5Aclose(attr_id);
	status = H5Sclose(dataspace_id); 
}

/********* BUFFER MAKING FUNCTIONS **********/

RayHeader *make_ray_header(Sweep * sweep){
	int nrays=sweep->h.nrays;
	int i,j,k;
	RayHeader * rays_header=(RayHeader *)malloc(sizeof(RayHeader)*nrays);
	if (rays_header==NULL){printf("Memory Allocation Error. I'm leaving you"); exit(-1);}
	for(i=0;i<nrays;i++){//
		rays_header[i].azimuth_start=sweep->ray[i]->h.azimuth;
		rays_header[i].azimuth_stop=sweep->ray[i]->h.azimuth;
		rays_header[i].elevation_start=sweep->ray[i]->h.elev;
		rays_header[i].elevation_stop=sweep->ray[i]->h.elev;
		rays_header[i].timestamp=unix_time(sweep->ray[i]->h);
	}
	return rays_header;
}

unsigned char *make_buffer(int moment,Sweep * sweep,double max, double min){
	int nrays,nbins,i,j;
	nrays=sweep->h.nrays;
	nbins=sweep->ray[0]->h.nbins;
	float value;
	unsigned char scaled;
//	printf("%i  %f\n",moment,sweep->h.elev);
	unsigned char * buffer=(unsigned char *)malloc(sizeof(char)*nrays*nbins);
	if (buffer==NULL){printf("Memory Allocation Error. I'm leaving you"); exit(-1);}
	for(i=0;i<nrays;i++){
		for(j=0;j<nbins;j++){
			value=sweep->ray[i]->h.f(sweep->ray[i]->range[j]);
			if(value<131000){ //actually 131072 but, as one says, when it floats go safe
				scaled=(value-min)*(254)/(max-min)+1;
				if(scaled>255) scaled=255;
			}
			else
				scaled=0;
			buffer[i*nbins+j]=scaled;
		}
	}
	//printf("max %f, min %f\n",max,min);
	return buffer;
			
}
/********* ACCESSORY FUNCTIONS **********/

int get_unfolding(Ray * ray){
	if (ray->h.prf2<1) return 0;
	float frac=ray->h.prf/ray->h.prf2;
	if(frac>17/12) return 1;
	if(frac>31/24) return 2;
	if(frac>9/8) return 3;
	return 0;
}

unsigned long long int unix_time(Ray_header header){
	time_t time;
	struct tm tmTime;
	double time2;
	unsigned long long int time3;
	
	tmTime.tm_sec=0;	
	tmTime.tm_min=header.minute	;
	tmTime.tm_hour=header.hour	;
	tmTime.tm_mday=header.day	;
	tmTime.tm_mon=header.month-1;
	tmTime.tm_year=header.year-1900;	
	tmTime.tm_wday=-1	;
	tmTime.tm_yday=-1;
	tmTime.tm_isdst=-1;
	
	time=mktime(&tmTime);
	time=time-timezone;
	time2=(((double)time)+header.sec);
	time3=time2*pow(10,6);
	return time3;

}



double get_max(Sweep *sweep){

	int i,j;
	float max,current=sweep->ray[0]->h.f(sweep->ray[0]->range[0]);
	double ret;
	max=-131000;
	for(i=0;i<sweep->h.nrays;i++){
		for(j=0;j<sweep->ray[i]->h.nbins;j++){
			current=sweep->ray[i]->h.f(sweep->ray[i]->range[j]);
			if(max<current && current<131000) max=current;
		}
	}
	ret=max;
	return ret;
}

double get_min(Sweep *sweep){

	int i,j;
	float min,current=sweep->ray[0]->h.f(sweep->ray[0]->range[0]);
	double ret;
	min=131000;
	for(i=0;i<sweep->h.nrays;i++){
		for(j=0;j<sweep->ray[i]->h.nbins;j++){
			current=sweep->ray[i]->h.f(sweep->ray[i]->range[j]);
			//if(current<131000) printf("current %f  sweep %f ray %i range %i \n",current,sweep->h.elev,i,j);
			if(min>current) min=current;
		}
	}
	ret=min;
	return min;
}

