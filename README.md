![Alt text](simepar.png)


# WRITE RSL/TRMM RADAR TO HDF5

(C) Copyleft Sistema Meteorológico Simepar (SIMEPAR) See LICENSE.txt for more details

In this folder you will find a code in C to write an [RSL/TRMM](http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/) radar extruture into an [HDF5](http://www.hdfgroup.org/HDF5) file following the GAMIC convention

## DEPENDENCIES

In order to compile this programs you will need already installed:

* an C compiler, preferentially GCC 
* RSL/TRMM version 1.41 or later:  [http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/](http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/)
* HDF5 C library version 1.8 or later: [http://www.hdfgroup.org/HDF5/release/obtain5.html#obtain](http://www.hdfgroup.org/HDF5/release/obtain5.html#obtain)

## BUILDING


Download the files in this folder to your machine.
 
You may use this code direct in our aplication or create a dynamic library as follow:

First of all some times RSL doesn't install itself in standard locations, if this is your case update the first line of MAKEFILE to the instalation place of RSL

In a terminal at the downloaded folder execute **make**. This will create dynamic library *librsl_to_hdf5gamic.so*, move this file to the standard library location of your system as well as the file *rsl_to_hdf5gamic.h*, or link direct to it when compiling your aplication.

## USING

When writing your aplication include the header

    #include <rsl_to_hdf5gamic.h>

and call 

    radar_to_hdf5gamic(Radar* radar, char *outfile)
    
to write an HDF5 file

When compiling your aplication:

- If you have the dynamic library, link to it with the flag *-lrsl_to_hdf5gamic*. You must also link RSL and HDF5 with *-lrsl -lhdf5*

- If want to use the pure code add the file *rsl_to_hdf5gamic.c* in your compilation line and link to RSL and HDF5  with  *-lrsl -lhdf5*

