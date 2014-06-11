![Alt text](simepar.png)


# WRITE RSL/TRMM RADAR TO HDF5

(C) Copyleft SIMEPAR - Sistema Meteorologico do Parana. See LICENSE.txt for more details

In this folder you will find a code in C to write a [RSL/TRMM](http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/) radar structure into an [HDF5](http://www.hdfgroup.org/HDF5) file according to the GAMIC convention

## DEPENDENCIES

In order to compile this programs you will need, already installed:

* a C compiler, preferentially GCC 
* RSL/TRMM version 1.41 or later:  [http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/](http://trmm-fc.gsfc.nasa.gov/trmm_gv/software/rsl/)
* HDF5 C library version 1.8 or later: [http://www.hdfgroup.org/HDF5/release/obtain5.html#obtain](http://www.hdfgroup.org/HDF5/release/obtain5.html#obtain)

## BUILDING


Download the files in this folder to your machine.
 
You may use this code directly in our application or create a dynamic library as follows:

First of all, sometimes RSL does not install itself in standard locations. if this is your case, update the first line of MAKEFILE to the installation location of RSL

In a terminal at the downloaded folder execute **make**. This will create dynamic library *librsl_to_hdf5gamic.so*. Move this file to the standard library location of your system as well as the file *rsl_to_hdf5gamic.h*, or link directly to it when compiling your application.

## USING

When writing your own application include the header

    #include <rsl_to_hdf5gamic.h>

and call 

    radar_to_hdf5gamic(Radar* radar, char *outfile)
    
to write an HDF5 file

When compiling your application:

- If you have the dynamic library, link with the flag *-lrsl_to_hdf5gamic*. You must also link RSL and HDF5 with *-lrsl -lhdf5*

- If want to use the pure code add the file *rsl_to_hdf5gamic.c* in your compilation line and link to RSL and HDF5  with  *-lrsl -lhdf5*