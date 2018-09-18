# FLARE Computing Library

FLARE is a GPU-accelerated library for building spatial models on 
geographical grids using NetCDF file format for data 
input/output.

This library includes functions for for performing operations
on latitude/longitude/levels and time, along with smaller
utility functions for printing data to an output stream.

This library provides the following features:

* Gridded Variables - Objects of this class include 
storage for gridded data along with meta-data, such as lat-lon
arrays and time. This library also provides operators on gVars
such as addition, multiplication etc.

* NetCDF input/output streams - These are data streams that 
can read from and write to netcdf files, 
automatically handling data distributed over multiple files and
interpolation of data into the desired grid.

* 2D Regridding (linear interpolation, coarse-graining), 
masking functions

* Time arithmatic, vector arithmatic utilities

* Initializer class - To quickly read configuration files to 
initialize model variables

* GPU-accelerated Turbulence Engine - A synthetic turbulence 
generator to create spatial heterogeneity and simulate fluid flow

* GPU-accelerated logistic resource dynamics - A spatial diffusible
resource that grows at a logistic growth rate 

* GPU-accelerated versions of grid functions. Currently, the trend 
function has a GPU version. Other functions are under 
development.



# Installation


You must have the following libraries installed to use FLARE.

* NetCDF-C++ legacy version (v4.2) 
* libgsl 
* CUDA toolkit.
=======
You must have the NetCDF-C++ legacy version (v4.2) installed before 
you proceed with installation. You will also need the CUDA toolkit.
>>>>>>> 1928751c31843ed554d25cc92f0a4a0e69b01dd8

```
mkdir builddir lib
make LIBPATH="-L/path/to/netcdf-cxx-legacy/lib -L/path/to/cuda/lib64" INCPATH="-I/path/to/netcdf-cxx-legacy/include -I/path/to/netcdf-c/include -I/path/to/cuda/include"
```
then copy the libs folder to desired location

# Acknowledgement

While using libGSM, kindly acknowledge the author Jaideep Joshi. 
The project can be found at the following link:
https://github.com/jaideep777/Flare


## Open Source Components Used

This application uses Open Source components. You can find the 
source code of their open source projects along with license information 
below. We acknowledge and are grateful to these developers for their 
contributions to open source.

Project: netCDF-c https://github.com/Unidata/netcdf-c
Copyright 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014,
University Corporation for Atmospheric Research/Unidata.


Project: netCDF-cxx4 https://github.com/Unidata/netcdf-cxx4
Copyright 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 University
Corporation for Atmospheric Research/Unidata.







