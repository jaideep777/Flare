Gridded Spatial Modelling library
=================================

libGSM is a simple library for building spatial models on 
geographical grids using NetCDF file format for data 
input/output.

This library includes functions for for performing operations
on latitude/longitude/levels and time, along with smaller
utility functions for printing data to an output stream.

This library provides two classes:

1) Gridded Variable (gVar) - Objects of this class include 
storage for gridded data along with meta-data, such as lat-lon
arrays and time. This library also provides operators on gVars
such as addition, multiplication etc.

2) A handle for netCDF files (NcFile_handle) - This class
provides wrappers for easily reading netcdf files into gridded
variables.

Spatial models with NetCDF can be quickly and easily written 
using libGSM. 



Acknowledgement
===============

While using libGSM, kindly acknowledge the author Jaideep Joshi. 
The project can be found at the following link:
https://github.com/jaideep777/Gridded-Spatial-Modelling-library 


Open Source Components Used
===========================

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





