# makefile for libgsm

TARGET = 
LIBPATH = -L/home/jaideep/codes/Flare/lib \
		  -L/usr/local/netcdf-cxx-legacy/lib \
		  -L/usr/local/cuda/lib64	# Netcdf-c++ libaray path
		  
INCPATH = -I/usr/local/netcdf-c/include \
		  -I/usr/local/netcdf-cxx-legacy/include \
		  -I/home/jaideep/codes/Flare/include \
		  -I/usr/local/cuda/include  # need paths to netcdf-c as well as c++ includes
		  
CPPFLAGS = 

LIBS = -l:libflare.so.3 -lnetcdf_c++ -lgsl -lgslcblas 
CUDA_LIBS = -lcudart -lcurand -lcufft

SOURCES = $(TARGET).cpp

all: clean $(TARGET)

$(TARGET): 
	g++ $(INCPATH) $(LIBPATH) $(CPPFLAGS) -o $(TARGET) $(SOURCES) $(LIBS) $(CUDA_LIBS)

clean:
	rm -f $(TARGET)



