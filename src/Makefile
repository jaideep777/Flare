# makefile for GSM library

OUTLIBPATH = ../lib
OUTINCPATH = ../include

#OUTLIBPATH = /usr/local/lib
#OUTINCPATH = /usr/local/include

TARGET = libgsm
LIBPATH = -L/usr/local/netcdf-cxx-legacy/lib 	# Netcdf-c++ libaray path
INCPATH = -I/usr/local/netcdf-cxx-legacy/include -I/usr/local/netcdf-c-4.3.2/include  # need paths to netcdf-c as well as c++ includes
LDFLAGS =  
CPPFLAGS = -fPIC 

CPPFILES = globals.cpp arrayutils.cpp time.cpp gvar.cpp grid.cpp ncio.cpp 
OBJFILES = $(CPPFILES:.cpp=.o)

LIBS = -lnetcdf_c++


all: install

install: $(OBJFILES)
	g++ -shared $(LIBPATH) $(LDFLAGS) -o $(TARGET).so $(OBJFILES) $(LIBS)
	cp $(TARGET).so $(OUTLIBPATH)

.cpp.o:
	g++ -c $(CPPFLAGS) $(INCPATH) $< 

clean:
	rm -f *.o *.so
	rm -f $(OUTLIBPATH)/$(TARGET).so
#	rm -f $(INCPATH)/gsm.h

re: clean install


