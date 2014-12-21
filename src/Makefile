# makefile for GSM library

OUTLIBPATH = /home/jaideep/libgsm/lib
OUTINCPATH = /home/jaideep/libgsm/include

#OUTLIBPATH = /usr/local/lib
#OUTINCPATH = /usr/local/include

TARGET = libgsm
LIBPATH = -L/usr/local/netcdf-4.1.3/lib
INCPATH = -I/usr/local/netcdf-4.1.3/include
LDFLAGS =  
CPPFLAGS = -fPIC 

CPPFILES = globals.cpp arrayutils.cpp time.cpp gvar.cpp grid.cpp ncio.cpp 
OBJFILES = $(CPPFILES:.cpp=.o)

LIBS = -lnetcdf_c++


all: install

install: $(OBJFILES)
	g++ -shared -o $(LDFLAGS) $(TARGET).so $(OBJFILES) $(LIBPATH) $(LIBS)
	cp $(TARGET).so $(OUTLIBPATH)

.cpp.o:
	g++ -c $(CPPFLAGS) $(INCPATH) $< 

clean:
	rm -f *.o *.so
	rm -f $(OUTLIBPATH)/$(TARGET).so
#	rm -f $(INCPATH)/gsm.h

re: clean install


