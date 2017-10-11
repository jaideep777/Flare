# makefile for libgsm

TARGET = libgsm
VERSION = 2
LIBPATH = -L/usr/local/netcdf-cxx-legacy/lib 	# Netcdf-c++ libaray path
INCPATH = -I/usr/local/netcdf-cxx-legacy/include -I/usr/local/netcdf-c-4.3.2/include  # need paths to netcdf-c as well as c++ includes
LDFLAGS =  
CPPFLAGS = -fPIC 

LIBS = -lnetcdf_c++

SOURCEDIR = src
BUILDDIR = build
OUTLIBPATH = lib
INSTALLDIR = 

SOURCES = $(wildcard $(SOURCEDIR)/*.cpp)
OBJECTS = $(patsubst $(SOURCEDIR)/%.cpp, $(BUILDDIR)/%.o, $(SOURCES))


all: dir $(TARGET)

dir:
	mkdir -p $(BUILDDIR)

$(TARGET): $(OBJECTS)
	g++ -shared $(LIBPATH) $(LDFLAGS) -o $(OUTLIBPATH)/$(TARGET).so.$(VERSION) $(OBJECTS) $(LIBS)

$(OBJECTS): $(BUILDDIR)/%.o : $(SOURCEDIR)/%.cpp
	g++ -c $(CPPFLAGS) $(INCPATH) $< -o $@ 

install:
	cp $(OUTLIBPATH)/$(TARGET).so.$(VERSION) $(INSTALLDIR)
	
clean:
	rm -f $(BUILDDIR)/*.o $(OUTLIBPATH)/$(TARGET)*



