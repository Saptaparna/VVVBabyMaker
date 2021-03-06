#
# stuff to make
#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

all: 	build
help:
	@echo "Available Targets:";\
	cat Makefile | perl -ne 'printf("\t%-15s %s\n",$$1,$$2) if(/^(\S+):[^#]+(#.*)$$/)'

ifndef VERBOSE
  QUIET := @
endif

#ROOFITINCLUDE = 
#ifdef CMSSW_VERSION
#	ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
#endif

PACKAGES = coreutil

CC = g++
CMSROOT = ./
INCLUDE = $(shell root-config --cflags) -I$(CMSROOT) -I$(COREDIR) -I/cvmfs/cms.cern.ch/$(SCRAM_ARCH)/external/boost/1.57.0-ikhhed/include/ 
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(INCLUDE) $(EXTRACFLAGS) -fno-var-tracking
ROOTLIBS = $(shell root-config --ldflags --cflags --libs) -lTMVA #-lEG -lGenVector
#COREDIR = ../CORE
ifndef COREDIR
    $(error COREDIR is not set. Please export COREDIR=/abspath/to/CORE)
endif

DICTINCLUDE = $(ROOTSYS)/include/Math/QuantFuncMathCore.h $(ROOTSYS)/include/TLorentzVector.h $(ROOTSYS)/include/Math/Vector4D.h

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags --libs) -lEG -lGenVector -lTMVA -O2
# ifeq ($(shell root-config --platform),macosx)
# 	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
# endif

DIR = ./

SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:.cc=.o)
LIB = libScanChain.so

CORESOURCES = $(wildcard $(COREDIR)/*.cc)
COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = libBabymakerCORE.so

TOOLSSOURCES = $(wildcard $(COREDIR)/Tools/*.cc) $(wildcard $(COREDIR)/Tools/MT2/*.cc) $(wildcard $(COREDIR)/Tools/btagsf/*.cc) $(wildcard $(COREDIR)/Tools/datasetinfo/*.cc) $(wildcard $(COREDIR)/Tools/JetResolution/JetResolution.cxx)
TOOLSOBJECTS = $(addsuffix .o,$(basename $(TOOLSSOURCES)))
TOOLSLIB = libBabymakerTools.so

# #FWLIB = libMiniFWLite.so
# FWLIB = ../Software/MiniFWLite/libMiniFWLite.so

DICT = LinkDef_out.o

LIBS = $(LIB) $(CORELIB) $(TOOLSLIB) $(FWLIB) 

EXE = processBaby

.PHONY: all help compile clean cms2env $(PACKAGES)

#
# how to make it
#

libs:	$(LIBS)

# $(LIB):	$(OBJECTS) 
# 	echo "Linking $(LIB)"; \
# 	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@

# $(CORELIB):	$(COREOBJECTS) 
# 	echo "Linking $(CORELIB)"; \
# 	$(LINKER) $(LINKERFLAGS) -shared $(COREOBJECTS) -o $@

# $(TOOLSLIB):	$(TOOLSOBJECTS) 
# 	echo "Linking $(TOOLSLIB)"; \
# 	$(LINKER) $(LINKERFLAGS) -shared $(TOOLSOBJECTS) -o $@


$(CORELIB): $(DICT) $(COREOBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) -shared -o $@ $(COREOBJECTS) $(DICT) $(LINKERFLAGS)"; \
	$(LINKER) -shared -o $@ $(COREOBJECTS) $(DICT) $(LINKERFLAGS)

$(TOOLSLIB): $(DICT) $(TOOLSOBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) -shared -o $@ $(TOOLSOBJECTS) $(DICT) $(LINKERFLAGS)"; \
	$(LINKER) -shared -o $@ $(TOOLSOBJECTS) $(DICT) $(LINKERFLAGS)

$(LIB):	$(DICT) $(OBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) -shared -o $@ $(OBJECTS) $(DICT) $(LINKERFLAGS)"; \
	$(LINKER) -shared -o $@ $(OBJECTS) $(DICT) $(LINKERFLAGS)

$(PACKAGES):
	make -j 15 -C $@
	cp $@/lib$@.so .

# $(FWLIB):
# 	echo "making MiniFWLite"; \
#         cd ../Software/MiniFWLite; \
#         $(MAKE) -f Makefile; cd -; \
# 	cp $(FWLIB) .;

LinkDef_out.cxx: LinkDef.h
	$(QUIET) echo "Making CINT dictionaries"; \
	rootcint -f LinkDef_out.cc -c -p $(DICTINCLUDE)  LinkDef.h; \
	cat LinkDef_out.cc > LinkDef_out.cxx; rm LinkDef_out.cc

# LinkDef_out.cxx: LinkDef.h ScanChain.h
# 	rootcint -f $@ -c $(INCLUDE) ScanChain.h $<

# General rule for making object files
# %.d:	%.cc
# 	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
#                      [ -s $@ ] || rm -f $@
# %.d:	%.cxx
# 	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
#                      [ -s $@ ] || rm -f $@

# %.o: 	%.cc 
# 	$(CC) $(CFLAGS) $< -c -o $@

# %.o: 	%.cxx
# 	$(CC) $(CFLAGS) $< -c -o $@

# .PHONY: all
# all:	$(LIBS)  

# .PHONY: clean

# the option "-Wl,-rpath,./" adds ./ to the runtime search path for libraries
$(EXE): $(PACKAGES) $(LIBS)
	$(QUIET) echo "Building $@"; \
	echo "$(CC) -o $@ $(LIBS) $(ROOTLIBS) -Wl,-rpath,./ $(addprefix -L,$(PACKAGES)) $(addprefix -l,$(PACKAGES))"; \
	$(CC) -o $@ $(LIBS) $(ROOTLIBS) -Wl,-rpath,./  $(addprefix -L,$(PACKAGES)) $(addprefix -l,$(PACKAGES))

%.o: 	%.cc %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.C %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cc
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.C
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o:    %.cxx 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $(ROOTLIBS) $< -c -o $@

libs:	$(LIBS)

build:  $(EXE)

b: build

clean:
	rm -f \
	LinkDef_out* \
	*.o \
	$(LIBS) \
	$(EXE) \
	roccor.2017.v0/*.o \
	coreutil/*.o \
	coreutil/*.so \
	$(COREDIR)/*.o \
	$(COREDIR)/*.d \
	$(COREDIR)/*.so \
	$(COREDIR)/Tools/*.o \
	$(COREDIR)/Tools/*.d \
	$(COREDIR)/Tools/*.so \
	$(COREDIR)/Tools/MT2/*.o \
	$(COREDIR)/Tools/MT2/*.d \
	$(COREDIR)/Tools/MT2/*.so \
	$(COREDIR)/Tools/JetResolution/*.o \
	$(COREDIR)/Tools/JetResolution/*.d \
	$(COREDIR)/Tools/JetResolution/*.so

localclean:
	rm -f \
	LinkDef_out* \
	*.o \
	*.so 

# test: 
# 	@echo $(INCLUDE)

# -include $(LIBDIR)/LinkDef_out.d

endif
