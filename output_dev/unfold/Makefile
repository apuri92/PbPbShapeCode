ROOUNFOLDDIR = /Applications/root_v6.10.08/RooUnfold

CXX      = g++ 
CXXFLAGS = -Wall `root-config --cflags` -g -o2
LDFLAGS  = -Wall `root-config --glibs`
INCLUDEDIRS = -I$(ROOUNFOLDDIR)/src

LIBDIRS     = -L$(ROOUNFOLDDIR) -L/usr/local/lib/
LIBS        = -lRooUnfold

OBJ         = analysis
SOURCES     = main.cpp
INCLUDES    = 

DSOURCES    = src/
DINCLUDES   = 

CSOURCES    = $(addprefix $(DSOURCES),$(SOURCES))
CINCLUDES   = $(addprefix $(DINCLUDES),$(INCLUDES))

all: $(CSOURCES) $(CINCLUDES) 
	$(CXX) $(CXXFLAGS) $(CSOURCES) $(INCLUDEDIRS) \
	$(LIBDIRS) $(LIBS) $(LDFLAGS) -o $(OBJ)

clean : 
	rm $(OBJ)
