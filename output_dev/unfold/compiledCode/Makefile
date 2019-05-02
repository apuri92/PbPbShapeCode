CXX      = g++
CXXFLAGS = -Wall `root-config --cflags` -g -o2
LDFLAGS  = -Wall `root-config --glibs`

LIBDIRS     = -L/usr/local/lib/


OBJ         = conf_plots
SOURCES     = ManyHist.cpp drawingClass.cpp main.cpp
INCLUDES    = 

DSOURCES    = ./
DINCLUDES   = 

CSOURCES    = $(addprefix $(DSOURCES),$(SOURCES))
CINCLUDES   = $(addprefix $(DINCLUDES),$(INCLUDES))

all: $(CSOURCES) $(CINCLUDES) 
	$(CXX) $(CXXFLAGS) $(CSOURCES) $(INCLUDEDIRS) \
	$(LIBDIRS) $(LDFLAGS) -o $(OBJ)

clean : 
	rm $(OBJ)
