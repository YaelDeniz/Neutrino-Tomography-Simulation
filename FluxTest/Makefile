#MACROS:  Expression which is likely to be used repeatedl
CC= g++

# Source and build output
SRC_DIR = .
MAIN_DIR= .
OBJ_DIR = .
BIN_DIR = .

#INCLUDE PATHS
ROOT_LIB=/home/dehy0499/root/lib
ROOTMATH_INCLUDE=/home/dehy0499/root/include/Math
#OSCPROB_PATH = /home/dehy0499/Desktop/Tomography/Simulation/OscProb
#NUTOMO_INCLUDE = /home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/inc


#Flags
ROOT_FLAGS = `root-config  --cflags --glibs --libs`
CPPFLAGS = -I$(ROOT_LIB) -I$(ROOTMATH_INCLUDE) -I$(OSCPROB_PATH) -I$(NUTOMO_INCLUDE) # Where to find /inc
CFLAGS   := -Wall
LDFLAGS  := -L$(OSCPROB_PATH) -L$(ROOT_LIB) # Where to find /lib
LDLIBS   := -lMathMore  -lGeom# load specific lib


#OBJECTS = $(SOURCES:.c=.o)

OBJECTS_Test = NuFlux.o 



#OBJECTSProb = GetProbPREM.o


.PHONY: all clean

#TARGETS

all: Test #DetectorResolution   # All targets


Test: $(MAIN_DIR)/FluxMain.cpp $(OBJECTS_Test) 
	@echo " "
	@echo "Generating Application for FluxTest"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "Done"


#Build Objects

NuFlux.o: $(SRC_DIR)/NuFlux.C  
	@echo "  "
	@echo "Creating Objects"
	$(CC) -g -c $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@




clean:
	@echo "Removing Objects"
	rm *.o  Test
