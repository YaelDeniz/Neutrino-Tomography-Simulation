#MACROS:  Expression which is likely to be used repeatedl
CC= g++

# Source and build output
SRC_DIR = ./src
MAIN_DIR= ./main
OBJ_DIR = ./obj
BIN_DIR = .

#INCLUDE PATHS
ROOT_LIB=/home/dehy0499/root/lib
ROOTMATH_INCLUDE=/home/dehy0499/root/include/Math
OSCPROB_PATH = /home/dehy0499/Desktop/Tomography/Simulation/OscProb
NUTOMO_INCLUDE = /home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/inc


#Flags
ROOT_FLAGS = `root-config  --cflags --glibs --libs`
CPPFLAGS = -I$(ROOT_LIB) -I$(ROOTMATH_INCLUDE) -I$(OSCPROB_PATH) -I$(NUTOMO_INCLUDE) # Where to find /inc
CFLAGS   := -Wall
LDFLAGS  := -L$(OSCPROB_PATH) -L$(ROOT_LIB) # Where to find /lib
LDLIBS   := -lOscProb -lMathMore  -lGeom# load specific lib


#OBJECTS = $(SOURCES:.c=.o)

OBJECTS_AsivTrue = MathTools.o PhyTools.o Earth3DModel.o AsimovDataTrue.o 
OBJECTS_SensitivityTrue = MathTools.o PhyTools.o Earth3DModel.o AsimovDataTrue.o 
OBJECTS_AsivObs  = MathTools.o PhyTools.o Earth3DModel.o AsimovDataTrue.o AsimovDataObs.o 
OBJECTS_SensitivityObs  = MathTools.o PhyTools.o Earth3DModel.o AsimovDataTrue.o AsimovDataObs.o 

OBJECTS_Res  = MathTools.o PhyTools.o  DetectorResolution.o 
OBJECTS_OscProb  = MathTools.o PhyTools.o  OscProbEarth.o 
OBJECTS_Flux = PhyTools.o
OBJECTS_EARTHMODEL = Earth3DModel.o OscProbEarth.o 


#OBJECTSProb = GetProbPREM.o


.PHONY: all clean

#TARGETS

all: ObservedEvents TrueEvents  OscProbEarth HondaFlux Earth3D EarthSensitivity SensitivityTrue SensitivityObs #DetectorResolution   # All targets #Add StatsAnaTrue laters




TrueEvents: $(MAIN_DIR)/TrueEvents.cpp $(OBJECTS_AsivTrue) 
	@echo " "
	@echo "Generating Application for True Events"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "True event generator Done"

SensitivityTrue: $(MAIN_DIR)/TrueSensitivity.cpp $(OBJECTS_SensitivityTrue) 
	@echo " "
	@echo "Generating application for true sensitivity simulations"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "Application for true sensitivity simulation ready"

EarthSensitivity: $(MAIN_DIR)/EarthTrueSensitivity.cpp $(OBJECTS_AsivTrue) 
	@echo " "
	@echo "Generating Application for True Events"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "True event generator Done"

ObservedEvents: $(MAIN_DIR)/ObservedEvents.cpp $(OBJECTS_AsivObs) 
	@echo " "
	@echo "Generating Application for True Events"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "Observed event generator Done"

SensitivityObs: $(MAIN_DIR)/ObsSensitivity.cpp $(OBJECTS_SensitivityObs) 
	@echo " "
	@echo "Generating application for observed sensitivity simulations"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "Application for observed sensitivity simulation ready"

#StatsAnaTrue: $(MAIN_DIR)/StatsAnaTrue.cpp $(OBJECTS_AsivTrue)
#	@echo " "
#	@echo "Generating Application for Statistical Analysis of True Events"
#	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
#	@echo " "
#	@echo "StatsTrue Done"

#StatsAnaObs: $(MAIN_DIR)/StatsAnaObs.cpp $(OBJECTS_AsivObs)
#	@echo " "
#	@echo "Generating Application for Statistical Analysis of Observed Events"
#	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
#	@echo " "
#	@echo "StatsObs Done"

OscProbEarth: $(MAIN_DIR)/OscProbEarth.cpp $(OBJECTS_OscProb)
	@echo " "
	@echo "Generating Application for Oscillation Probabilities"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "Probabilities Done"

HondaFlux: $(MAIN_DIR)/HondaFlux.cpp $(OBJECTS_Flux) 
	@echo " "
	@echo "Generating Application for Neutrino Flux visualization"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "Neutrino Flux vis generator Done"

Earth3D: $(MAIN_DIR)/Earth3DVis.cpp $(OBJECTS_EARTHMODEL) 
	@echo " "
	@echo "Generating Application for Model visualization"
	$(CC) -g $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@
	@echo " "
	@echo "3D Earth Modelvis generator Done"

#Build Objects

AsimovDataTrue.o: $(SRC_DIR)/AsimovDataTrue.C MathTools.o PhyTools.o 
	@echo "  "
	@echo "Creating Asimov True Events Objects"
	$(CC) -g -c $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@


AsimovDataObs.o: $(SRC_DIR)/AsimovDataObs.C MathTools.o PhyTools.o 
	@echo "  "
	@echo "Creating Asimov Observed Events Objects"
	$(CC) -g -c $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@

OscProbEarth.o: $(SRC_DIR)/OscProbEarth.C MathTools.o PhyTools.o
	@echo "  "
	@echo "Creating Oscillation Probabilities  Objects"
	$(CC) -g -c $^ $(ROOT_FLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS) -o $@

MathTools.o: $(SRC_DIR)/MathTools.C
	@echo "  "
	@echo "Creating MathTools Objects"
	$(CC) -c -g $^ $(ROOT_FLAGS) $(CPPFLAGS) -o $@


PhyTools.o: $(SRC_DIR)/PhyTools.C
	@echo "  "
	@echo "Creating PhyTools Objects"
	$(CC) -c -g $^ $(ROOT_FLAGS) $(CPPFLAGS) -o $@
	

Earth3DModel.o: $(SRC_DIR)/Earth3DModel.C
	@echo "  "
	@echo "Creating a 3D Earth from a Geological models Objects"
	$(CC) -c -g $^ $(ROOT_FLAGS) $(CPPFLAGS) -o $@





clean:
	@echo "Removing Objects"
	rm *.o  ObservedEvents TrueEvents   OscProbEarth HondaFlux EarthSensitivity Earth3D SensitivityTrue SensitivityObs
