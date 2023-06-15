#Variable names
CC= g++
CFLAGS=`root-config  --cflags`
LDFLAGS=`root-config --glibs `
ROOTLIBS = `root-config --libs `

ROOT_FLAGS = `root-config  --cflags --glibs --libs`

CFILES = MathTools.C PhyTools.C  GetObservedEvents.C GetTrueEvents.C

OBJECTStrue = MathTools.o PhyTools.o GetTrueEvents.o
OBJECTSobs  = MathTools.o PhyTools.o GetObservedEvents.o 
OSCPROB_PATH = /home/dehy0499/Desktop/Tomography/Simulation/OscProb

ROOTMATH_LIB=/home/dehy0499/root/lib
ROOTMATH_INCLUDE=/home/dehy0499/root/include/Math

all: ObservedEvents TrueEvents


ObservedEvents: ObservedEvents.cpp $(OBJECTSobs)
	@echo " "
	@echo "Generating Application for Observed Events"
	@ #$(CC)  $^ $(CFLAGS) $(LDFLAGS) $(ROOTLIBS) -L$(OSCPROB_PATH) -lOscProb -I$(OSCPROB_PATH) -L$(ROOTMATH_LIB) -lMathMore -I$(ROOTMATH_INCLUDE) -o $@
	$(CC) -g $^ $(ROOT_FLAGS) -L$(OSCPROB_PATH) -lOscProb -I$(OSCPROB_PATH) -L$(ROOTMATH_LIB) -lMathMore -I$(ROOTMATH_INCLUDE) -o $@
	@echo " "
	@echo "Done"

TrueEvents: TrueEvents.cpp $(OBJECTStrue)
	@echo " "
	@echo "Generating Application for True Events"
	@ #$(CC)  $^ $(CFLAGS) $(LDFLAGS) $(ROOTLIBS) -L$(OSCPROB_PATH) -lOscProb -I$(OSCPROB_PATH) -L$(ROOTMATH_LIB) -lMathMore -I$(ROOTMATH_INCLUDE) -o $@
	$(CC) -g $^ $(ROOT_FLAGS) -L$(OSCPROB_PATH) -lOscProb -I$(OSCPROB_PATH) -L$(ROOTMATH_LIB) -lMathMore -I$(ROOTMATH_INCLUDE) -o $@
	@echo " "
	@echo "Done"



GetObservedEvents.o: GetObservedEvents.C MathTools.o PhyTools.o
	@echo "  "
	@echo "Creating Observed Events Generator Objects"
	$(CC)  -c $^ $(ROOT_FLAGS) -L$(OSCPROB_PATH) -lOscProb -I$(OSCPROB_PATH) -L$(ROOTMATH_LIB) -lMathMore -I$(ROOTMATH_INCLUDE) -o $@


GetTrueEvents.o: GetTrueEvents.C MathTools.o PhyTools.o
	@echo "  "
	@echo "Creating True Events Objects"
	$(CC)  -c $^ $(ROOT_FLAGS) -L$(OSCPROB_PATH) -lOscProb -I$(OSCPROB_PATH) -L$(ROOTMATH_LIB) -lMathMore -I$(ROOTMATH_INCLUDE) -o $@

## SYNTAX:  "@" indicates the parameters before ":", "@" indicates the parameters after ":"
MathTools.o: MathTools.C
	@echo "  "
	@echo "Creating MathTools Objects"
	@ #$(CC) -c $^ $(CFLAGS) $(LDFLAGS) $(ROOTLIBS) -o $@
	$(CC) -c $^ $(ROOT_FLAGS) -o $@


PhyTools.o: PhyTools.C
	@echo "  "
	@echo "Creating PhyTools Objects"
	@ #$(CC) -c $^ $(CFLAGS) $(LDFLAGS) $(ROOTLIBS) -o $@
	$(CC) -c $^ $(ROOT_FLAGS) -o $@





clean:
	@echo "Removing Objects"
	rm *.o ObservedEvents TrueEvents

