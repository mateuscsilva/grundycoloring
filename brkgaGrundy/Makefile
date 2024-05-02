# Compiler binary:
CPP= g++

# Recommended compiler flags for speed:
#	OpenMP enabled
#	full binary code optimization
#	full error and warning reports
#	no range checking within BRKGA:
CFLAGS= -O3 -fpermissive -fopenmp -Wextra -Wall -Weffc++ -ansi -pedantic -Woverloaded-virtual -Wcast-align -Wpointer-arith

# Compiler flags for debugging; uncomment if needed:
#	range checking enabled in the BRKGA API
#	OpenMP disabled
#	no binary code optimization
#CFLAGS= -DRANGECHECK -Wextra -Wall -Weffc++ -ansi -pedantic -Woverloaded-virtual -Wcast-align -Wpointer-arith

OBJECTS= GrundyInstance.o GrundySolver.o GrundyDecoder.o GrundyVerifier.o Population.o brkga-grundy.o

# Targets:
all: brkga-grundy

brkga-grundy: $(OBJECTS)
	$(CPP) $(CFLAGS) $(OBJECTS) -o brkga-grundy

brkga-grundy.o:
	$(CPP) $(CFLAGS) -c brkga-grundy.cpp

GrundyInstance.o:
	$(CPP) $(CFLAGS) -c GrundyInstance.cpp

GrundySolver.o:
	$(CPP) $(CFLAGS) -c GrundySolver.cpp

GrundyDecoder.o:
	$(CPP) $(CFLAGS) -c GrundyDecoder.cpp

GrundyVerifier.o:
	$(CPP) $(CFLAGS) -c GrundyVerifier.cpp

Population.o:
	$(CPP) $(CFLAGS) -c brkgaAPI/Population.cpp

# Remove:
clean:
	rm -f brkga-grundy $(OBJECTS)
