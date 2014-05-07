
ROOTFLAGS=`root-config --cflags`
ROOTLIBDIR=`root-config --libdir`
ROOTLIBFLAGS=`root-config --libs`
INCLUDEFLAGS=-Iinclude

COMP=g++ $(INCLUDEFLAGS) -ggdb -fpermissive

all: bin/clustersim

lib/Palette.o: src/Palette.cpp include/Palette.h
	$(COMP) $(ROOTFLAGS) -c -o lib/Palette.o src/Palette.cpp

lib/sensor.o: src/sensor.cpp include/sensor.h
	$(COMP) $(ROOTFLAGS) -c -o lib/sensor.o src/sensor.cpp 

bin/clustersim: src/clustersim.cpp lib/sensor.o lib/Palette.o
	$(COMP) $(ROOTFLAGS) lib/Palette.o lib/sensor.o src/clustersim.cpp $(ROOTLIBFLAGS) -o bin/clustersim

clean:
	rm -f include/*~ *~ lib/* bin/clustersim

cleandata:
	rm -f *.png *.root

install: bin/clustersim
	cp bin/clustersim ~/bin
