CXX=g++
CPPFLAGS= -O3

all: 
	g++ -O3 -o alleluia main.cpp

clean:
	rm -f alleluia

