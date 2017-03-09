CXXFLAGS = -Wall -std=c++14

.PHONY: clean
all: mesher

types.o: types.cpp
	g++ ${CXXFLAGS} -c types.cpp

options.o: options.cpp
	g++ ${CXXFLAGS} -c options.cpp

functions.o: functions.cpp
	g++ ${CXXFLAGS} -c functions.cpp

mesher.o: mesher.cpp
	g++ ${CXXFLAGS} -c mesher.cpp

mesher: types.o options.o functions.o mesher.o
	g++ -lgamer -o mesher types.o options.o functions.o mesher.o

clean:
	rm -f *.o mesher
