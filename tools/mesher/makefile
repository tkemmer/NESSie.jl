CXXFLAGS = -Wall -std=c++14

.PHONY: clean all
all: mesher sphere

types.o: types.cpp
	g++ ${CXXFLAGS} -c types.cpp

options.o: options.cpp
	g++ ${CXXFLAGS} -c options.cpp

functions.o: functions.cpp
	g++ ${CXXFLAGS} -c functions.cpp

mesher.o: mesher.cpp
	g++ ${CXXFLAGS} -c mesher.cpp

sphere.o: sphere.cpp
	g++ ${CXXFLAGS} -c sphere.cpp

mesher: types.o options.o functions.o mesher.o
	g++ -lgamer -o mesher types.o options.o functions.o mesher.o

sphere: types.o options.o functions.o sphere.o
	g++ -lgamer -o sphere types.o options.o functions.o sphere.o

clean:
	rm -f *.o mesher sphere
