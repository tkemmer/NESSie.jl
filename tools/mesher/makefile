.PHONY: clean
all: mesher

types.o: types.cpp
	g++ -Wall -std=c++14 -c types.cpp

options.o: options.cpp
	g++ -Wall -std=c++14 -c options.cpp

functions.o: functions.cpp
	g++ -Wall -std=c++14 -c functions.cpp

mesher.o: mesher.cpp
	g++ -Wall -std=c++14 -c mesher.cpp

mesher: types.o options.o functions.o mesher.o
	g++ -Wall -std=c++14 -lgamer -o mesher types.o options.o functions.o mesher.o

clean: 
	rm -f *.o mesher

