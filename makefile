CC = g++
CCFLAGS = -c --std=c++14 -frtti -O2

all : force.o Body_t.o basic_functions.o main.o
	$(CC) force.o Body_t.o basic_functions.o main.o --std=c++14 -frtti -O2 -o main

main.o: main.cpp
	$(CC) $(CCFLAGS) main.cpp

force.o: force.cpp
	$(CC) $(CCFLAGS) force.cpp

Body_t.o: Body_t.cpp
	$(CC) $(CCFLAGS) Body_t.cpp

basic_functions.o: basic_functions.cpp
	$(CC) $(CCFLAGS) basic_functions.cpp

clean:
	rm -rf *.o main