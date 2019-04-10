CC = g++
CCFLAGS = -c --std=c++1z -frtti -g -lm -O2

all : force.o Body_t.o basic_functions.o main.o
	$(CC) force.o Body_t.o basic_functions.o main.o --std=c++1z -frtti -O2 -o main

main.o: main.cpp Body_t.hpp basic_types.hpp
	$(CC) $(CCFLAGS) main.cpp

force.o: force.cpp Body_t.hpp 
	$(CC) $(CCFLAGS) force.cpp

Body_t.o: Body_t.cpp Body_t.hpp basic_types.hpp
	$(CC) $(CCFLAGS) Body_t.cpp

basic_functions.o: basic_functions.cpp
	$(CC) $(CCFLAGS) basic_functions.cpp

Body_t.hpp : functions.hpp runge-kutta.hpp

basic_types.hpp : functions.hpp runge-kutta.hpp

clean:
	rm -rf *.o main

result: all
	./main
	