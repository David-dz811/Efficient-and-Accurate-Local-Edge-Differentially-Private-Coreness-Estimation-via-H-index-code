# CC = g++
# CPPFLAGS = -I. -I/usr/local/include  
# LDFLAGS = -g -fopenmp -L/usr/local/lib  
# LDLIBS = -lgsl -lgslcblas -lm -lgomp  # Link against GSL libraries and libgomp


# DEPS = graph.h utility.h counting.h
# OBJ = graph.o main.o utility.o counting.o

# %.o: %.cpp $(DEPS)
# 	$(CC) -std=c++1y -c -O3 -o $@ $< $(CPPFLAGS) $(LDFLAGS)

# counting: $(OBJ)
# 	$(CC) -std=c++1y -O3 -pthread -o $@ $^ $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)

# clean:
# 	-rm -f counting *.o







# Portable compiler selection: prefer environment's CXX, fallback to g++/clang++
CXX ?= g++
CPPFLAGS = -I. -I/usr/local/include

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	# macOS defaults
	CXX := clang++
	CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -pthread
	LDFLAGS  = -g -L/usr/local/lib
	LDLIBS   = -lm -pthread
else
	# Linux / others
	CXX := $(CXX)
	CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -march=native -mtune=native -DNDEBUG
	LDFLAGS  = -g -L/usr/local/lib
	LDLIBS   = -lm
endif

DEPS = graph.h utility.h counting.h
OBJ  = graph.o main.o utility.o counting.o

all: counting

%.o: %.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

counting: $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	-rm -f counting *.o
