SRC := $(wildcard *.cpp)
CXX = g++
CXXFLAGS = -I. -O2 -Wall -pedantic -pthread -std=c++17

all: main

main: $(SRC)
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	rm main
