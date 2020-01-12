SRC := $(wildcard *.cpp)
CXX = g++
CXXFLAGS = -I. -Wall -std=c++17 -lstdc++fs

all: main

main: $(SRC)
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	rm main
