CXX = g++
CXXFLAGS = -g -O2 -std=c++0x

all: imm

imm: imm.cpp tools.cpp models.cpp mt19937ar.c
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean

clean:
	rm -f do
