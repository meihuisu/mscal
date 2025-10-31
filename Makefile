
CXX = g++
CC = gcc
CXXFLAGS = -O2 -std=c++11
CFLAGS = -O2

ALGLIB_SRC = alglib/interpolation.cpp
WRAPPER_SRC = alglib_rbf_wrapper.cpp
C_SRC = main.c

OBJS = interpolation.o rbf_wrapper.o main.o

all: rbf_app

interpolation.o: $(ALGLIB_SRC)
    $(CXX) $(CXXFLAGS) -c $< -o $@

rbf_wrapper.o: $(WRAPPER_SRC)
    $(CXX) $(CXXFLAGS) -c $< -o $@

main.o: $(C_SRC)
    $(CC) $(CFLAGS) -c $< -o $@

rbf_app: $(OBJS)
    $(CXX) $(OBJS) -o $@

clean:
    rm -f *.o rbf_app
