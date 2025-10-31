# Define the name of the final executable
bin_PROGRAMS = rbf_app

# Source files for the executable
rbf_app_SOURCES = \
    main.c \
    alglib_rbf_wrapper.cpp \
    alglib/interpolation.cpp

# Compiler flags
AM_CXXFLAGS = -O2 -std=c++11
AM_CFLAGS = -O2

# Include path if needed
INCLUDES = -I$(srcdir)/alglib

