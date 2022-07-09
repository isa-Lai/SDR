# this makefile is intended for g++ on Linux

CC = g++
CFLAGS = -c -Wall -O3 -pthread
# CFLAGS = -c
LDFLAGS =
INCLUDE = -I../include/
SOURCES = project.cpp iofunc.cpp filter.cpp  fourier.cpp  genfunc.cpp  logfunc.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = project

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) -pthread

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDE) $< -o $@

clean:
	-rm $(OBJECTS)
