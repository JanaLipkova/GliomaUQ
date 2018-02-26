SHELL := /bin/bash

CC = g++
OPENMP_FLAG= -fopenmp
		
CPPFLAGS+= -I..
#CPPFLAGS+= -D DoubleSigmoid
#CPPFLAGS+= -D HIERARCHICAL -D PET
#CPPFLAGS+= -D HIERARCHICAL -D MRI_ALL
#CPPFLAGS+= -D HIERARCHICAL -D MRI_T1
#CPPFLAGS+= -D HIERARCHICAL -D MRI_T2

ifeq "$(config)" "production"
CPPFLAGS += -DNDEBUG
endif

##################	


VPATH := ../source/  
.DEFAULT_GOAL := likelihood


OBJECTS = Likelihood_main.o HGG_Likelihood.o

all: likelihood

likelihood: $(OBJECTS)
	$(CC) $^ -o $@ 

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o likelihood

