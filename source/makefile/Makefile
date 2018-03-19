SHELL := /bin/bash

hn ?= $(shell hostname)

ifneq "$(findstring brutus,$(shell hostname))" ""
CPPFLAGS += -D BRUTUS
include make.brutus
endif

ifneq "$(findstring kraken,$(shell hostname))" ""
CPPFLAGS += -D KRAKEN
include make.kraken
endif

ifneq "$(findstring mpp,$(shell hostname))" ""
CPPFLAGS += -D LRZ_CLUSTER
include make.lrz
endif

ifneq "$(findstring lxa,$(shell hostname))" ""
CPPFLAGS += -D LRZ_CLUSTER
include make.lrz
endif

ifeq "$(helmholtz)" "hypre"
	CPPFLAGS += -D_MULTIGRID_
endif

ifeq "$(config)" "production"
CPPFLAGS += -DNDEBUG
CPPFLAGS += -D_FMMSILENT
endif

#set parameters
dim ?= 3 
bpd ?= 8
bs ?= 16 #16
#bpd = 8 for lower resolution 
CPPFLAGS += -D_DIM=$(dim) -D_BLOCKSIZE_=$(bs) -D_BPD_=$(bpd) 
CPPFLAGS += -D HGG

#CPPFLAGS += -D Propagation
#CPPFLAGS += -D Synthetic
#CPPFLAGS += -D HGG_UQ_MALA
#CPPFLAGS += -D Visualisation
#CPPFLAGS += -D Necrosis
#CPPFLAGS += -D VP
#CPPFLAGS += -D Turing


ifeq "$(dim)" "2"
	bsz ?= 1
	CPPFLAGS += -D_BLOCKSIZE_Z_=$(bsz)
endif

ifeq "$(bpd)" "16"
	ml ?= 4
	CPPFLAGS += -D _MAXLEVEL_=$(ml)
endif

ifeq "$(bpd)" "8"
        ml ?= 3
        CPPFLAGS += -D _MAXLEVEL_=$(ml)
endif

	
LIBS += \
	-L$(TBB_LIB_DIR) \
	-ltbb \
	-ltbbmalloc \
	-L$(VTK_LIB_DIR) \
	-lvtkViews \
	-lvtkInfovis \
	-lvtkWidgets \
	-lvtkHybrid \
	-lvtkRendering \
	-lvtkGraphics \
	-lvtkverdict \
	-lvtkImaging \
	-lvtkftgl \
	-lvtkfreetype \
	-lvtkIO \
	-lvtkFiltering \
	-lvtkCommon \
	-lm \
	-lvtkDICOMParser \
	-lvtkmetaio \
	-lvtksqlite \
	-lvtkpng \
	-lvtktiff \
	-lvtkjpeg \
	-lvtkexpat \
	-lvtksys \
	-lvtkexoIIc \
	-lvtkNetCDF \
	-lvtklibxml2 \
	-lvtkzlib \
	-lpthread \
	-ldl \
	$(OPENMP_FLAG)
	
	
##################

ifeq "$(helmholtz)" "hypre"
	CC=mpic++
	LD=$(CC)
	CPPFLAGS += -I$(HYPRE_INC_DIR)
	LIBS += -L$(HYPRE_LIB_DIR) -lHYPRE
endif	


VPATH := ../Glioma/ ../MRAGcore/
.DEFAULT_GOAL := brain

########################################################################
# List of all programs to be compiled
########################################################################
	


#HGG PART
OBJECTS = \
Glioma_main.o \
Glioma_HG_UQ.o\
Glioma_Necrosis.o\
Glioma_ProcessSyntheticData.o\
Glioma_HG_ProcessPatientData.o\
Glioma_HG_Propagation.o\
Glioma_HG_Visualizations.o\
Glioma_dat2VP.o\
MRAGBoundaryBlockInfo.o \
MRAGProfiler.o \
MRAGWavelets_StaticData.o \
Test.o\

#all: brain

brain: $(OBJECTS)
	$(CC) $^ -o $@ $(LIBS)
	@echo done

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o 
	rm -f brain 


