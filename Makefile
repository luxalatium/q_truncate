TARGET_NAME=qtr_serial

#-----------------------------------------------------------

#C compiler settings

CC=gcc
CXX=g++
LD=ld
AR=ar cru
RANLIB=ranlib
#DEPFLAGS=-MMD
#OPT=-O3 -fopenmp -march=knl

ifndef DEBUG
        NDEBUG=1
        #OPT=-O3 -fopenmp -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512
        OPT=-O3 -fopenmp -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl -mavx512ifma -mavx512vbmi
endif

#CXXFLAGS +=  -Wall -Wno-unused-function -Wfatal-errors $(OPT) -std=c++11 -I/opt/apps/intel18/impi18_0/boost/1.66/include
CXXFLAGS +=  -Wall -Wno-unused-local-typedefs -Wno-sign-compare -Wno-unused-function -Wfatal-errors $(OPT) -std=c++11 -I/opt/apps/gcc7_1/boost/1.64/include
LDFLAGS += -fopenmp

#-----------------------------------------------------------

#MPI settings

ifdef QTRMPI
        CXX=mpicxx
        CXXFLAGS += -DQTRMPI
        TARGET_NAME=qtr_mpi
endif

#-----------------------------------------------------------

#Ensures all code is statically linked on a Linux machine

UNAME := $(shell uname)
MACHINE := $(shell uname -m)

ifeq ($(UNAME), Linux)

        LINUX=1

        ifeq ($(MACHINE), x86_64)
                LINUX64=1
        else
                LINUX32=1
        endif

#MPI often doesn't mix well with static linking...

        ifndef QTRMPI
                #LDFLAGS += -static
        endif
endif

ifeq ($(UNAME), Darwin)
        OSX=1
endif

ifeq ($(UNAME), CYGWIN_NT-5.1)
        WIN32=1
endif

ifdef LINUX
        CXXFLAGS += -DLINUX
endif

ifdef OSX
        CXXFLAGS += -DOSX
endif

ifdef WIN32
        CXXFLAGS += -DWIN32
endif

include Rules.mk
