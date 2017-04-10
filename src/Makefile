
program = mtxuss
subdirs = umfpack pardiso

include $(JIVEDIR)/makefiles/prog.mk
include $(JIVEDIR)/makefiles/packages/*.mk

MY_LIBS 	= umfpack amd suitesparseconfig pardiso500-GNU481-X86-64 lapack openblas gfortran
MY_INCDIRS      = $(subdirs) /opt/include /usr/include/suitesparse
MY_LIBDIRS	= . /opt/lib
MY_CXX_STD_FLAGS    = -DWITH_UMFPACK -DWITH_PARDISO -fopenmp
MY_LD_FLAGS     = -static
