PRGM		=test_f
FC		=ifort
CC		=icc
FFLAGS		=-O2
CFLAGS 		=
CPPFLAGS	=
LDFLAGS 	=  
NETCDFDIR	=
INCS		=
LIBS		=
FCOMPILE 	= $(FC) $(FFLAGS)
CCOMPILE 	= $(CC) $(CPPFLAGS) $(CFLAGS)
LINKCC 		= $(FC) $(LDFLAGS)

SRCS := $(wildcard *.f90)
OBJS := $(patsubst %.f90,%.o,$(SRCS))

all :  $(PRGM)

$(PRGM):$(OBJS)
	$(LINKCC) $(OBJS) $(LIBS)  -o $(PRGM)

%.o:%.cpp
	$(CCOMPILE) -c $(INCS) $< -o $@

%.o: %.f90
	$(FCOMPILE) -c $(INCS) $< -o $@

clean:
	rm -f $(OBJS) $(PRGM) 
clean-data:
	rm -f *~ *.txt *.dat *.nc errfile outfile
