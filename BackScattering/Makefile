PRGM		=main
CC			=icc
CFLAGS		=-O3 -std=c++11
CCOMPILE	=$(CC) $(CFLAGS)
CLINK		=$(CC)


SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp, %.o, $(SRCS))

all : $(PRGM)

$(PRGM): $(OBJS)
	$(CLINK) $(OBJS) -o $(PRGM)

%.o: %.cpp
	$(CCOMPILE) -c $< -o $@

clean:
	rm -f $(OBJS) $(PRGM)
clean-data:
	rm -f *~ *.txt *.dat *.nc errfile outfile
