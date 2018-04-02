SHELL=/bin/tcsh
FFLAGS = -g +es
INCLUDES = -I. -I/usr/include
LIBS =

FF= g77

#Rule to create .o files from .f files
.f.o:
	 $(FF) -c $(FFLAGS) $(INCLUDES) $*.f

smpsi : smpsi.o physput.o chmsput.o ressput.o thermev.o \
	heatdif.o bshevi.o
	$(FF) smpsi.o  physput.o chmsput.o ressput.o thermev.o \
	heatdif.o bshevi.o $(LIBS) -o smpsi

smpsi.o: smpsi.f physput.f chmsput.f ressput.f thermev.f \
        heatdif.f bshevi.f
physput.o: physput.f
chmsput.o: chmsput.f
ressput.o: ressput.f
thermev.o: thermev.f
heatdif.o: heatdif.f
bshevi.o: bshevi.f

clean : 
	/bin/rm *.o smpsi

run :
	./smpsi
