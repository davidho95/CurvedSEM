MOD=../../mod/
LOBJ=../../obj/
BIN=../../bin/
MV=mv -f
F95=gfortran


FFLAGS =  -I../../mod -J../../mod -O5           \
          -fbounds-check                  \
	  -fbacktrace                     


SRCS =  ../gll_library.f90                               \
        ../lagrange_poly.f90                             \
        OneDimWaveModule.f90 \
        OneDimParameters.f90

OBJS =	$(LOBJ)gll_library.o                          \
		$(LOBJ)lagrange_poly.o   \
		$(LOBJ)OneDimWaveModule.o

MODS =	$(MOD)OneDimWaveModule.mod\

all:    OneDimWaveEqn

OneDimWaveEqn: $(BIN)OneDimWaveEqn $(OBJS)
$(BIN)OneDimWaveEqn: $(OBJS) $(LOBJ)OneDimWaveEqn.o
	$(F95) $(FFLAGS) $(LOBJ)OneDimWaveEqn.o $(OBJS) -o $(BIN)OneDimWaveEqn

$(LOBJ)OneDimWaveEqn.o: OneDimWaveEqn.f90 
	$(F95) $(FFLAGS) OneDimWaveEqn.f90 -c
	$(MV) OneDimWaveEqn.o $(LOBJ)


$(OBJS): $(SRCS)	
	$(F95)   $(FFLAGS)  -c $(SRCS)
	$(MV)  *.o $(LOBJ)
