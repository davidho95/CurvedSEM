MOD=../mod/
LOBJ=../obj/
BIN=../bin/
MV=mv -f
F95=gfortran


FFLAGS =  -I../mod -J../mod -O5           \
          -fbounds-check                  \
	  -fbacktrace                     


SRCS =  gll_library.f90                               \
        lagrange_poly.f90                             \
        OneDimParameters.f90 \
        OneDimWaveModule.f90 \
        TwoDimWaveModule.f90 \
        get_global.f90

OBJS =	$(LOBJ)gll_library.o                          \
		$(LOBJ)lagrange_poly.o   \
		$(LOBJ)OneDimParameters.o \
		$(LOBJ)OneDimWaveModule.o                   \
		$(LOBJ)TwoDimWaveModule.o \
		$(LOBJ)get_global.o

MODS = $(MOD)OneDimParameters.mod \
	   $(MOD)OneDimWaveModule.mod \
	   $(MOD)TwoDimWaveModule.mod

all:    TwoDimWaveEqn

TwoDimWaveEqn: $(BIN)TwoDimWaveEqn $(OBJS)
$(BIN)TwoDimWaveEqn: $(OBJS) $(LOBJ)TwoDimWaveEqn.o
	$(F95) $(FFLAGS) $(LOBJ)TwoDimWaveEqn.o $(OBJS) -o $(BIN)TwoDimWaveEqn

$(LOBJ)TwoDimWaveEqn.o: TwoDimWaveEqn.f90 
	$(F95) $(FFLAGS) TwoDimWaveEqn.f90 -c
	$(MV) TwoDimWaveEqn.o $(LOBJ)


$(OBJS): $(SRCS)	
	$(F95)   $(FFLAGS)  -c $(SRCS)
	$(MV)  *.o $(LOBJ)
