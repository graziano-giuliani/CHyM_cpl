#FC = mpiifort
FC = mpifort
PRESCRIBED = -DNILE # -DAZOV -DBLACKSEA
# PRODUCTION
#FCFLAGS = -O3 -fp-model precise -fp-model source -DCPL $(PRESCRIBED)
FCFLAGS = -O3 -DCPL $(PRESCRIBED)
# DEBUG
#FCFLAGS = -O0 -g -Warn all -check all -traceback -DCPL $(PRESCRIBED)
#FCFLAGS = -O0 -g -Wall -pedantic -fcheck=all -fbacktrace -DCPL $(PRESCRIBED)
CPPFLAGS = `nf-config --fflags`
LIBS = `nf-config --flibs`
APP = chymmain.x

SRC = mod_chym_param.F90 \
      mod_chym_io.F90    \
      mod_chym_model.F90 \
      mod_chym_iface.F90 \
      chymmain.f90

OBJ = $(SRC:.f90=.o)

$(APP): $(OBJ) mod_chym_iface.o 
	$(FC) $(CPPFLAGS) $(FCFLAGS) -o $(APP) $(OBJ) $(LIBS)

%.o: %.F90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.f
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

install: $(APP)

clean:
	rm -f $(APP) *.o *.mod

mod_chym_param.o: mod_chym_param.F90
mod_chym_io.o: mod_chym_io.F90 mod_chym_param.o
mod_chym_model.o: mod_chym_model.F90 mod_chym_param.o
mod_chym_iface.o: mod_chym_iface.F90 mod_chym_param.o mod_chym_io.o mod_chym_model.o
chymmain.o: chymmain.f90 mod_chym_param.o mod_chym_iface.o
