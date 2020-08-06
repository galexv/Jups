CC=gcc
OPTFLAGS:=-g -O3
# OPTFLAGS:=-g -Ofast # to allow "cutting corners" in FP operations
OMPFLAGS:=-fopenmp
CFLAGS=-I. $(OPTFLAGS) $(OMPFLAGS)
LDFLAGS:=-lm $(OMPFLAGS)
DEPS=atmos.h constant.h include.h input.h nrutil.h opac.h
OBJ=geometry.o interpol.o main_rt_no_scat.o planck.o read_t_p_doppler.o readchemtable.o readopactable.o totalopac.o rt_emission_multiphase.o utils.o nrutil.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

rt_emission_aerosols.exe: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

clean:
	rm -f rt_emission_aerosols.exe $(OBJ)
