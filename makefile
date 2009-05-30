#
#
CFLAGS = -fast -arch ev56
CC = ccc 

SRCS = \
bounds.c diag.c dudp_calc.c dump.c fixup.c gaussj.c grbondi.c \
image.c init.c interp.c lubksb.c ludcmp.c main.c metric.c \
mnewt.c nrutil.c phys.c ranc.c restart.c set_arrays.c set_grid.c \
step_ch.c tensor.c utoprim.c vchar.c zbrent.c tetrad.c \
bltoks.c 
 
OBJS = \
bounds.o diag.o dudp_calc.o dump.o fixup.o gaussj.o grbondi.o \
image.o init.o interp.o lubksb.o ludcmp.o main.o metric.o \
mnewt.o nrutil.o phys.o ranc.o restart.o set_arrays.o set_grid.o \
step_ch.o tensor.o utoprim.o vchar.o zbrent.o tetrad.o \
bltoks.o 

grmhd: $(OBJS) makefile
	$(CC) $(CFLAGS) -o grmhd $(OBJS) \
                $(LAPACKLIB) $(BLASLIB) $(F2CLIB) -lm -lcxml

# dependencies
$(OBJS) : defs.h decs.h nrutil.h makefile

clean:
	rm *.o

OBJB = postmort.o \
bounds.o diag.o dudp_calc.o dump.o fixup.o gaussj.o grbondi.o \
image.o init.o interp.o lubksb.o ludcmp.o metric.o \
mnewt.o nrutil.o phys.o ranc.o restart.o set_arrays.o set_grid.o \
step_ch.o tensor.o utoprim.o vchar.o zbrent.o tetrad.o \
bltoks.o

postmort: $(OBJB) makefile
	$(CC) $(CFLAGS) -o postmort $(OBJB) -lcxml


