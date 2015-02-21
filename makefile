FC = gfortran
FFLAGS = -ffast-math #-Wall -Wextra -march=native -O3

LDFLAGS =
LIBS = -llapack

FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS += $(shell pkg-config --libs plplotd-f95)

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

OBJS = 
OBJS += argon_box_dynamics.o
OBJS += argon_box_init.o
OBJS += md_plot.o
OBJS += argon_box.o


all: argon_box

argon_box: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) argon_box $(OBJS) *.mod





