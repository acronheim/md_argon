FC = gfortran
FFLAGS = -ffast-math -march=native -O3 -Wall -Wextra -Wtabs -mno-avx 

LDFLAGS =
LIBS = -llapack -lblas

FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS += $(shell pkg-config --libs plplotd-f95)

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)
MV = mv

OBJDIR = obj/
MODDIR = mod/

OBJS = 
OBJS += argon_box_dynamics.o
OBJS += argon_box_init.o
OBJS += md_plot.o
OBJS += argon_box_results.o
OBJS += argon_box.o
OBJS2 += test.o


all: $(OBJDIR) argon_box test move

$(OBJDIR):
	mkdir -p $(OBJDIR) $(MODDIR)

argon_box: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

test: $(OBJS2)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<

move:           
	$(MV) $(OBJS) $(OBJS2) $(OBJDIR)
	$(MV) *.mod $(MODDIR)

.PHONY: clean
clean:
	$(RM) argon_box test $(OBJS) *.mod 
	rm -rf obj/ mod/





