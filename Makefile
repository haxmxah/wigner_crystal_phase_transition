.SUFFIXES:

FC = gfortran
FCFLAGS = -O3 -Wall -Wextra -J$(MODDIR) -fcheck=all -march=native
FLFLAGS =
MODDIR = modules

PROG = mc_annealing_wigner_crystal
SRC = ${PROG:=.f90}
MOD = ${MODDIR}/config.f90 ${MODDIR}/montecarlo.f90 \
	  ${MODDIR}/physics.f90 ${MODDIR}/io.f90 \
	  ${MODDIR}/analysis.f90
OBJ = ${MOD:.f90=.o} ${SRC:.f90=.o}
ANC = ${OBJ:.o=.anc}

.PHONY: all clean plot run

all: ${PROG}

# Main program compilation

$(PROG): ${OBJ}
	$(FC) $(FCFLAGS) -o $@ ${@:=.o} ${MOD:.f90=.o}

%.anc: %.f90
	$(FC) $(FCFLAGS) -fsyntax-only -c -o $@ $<
	@touch $@

%.o : %.anc
	$(FC) $(FCFLAGS) -c -o $@ $(<:.anc=.f90)
	@touch $@

${MODDIR}/config.anc: ${MODDIR}/config.mod
${MODDIR}/montecarlo.anc: ${MODDIR}/config.anc ${MODDIR}/physics.anc \
	${MODDIR}/montecarlo.mod
${MODDIR}/physics.anc: ${MODDIR}/config.anc ${MODDIR}/physics.mod
${MODDIR}/io.anc: ${MODDIR}/config.anc ${MODDIR}/io.mod
${MODDIR}/analysis.anc: ${MODDIR}/config.anc ${MODDIR}/analysis.mod
${PROG:=.anc}: ${MOD:.f90=.anc}

${MODDIR}/config.mod:
${MODDIR}/montecarlo.mod:
${MODDIR}/physics.mod:
${MODDIR}/io.mod:
${MODDIR}/analysis.mod:

clean:
	rm -f ${PROG} ${OBJ} ${OBJ:.o=.mod} ${ANC} *.o
