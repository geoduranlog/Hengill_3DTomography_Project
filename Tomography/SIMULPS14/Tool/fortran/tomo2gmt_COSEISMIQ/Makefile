# Makefile for tomo2gmt
# s.h. 99-10-11
#
#  modified for IBM G5 compiler 03-11-04
#

# standard SUN compiler
#F77 = f77
# IBM compiler
#F77 = xlf
F77=gfortran

BINDIR = .

#OFFLAGS= -O4 -std=legacy -ffixed-line-length-none  -static-libgfortran -finit-local-zero

OFFLAGS= -ffixed-line-length-none

# optimization flags for IBM
#IBMFLAGS = -O4 -qarch=ppc970 -qtune=ppc970 -qalign=4k:struct=natural -qfixed=132 -qxlf77=leadzero


# flags for SUN
#SUNFLAGS = -lm


tomo2gmt: tomo2gmt.o 
	${F77} tomo2gmt.o ${OFFLAGS} -o ${BINDIR}/tomo2gmt
	
tomo2gmt.o : tomo2gmt.f tomo2gmt_common.inc
	${F77} -c ${OFFLAGS} tomo2gmt.f
