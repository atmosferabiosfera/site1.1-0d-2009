#
#   Makefile for Objective Analysis related programs.
#

#....................................................................

#   The parameters below may vary from machine to machine:
#
#   F77     the Fortran 77 compiler (on most machines it is simply
#           f77, but machines like the Cray or Convex have their
#           own compiler name)
#

F77 = f77

#....................................................................

.SUFFIXES: .F .o

.f.o:
	$(F77) -c $*.f

#....................................................................

site: site.o
	$(F77) -o site site.o

clean:
	rm -rf *.o site
