# VDRC labsource PolynomialSolver Makefile

CXX= g++
AR=ar rcs

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -I../../include \
#	   -I../../../VoronoiDiagram3D/include 


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH =

LDFLAGS = \
#	-Xlinker --start-group \
#          -L../../LIB \
#           -lvdrc \
#          -L../../../VoronoiDiagram3D/lib \
#            -lVoronoi
#	-Xlinker --end-group


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#



#=================================================
# Define target files (put all your files here)
#=================================================

TARGETS= COMPLEX.o MULLER.o NEWTON.o NULL.o TOOLS.o


all: $(TARGETS)



#=================================================
# Rules
#=================================================
%: %.o
	${CXX} $(LIBPATH) $(LDFLAGS) $< -o $@

.cpp.o:
	${CXX} -c $(CXXFLAGS) $< -o $@


#=================================================
# Clean object files
#=================================================
clean:
	-@test -z "*.o" || rm -f *.o

#=================================================
# Remove executable files
#=================================================
EXEPROGS=$(TARGETS:=$(EXETYPE))

veryclean: clean
	-@test -z "$(EXEPROGS)" || rm -f $(EXEPROGS)

