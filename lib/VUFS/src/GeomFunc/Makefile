# VDRC labsource GeomFunc Makefile

CXX= g++
AR=ar rcs

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -I../include

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = 

LDFLAGS = \
#	-Xlinker --start-group \
#          -L/../LIB \
#            -lvdrc
#	-Xlinker --end-group

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#



#=================================================
# Define target files (put all your files here)
#=================================================

TARGETS= rg_BSplineIntersector.o rg_BzIntersector.o rg_CurveSurfaceFunc.o rg_IntsecImplicitNotUsingMathematica.o rg_GeoFunc.o rg_NUBSplineIntersector.o rg_IntersectFunc.o rg_NURBSCurveIntersector.o rg_IntervalIntersector.o rg_TMatrix2D.o rg_IntsecImplicit3.o rg_TMatrix3D.o rg_IntsecImplicit4.o


all: $(TARGETS)



#=================================================
# Rules
#=================================================
%: %.o
	${CXX} -O $(LIBPATH) $(LDFLAGS) $< -o $@

.cpp.o:
	${CXX} -c -O $(CXXFLAGS) $< -o $@


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

