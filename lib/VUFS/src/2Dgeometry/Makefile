# VDRC labsource 2Dgeometry Makefile

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

TARGETS= rg_BoundingBox2D.o rg_Line2D.o rg_Point2D.o rg_RBzCurve2D.o rg_BzCurve2D.o rg_Point2D.o rg_RPolynomialCurve2D.o rg_CBzCurve2D.o rg_PolynomialCurve2D.o rg_RQBzCurve2D.o rg_Circle2D.o rg_QBZCurve2D.o

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

