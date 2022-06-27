# VDRC labsource MathUtils Makefile

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


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#



#=================================================
# Define target files (put all your files here)
#=================================================

TARGETS= rg_CubicPolynomial.o rg_QuadraticImplicitEq.o EigenDecomposition.o rg_ImplicitEqMatrix.o rg_QuarticPolynomial.o nrutil.o rg_ImplicitEquation.o rg_RandomNumber.o rg_MathFunc.o rg_RationalPolynomial.o rg_Matrix.o rg_RelativeOp.o rg_ComplexNumber.o rg_Polynomial.o rg_CubicImplicitEq.o rg_PolynomialWithBound.o rg_LinearPolynomial.o rg_QuadraticPolynomial.o 


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

