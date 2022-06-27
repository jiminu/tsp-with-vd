CXX= g++
LDFLAGS=
CXXFLAGS=
LIB=
INC=-I../include/

#=================================================
# Define target files (put all your files here)
#=================================================
TARGETS= linkedList.o rg_Curve.o rg_Geometry.o rg_Surface.o sortFunc.o


all: $(TARGETS)


#=================================================
# Rules
#=================================================
%: %.o
	${CXX} -O $(LDFLAGS) $< $(LIB) -o $@

.cpp.o:
	${CXX} -c -O $(CXXFLAGS) $(INC) $< -o $@

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

