# VDRC labsource 3Dgeometry Makefile

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

TARGETS= rg_Boundary3D.o rg_NURBSplineSurface3D.o rg_BoundingBox3D.o rg_Outline3D.o rg_Box3D.o rg_PartialCylinder3D.o rg_BSplineCurve3D.o rg_Plane3D.o rg_BSplineSurface3D.o rg_Point3D.o rg_BzCurve3D.o rg_Polyline2D.o rg_BzSurface3D.o rg_Polyline3D.o rg_Cylinder3D.o rg_PrimitiveBox3D.o rg_CylindericalSurface3D.o rg_PrimitiveCylinder3D.o rg_EllipticArc3D.o rg_PrimitiveExtrusion3D.o rg_Extrusion3D.o rg_PrimitiveFrustum3D.o rg_Frustum3D.o rg_RBzCurve3D.o rg_Line3D.o rg_RBzSurface3D.o rg_Link3D.o rg_RevolvedSurface3D.o rg_Loop3D.o rg_Sphere3D.o rg_NUBSplineCurve3D.o rg_ToroidalSurface3D.o rg_NUBSplineSurface3D.o rg_TrimmedNURBSSurface3D.o rg_NURBSplineCurve3D.o rg_TrimmedSurface3D.o

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

