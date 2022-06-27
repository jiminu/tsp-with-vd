#include "rg_RBzSurface3D.h"
#include "rg_RelativeOp.h"
#include "rg_MathFunc.h"

rg_RBzSurface3D::rg_RBzSurface3D()
:rg_BzSurface3D()
{
    weight_vector=rg_NULL;
}

rg_RBzSurface3D::rg_RBzSurface3D( const rg_DEGREE& uDegree, const rg_DEGREE& vDegree )
:rg_BzSurface3D(uDegree,vDegree)
{
    weight_vector= new rg_REAL*[uDegree+1];

    for( rg_INT i=0; i <= uDegree; i++ )
    {
        weight_vector[i]= new rg_REAL[vDegree+1];
        for ( rg_INT j=0 ; j <= vDegree; j++ )
        {
            weight_vector[i][j]=1.0;
        }
    }
}

            
rg_RBzSurface3D::rg_RBzSurface3D( const rg_DEGREE& uDegree,
                        const rg_DEGREE& vDegree,
		                const rg_Point3D** ctrlPts )
:rg_BzSurface3D(uDegree,vDegree,ctrlPts)
{
    weight_vector= new rg_REAL*[uDegree+1];

    for( rg_INT i=0; i <= uDegree; i++ )
    {
        weight_vector[i]= new rg_REAL[vDegree+1];
        for ( rg_INT j=0 ; j <= vDegree; j++ )
        {
            weight_vector[i][j]=1.0;
        }
    }
}

rg_RBzSurface3D::rg_RBzSurface3D( const rg_DEGREE& uDegree,
                        const rg_DEGREE& vDegree,
		                const rg_Point3D** ctrlPts ,
                        const rg_REAL**  weightVector)
:rg_BzSurface3D(uDegree,vDegree,ctrlPts)
{
    weight_vector= new rg_REAL*[uDegree+1];

    for( rg_INT i=0; i <= uDegree; i++ )
    {
        weight_vector[i]= new rg_REAL[vDegree+1];
        for ( rg_INT j=0 ; j <= vDegree; j++ )
        {
            weight_vector[i][j]=weightVector[i][j];
        }
    }
}


rg_RBzSurface3D::rg_RBzSurface3D( const rg_RBzSurface3D& surface )
:rg_BzSurface3D(surface)
{
    rg_DEGREE uDegree=rg_BzSurface3D::getDegreeOfU();
    rg_DEGREE vDegree=rg_BzSurface3D::getDegreeOfV();

    weight_vector= new rg_REAL*[uDegree+1];

    for( rg_INT i=0; i <= uDegree; i++ )
    {
        weight_vector[i]= new rg_REAL[vDegree+1];
        for ( rg_INT j=0 ; j <= vDegree; j++ )
        {
            weight_vector[i][j]=surface.weight_vector[i][j];
        }
    }
}

rg_RBzSurface3D::~rg_RBzSurface3D()
{
    rg_DEGREE uDegree=rg_BzSurface3D::getDegreeOfU();

    if ( weight_vector != rg_NULL )
    {
        for( rg_INT i=0 ; i <= uDegree; i++ )
        {
            delete[] weight_vector[i];
        }
        delete[] weight_vector;
    }
}

////  Access elements
rg_REAL   rg_RBzSurface3D::getWeight( const rg_INDEX& i, const rg_INDEX& j ) const
{
    return weight_vector[i][j];
}

//return pointer to the copy of control points.
rg_REAL** rg_RBzSurface3D::getWeightVector() const
{
    rg_DEGREE uDegree=rg_BzSurface3D::getDegreeOfU();
    rg_DEGREE vDegree=rg_BzSurface3D::getDegreeOfV();

	rg_REAL** output = new rg_REAL*[uDegree + 1];

	rg_INT i = 0;
	for( i=0; i <= uDegree; i++ )
    {
		output[i] = new rg_REAL[vDegree + 1];
    }

    for ( i=0; i <= uDegree; i++ )
	{
		for( rg_INDEX j=0; j <= vDegree; j++ )
        {
			output[i][j] = weight_vector[i][j];
        }
	}

    return output;
}



void    rg_RBzSurface3D::setOneWeight( const rg_INDEX& i,
        	                      const rg_INDEX& j,
					              const rg_REAL& weight )
{
    weight_vector[i][j]=weight;
}

void    rg_RBzSurface3D::setWeightVector(const rg_REAL** weights)
{
    rg_DEGREE uDegree=rg_BzSurface3D::getDegreeOfU();
    rg_DEGREE vDegree=rg_BzSurface3D::getDegreeOfV();

    for ( rg_INT i=0; i <= uDegree; i++ )
    {
        for( rg_INT j=0; j <= vDegree; j++ )
        {
            weight_vector[i][j]=weights[i][j];
        }
    }
}    

void    rg_RBzSurface3D::setDegree( const rg_DEGREE& newUDegree,
                               const rg_DEGREE& newVDegree )
{
    rg_BzSurface3D::setDegree(newUDegree,newVDegree);
    weight_vector= new rg_REAL*[newUDegree+1];

    for( rg_INT i=0; i <= newUDegree; i++ )
    {
        weight_vector[i]= new rg_REAL[newVDegree+1];
        for ( rg_INT j=0 ; j <= newVDegree; j++ )
        {
            weight_vector[i][j]=1.0;
        }
    }

}

void    rg_RBzSurface3D::setSurface( const rg_RBzSurface3D& temp)
{
    rg_INT newUDegree=temp.getDegreeOfU();
    rg_INT newVDegree=temp.getDegreeOfV();

    rg_RBzSurface3D::setDegree(newUDegree,newVDegree);
    for ( rg_INT i=0; i < newUDegree+1; i++)
    {
        for ( rg_INT j=0; j < newVDegree+1;j++ )
        {
            rg_Point3D newCtrlPt=temp.getCtrlPt(i,j);
            rg_REAL  newWeight=temp.getWeight(i,j);

            rg_BzSurface3D::setCtrlPt(i,j,newCtrlPt);
            rg_RBzSurface3D::setOneWeight(i,j,newWeight);
        }
    }
}
rg_RBzSurface3D rg_RBzSurface3D::getUHodograph() const
{

    const rg_INT oldUDegree=getDegreeOfU();
    const rg_INT oldVDegree=getDegreeOfV();
    const rg_INT newUDegree=2*oldUDegree;
    const rg_INT newVDegree=2*oldVDegree;
    rg_RBzSurface3D output(newUDegree,newVDegree);

    for (rg_INT i=0; i <= newUDegree; i++ )
    {
        for ( rg_INT j=0; j <= newVDegree; j++ )
        {
            rg_Point3D ctrlPts( 0.0, 0.0, 0.0);
            rg_REAL  weight=0.0;

            for( rg_INT k= rg_Max( 0, i-oldUDegree); k <= rg_Min(oldUDegree,i); k++ )
            {
                for( rg_INT l= rg_Max(0, j-oldVDegree); l <= rg_Min(oldVDegree,j ); l++ )
                {
                    rg_Point3D temp1(0.0,0.0,0.0);
                    rg_Point3D temp2(0.0,0.0,0.0);
                    if ( k != 0 )
                    {
                        temp1 =      getWeight(i-k,j-l) 
                               * (   getWeight(k,l)* getCtrlPt(k,l)  
                                   - getWeight(k-1,l)* getCtrlPt(k-1,l) );
                        temp1 =  temp1
                               - getWeight(i-k,j-l)*getCtrlPt(i-k,j-l)
                                *(  getWeight(k,l) - getWeight(k-1,l)) ;
                        temp1 = k* temp1;

                    }

                    if ( oldUDegree-k != 0 )
                    {
                        temp2 =     getWeight(i-k,j-l)
                                * (   getWeight(k+1,l)* getCtrlPt(k+1,l)  
                                    - getWeight(k,l)* getCtrlPt(k,l) );
                        temp2 =  temp2
                               - getWeight(i-k,j-l)*getCtrlPt(i-k,j-l)
                                * (  getWeight(k+1,l) - getWeight(k,l)) ;
                        temp2 = (oldUDegree-k)* temp2;
                    }
    
                    rg_REAL coefficient =  (rg_REAL) rg_MathFunc::combination( oldUDegree, i-k) * (rg_REAL)rg_MathFunc::combination( oldVDegree, j-l)
                                      * (rg_REAL)rg_MathFunc::combination( oldUDegree,   k) * (rg_REAL)rg_MathFunc::combination( oldVDegree,   l)
                                      / (rg_REAL)rg_MathFunc::combination( newUDegree, i  ) / (rg_REAL)rg_MathFunc::combination( newVDegree, j  );


                    ctrlPts+=  ( temp1+ temp2 )*coefficient;
                    weight+=  getWeight(i-k,j-l)
                            * getWeight(  k,  l)
                            * coefficient;
                }
            }
            output.setOneWeight(i,j,weight);
            output.setCtrlPt(i,j,ctrlPts/weight);
        }
    }

    return output;
}
rg_RBzSurface3D rg_RBzSurface3D::getVHodograph() const
{

    const rg_INT oldUDegree=getDegreeOfU();
    const rg_INT oldVDegree=getDegreeOfV();
    const rg_INT newUDegree=2*oldUDegree;
    const rg_INT newVDegree=2*oldVDegree;
    rg_RBzSurface3D output(newUDegree,newVDegree);

    for (rg_INT i=0; i <= newUDegree; i++ )
    {
        for ( rg_INT j=0; j <= newVDegree; j++ )
        {
            rg_Point3D ctrlPts( 0.0, 0.0, 0.0);
            rg_REAL  weight=0.0;

            for( rg_INT k= rg_Max( 0, i-oldUDegree); k <= rg_Min(oldUDegree,i); k++ )
            {   
                for( rg_INT l= rg_Max(0, j-oldVDegree); l <= rg_Min(oldVDegree,j ); l++ )
                {
                    rg_Point3D temp1(0.0,0.0,0.0);
                    rg_Point3D temp2(0.0,0.0,0.0);
                    if ( l != 0 )
                    {
                        temp1 =  getWeight(i-k,j-l) 
                               * (  getWeight(k,l)* getCtrlPt(k,l)
                                  - getWeight(k,l-1)* getCtrlPt(k,l-1) );
                        temp1 =  temp1
                               - getWeight(i-k,j-l)*getCtrlPt(i-k,j-l)
                                *(  getWeight(k,l) - getWeight(k,l-1)) ;
                        temp1 = l* temp1;

                    }

                    if ( oldVDegree-l != 0 )
                    {
                        temp2 =  getWeight(i-k,j-l)
                                * (  getWeight(k,l+1)* getCtrlPt(k,l+1) 
                                   - getWeight(k,l)* getCtrlPt(k,l)    );
                        temp2 =  temp2
                               - getWeight(i-k,j-l)*getCtrlPt(i-k,j-l)
                                * (  getWeight(k,l+1) - getWeight(k,l)) ;
                        temp2 = (oldVDegree-l)* temp2;
                    }
    
                    rg_REAL coefficient =  (rg_REAL)rg_MathFunc::combination( oldUDegree, i-k) * (rg_REAL)rg_MathFunc::combination( oldVDegree, j-l)
                                      * (rg_REAL)rg_MathFunc::combination( oldUDegree,   k) * (rg_REAL)rg_MathFunc::combination( oldVDegree,   l)
                                      / (rg_REAL)rg_MathFunc::combination( newUDegree, i  ) / (rg_REAL)rg_MathFunc::combination( newVDegree, j  );


                    ctrlPts+=  ( temp1+ temp2 )*coefficient;
                    weight+=  getWeight(i-k,j-l)
                            * getWeight(  k,  l)
                            * coefficient;
                }
            }
            output.setOneWeight(i,j,weight);
            output.setCtrlPt(i,j,ctrlPts/weight);
        }
    }

    return output;
}

//  Operations
rg_Point3D rg_RBzSurface3D::evaluatePt( const rg_PARAMETER& u,
		                      const rg_PARAMETER& v )
{
    rg_DEGREE uDegree=rg_BzSurface3D::getDegreeOfU();
    rg_DEGREE vDegree=rg_BzSurface3D::getDegreeOfV();

    rg_Point3D numerator(0.0,0.0,0.0);
    rg_REAL  denominator=0.0;

    for ( rg_INDEX i=0; i <= uDegree; i++ )
	{
		for ( rg_INDEX j=0; j <= vDegree; j++ )
        {
            rg_REAL basisValue=bernstein(i, uDegree, u)*bernstein(j, vDegree, v);

            numerator+=rg_BzSurface3D::getCtrlPt(i,j)*weight_vector[i][j]*basisValue;
            denominator+=weight_vector[i][j]*basisValue;
        }
	}   
    rg_Point3D output=numerator/denominator;
    return output;
}

rg_RBzSurface3D& rg_RBzSurface3D::operator=(const rg_RBzSurface3D& temp)
{
    if ( this == & temp )
    {
        return *this;
    }

    rg_RBzSurface3D::setSurface(temp);
    
    return *this;
}


