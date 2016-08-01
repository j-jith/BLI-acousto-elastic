// #include "ElasticRHS.h"

template <int dim>
ElasticRHS<dim>::ElasticRHS () : Function<dim> (dim)
{
}

template <int dim> 
inline void ElasticRHS<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const
{
   Assert(values.size()==dim, ExcDimensionMismatch(values.size(),dim));
   Assert(dim>=2,ExcNotImplemented());

// Point<dim> p1(0.123,0.008);
//
// if((p-p1).square()<0.001)
//    values(1)=1;
// else
//    values(1)=0;
   
   for(int j=0; j<dim-1; j++)
      values(j)=0;

   values(dim-1)=9.8e3;

}

template <int dim> 
void ElasticRHS<dim>::vector_value_list(const std::vector<Point<dim> > &points, std::vector<Vector<double> > &value_list) const
{
   Assert(value_list.size()==points.size(),ExcDimensionMismatch(value_list.size(),points.size()));

   const unsigned int n_points = points.size();

   for (unsigned int p=0; p<n_points; ++p)
      ElasticRHS<dim>::vector_value(points[p],value_list[p]);
}
