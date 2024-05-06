#include "WRF.H"
#include "interpolation.H"

using namespace Foam;

template<typename Type>
Type WRF::interpolate(const point& pt, GeometricField<Type, fvPatchField, volMesh>& psi, const word& interpMethod)
{
  autoPtr<interpolation<Type>> interp
  (
    interpolation<Type>::New(interpMethod, psi)
  );
  point wrf_pt = transform(pt);
  Type interpVal = interp->interpolate
  (
    wrf_pt, 
    searcher_.findCell(wrf_pt)
  );
  return interpVal;
  
}

template<typename Type>
Field<Type> WRF::interpolate(const Field<point>& pts, GeometricField<Type, fvPatchField, volMesh>& psi, const word& interpMethod)
{
  autoPtr<interpolation<Type>> interp
  (
    interpolation<Type>::New(interpMethod, psi)
  );
  pointField wrf_pts = transform(pts);
  Field<Type> interpVals(pts.size());
  for(size_t i=0; i < pts.size(); i++)
  {
    interpVals[i] = interp->interpolate
    (
      wrf_pts[i],
      searcher_.findCell(wrf_pts[i])
    );
  }
  return interpVals;
}
