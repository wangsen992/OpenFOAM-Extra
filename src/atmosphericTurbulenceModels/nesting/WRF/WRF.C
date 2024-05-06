#include "WRF.H"
#include "interpolation.H"

using namespace Foam;

WRF::WRF
(
  const string& ncfile_path, 
  const string& wrfCaseRoot, 
  const string& wrfCaseName, 
  const Time& db,
  const string& foam_proj4
)
:
  IOobject("WRF", wrfCaseRoot, db),
  nc_(ncfile_path, netCDF::NcFile::read),
  runTime_(wrfCaseRoot, wrfCaseName),
  foamTime_(db),
  mesh_
  (
    fvmeshFromNc(nc_, runTime_)()
  ),
  searcher_(mesh_),
  ptransformer_(nullptr)
{
  Info << "Running WRF init script" << endl;
  // load the transformation
  WRF_PROJ_PARAMS params;
  getProjAtts(nc_, params);
  std::shared_ptr<OGRSpatialReference> pcrs_wrf = getCRS(params);
  std::shared_ptr<OGRSpatialReference> pcrs_foam;
  pcrs_foam->SetFromUserInput(foam_proj4.c_str());
  ptransformer_ = OGRCreateCoordinateTransformation(pcrs_foam.get(), pcrs_wrf.get());
  Info << "WRF init script complete" << endl;
}

point WRF::transform(const point& pt)
{
  double x(pt.x()), y(pt.y()), z(pt.z());
  ptransformer_->Transform(1, &x, &y);
  return point{x,y,z};
}


pointField WRF::transform(const pointField& pts)
{
  pointField outPts(pts.size());
  for(int i=0; i<pts.size(); i++)
  {
    outPts[i] = transform(pts[i]);
  }
  return outPts;
}

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
