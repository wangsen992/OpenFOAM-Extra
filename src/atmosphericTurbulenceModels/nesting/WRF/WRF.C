#include "WRF.H"

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
  nc_(ncfile_path, netCDF::NcFile::read),
  runTime_(wrfCaseRoot, wrfCaseName),
  foamTime_(db),
  pmesh_
  (
    fvmeshFromNc(nc_, runTime_)
  ),
  searcher_(pmesh_()),
  ptransformer_(nullptr),
  pitransformer_(nullptr)
{
  Info << "Running WRF init script" << endl;
  // load the transformation
  WRF_PROJ_PARAMS params;
  getProjAtts(nc_, params);
  std::shared_ptr<OGRSpatialReference> pcrs_wrf = getCRS(params);
  std::shared_ptr<OGRSpatialReference> pcrs_foam
  (
    new OGRSpatialReference
  );
  pcrs_foam->SetFromUserInput(foam_proj4.c_str());
  ptransformer_ = OGRCreateCoordinateTransformation(pcrs_foam.get(), pcrs_wrf.get());
  pitransformer_ = ptransformer_->GetInverse();
  pmesh_->movePoints
  (
    itransform(pmesh_->points())
  );
  pmesh_->setInstance(runTime_.constant());
  pmesh_->write();
  Info << "WRF init script complete (constructed & transformed to Foam CRS" << endl;
}

tmp<volVectorField> WRF::U(size_t it)
{
  tmp<volVectorField> pu = load_U(mesh(), nc_, it);
  return pu;
}

tmp<volScalarField> WRF::var(const word& name, size_t it)
{
  Info << "Start loading var" << endl;
  tmp<volScalarField> pvar = load_var(mesh(), nc_, name, it);
  return pvar;
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

point WRF::itransform(const point& pt)
{
  double x(pt.x()), y(pt.y()), z(pt.z());
  pitransformer_->Transform(1, &x, &y);
  return point{x,y,z};
}


pointField WRF::itransform(const pointField& pts)
{
  pointField outPts(pts.size());
  for(int i=0; i<pts.size(); i++)
  {
    outPts[i] = itransform(pts[i]);
  }
  return outPts;
}
