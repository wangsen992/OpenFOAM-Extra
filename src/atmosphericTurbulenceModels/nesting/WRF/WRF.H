#ifndef WRF_H
#define WRF_H

#include "fvMesh.H"
#include "GeometricField.H"
#include "meshSearch.H"
#include "nesting_utils.H"

namespace Foam
{
  class WRF
  {
    // Private member attributes
    netCDF::NcFile nc_;
    Time runTime_; // Time registry for the main of case
    const Time& foamTime_;
    autoPtr<fvMesh> pmesh_;
    meshSearch searcher_;
    OGRCoordinateTransformation* ptransformer_;

  public:

    WRF
    (
      const string& ncfile_path, 
      const string& wrfCaseRoot, 
      const string& wrfCaseName, 
      const Time& db,
      const string& foam_proj4
    );

    ~WRF(){};

    // Member functions for access
      
      netCDF::NcFile& file(){return nc_;};

      autoPtr<fvMesh> pmesh(){return pmesh_;};

      fvMesh& mesh(){return pmesh_();};

      Time& time(){return runTime_;};

      const Time& foamTime(){return foamTime_;};

      template<typename Type>
      inline bool hasVar(const string& name)
      {
        return pmesh_->foundObject<Type>(name);
      };

    // Load variables from netcdf file to mesh

      tmp<volVectorField> U(size_t it);

      tmp<volScalarField> var(const word& name, size_t it);

    // Coordinate transformation
      // Transform one point from input to WRF
      point transform(const point&);

      // Transform pointField from input to WRF
      pointField transform(const pointField&);

    // Interpolation functions
    // interpolate one point from source field
       template<typename Type>
       Type interpolate(const point& pt, GeometricField<Type, fvPatchField, volMesh>& psi, const word& interpMethod="pointMVC");
 
       // interpolate multiple points from source field
       template<typename Type>
       Field<Type> interpolate(const Field<point>& pts, GeometricField<Type, fvPatchField, volMesh>& psi, const word& interpMethod="pointMVC");
 
   };
}
#include "WRFI.C"

#endif
