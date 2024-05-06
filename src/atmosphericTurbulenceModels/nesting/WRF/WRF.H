#ifndef WRF_H
#define WRF_H

#include "fvMesh.H"
#include "meshSearch.H"
#include "nesting_utils.H"

namespace Foam
{
  class WRF
  :
    public IOobject
  {
    // Private member attributes
    netCDF::NcFile nc_;
    Time runTime_; // Time registry for the main of case
    const Time& foamTime_;
    fvMesh& mesh_;
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
      
      inline netCDF::NcFile& file(){return nc_;};

      inline fvMesh& mesh(){return mesh_;};

      inline Time& time(){return runTime_;};

      inline const Time& foamTime(){return foamTime_;};

      template<typename Type>
      inline bool hasVar(const string& name)
      {
        return mesh_.foundObject<Type>(name);
      };
    
    // Coordinate transformation
      // Transform one point from input to WRF
      point transform(const point&);

      // Transform pointField from input to WRF
      pointField transform(const pointField&);

    // Interpolation functions
      // interpolate one point from source field
      template<typename Type>
      Type interpolate(const point&, GeometricField<Type, fvPatchField, volMesh>&, const word& interpMethod="pointMVC");

      // interpolate multiple points from source field
      template<typename Type>
      Field<Type> interpolate(const Field<point>&, GeometricField<Type, fvPatchField, volMesh>&, const word& interpMethod="pointMVC");
  };
}

#endif
