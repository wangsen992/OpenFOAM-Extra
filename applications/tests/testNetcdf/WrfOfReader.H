#ifndef WrfOfReader_H
#define WrfOfReader_H

#include <netcdf>
#include "fvCFD.H"
#include "WrfTools.H"

class WrfOfReader
{
      
      fvMesh& mesh_;
      netCDF::NcFile& dataFile_;
      struct WrfCaseInfo wrfInfo_;
    
  public:
      
      WrfOfReader(fvMesh&, netCDF::NcFile&);

      autoPtr<volVectorField> load_U(int);

      autoPtr<volVectorField> load_theta();

      autoPtr<volScalarField> load_var(word varName);
};


#endif
