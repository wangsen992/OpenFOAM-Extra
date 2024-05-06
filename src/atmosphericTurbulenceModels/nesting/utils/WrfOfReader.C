#include "nesting_utils.H"
#include "fvCFD.H"
#include "GeometricField.H"

using namespace Foam;
std::shared_ptr<volVectorField> load_U(fvMesh& mesh, netCDF::NcFile& dataFile, size_t it)
{
      WrfCaseInfo wrfInfo;
      readWrfCaseInfo(&wrfInfo, dataFile);

      Info << "Retrieving U value at timestep index" << it << endl;
      // Get a variable to look how it behaves
      std::shared_ptr<volVectorField> pVar;
        
      pVar.reset
      (
        new volVectorField 
        (
          IOobject
          (
            "Utmp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
          ),
          mesh,
          dimVelocity
          )
      );

      volVectorField& Var(*pVar);

      float tmp_u, tmp_un,  tmp_v, tmp_vn, tmp_w, tmp_wn;
      netCDF::NcVar U_wrf = dataFile.getVar("U");
      netCDF::NcVar V_wrf = dataFile.getVar("V");
      netCDF::NcVar W_wrf = dataFile.getVar("W");

      int cc = 0;
      // int pcu = 0;
      // int pcv = 0;
      // int pcw = 0;
      for(size_t ibt = 0; ibt < wrfInfo.Ncellz; ibt++)
      {
        for(size_t isn = 0; isn < wrfInfo.Ncelly; isn++)
        {
          for(size_t iwe = 0; iwe < wrfInfo.Ncellx; iwe++)
          {
            cc = iwe + wrfInfo.Ncellx*isn + wrfInfo.Ncellx*wrfInfo.Ncelly*ibt;

            U_wrf.getVar(std::vector<size_t>{it, ibt, isn, iwe}, &tmp_u);
            U_wrf.getVar(std::vector<size_t>{it, ibt, isn, iwe+1}, &tmp_un);
            V_wrf.getVar(std::vector<size_t>{it, ibt, isn, iwe}, &tmp_v);
            V_wrf.getVar(std::vector<size_t>{it, ibt, isn+1, iwe}, &tmp_vn);
            W_wrf.getVar(std::vector<size_t>{it, ibt, isn, iwe}, &tmp_w);
            W_wrf.getVar(std::vector<size_t>{it, ibt+1, isn, iwe}, &tmp_wn);

            Var.primitiveFieldRef()[cc][0] = 0.5*(tmp_u+tmp_un);
            Var.primitiveFieldRef()[cc][1] = 0.5*(tmp_v+tmp_vn);
            Var.primitiveFieldRef()[cc][2] = 0.5*(tmp_w+tmp_wn);
          }
        }
      }
    return pVar;
}

std::shared_ptr<Foam::volScalarField> load_var(Foam::fvMesh& mesh, netCDF::NcFile& dataFile, const std::string& varname, size_t it)
{
      WrfCaseInfo wrfInfo;
      readWrfCaseInfo(&wrfInfo, dataFile);

      Info << "Retrieving U value at timestep index" << it << endl;
      // Get a variable to look how it behaves
      std::shared_ptr<volScalarField> pVar;
        
      pVar.reset
      (
        new volScalarField 
        (
          IOobject
          (
            varname+"tmp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimVelocity
          )
      );

      volScalarField& Var(*pVar);

      float tmp;
      netCDF::NcVar var_wrf = dataFile.getVar(varname);
      Info << "load " << varname << " from wrf" << endl;

      int cc = 0;
      // int pcu = 0;
      // int pcv = 0;
      // int pcw = 0;
      for(size_t ibt = 0; ibt < wrfInfo.Ncellz; ibt++)
      {
        for(size_t isn = 0; isn < wrfInfo.Ncelly; isn++)
        {
          for(size_t iwe = 0; iwe < wrfInfo.Ncellx; iwe++)
          {
            cc = iwe + wrfInfo.Ncellx*isn + wrfInfo.Ncellx*wrfInfo.Ncelly*ibt;

            var_wrf.getVar(std::vector<size_t>{it, ibt, isn, iwe}, &tmp);

            Var.primitiveFieldRef()[cc] = tmp;
          }
        }
      }
    return pVar;
}
