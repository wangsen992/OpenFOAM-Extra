#include "WrfOfReader.H"

WrfOfReader::WrfOfReader(fvMesh& mesh, netCDF::NcFile& dataFile)
:
  mesh_(mesh),
  dataFile_(dataFile)
{
  readWrfCaseInfo(&wrfInfo_, dataFile_);
}


autoPtr<volVectorField> WrfOfReader::load_U(int nt)
{

      Info << "Retrieving U value at timestep " << nt << endl;
      // Get a variable to look how it behaves
      autoPtr<volVectorField> pVar;
        
      pVar.set
      (
        new volVectorField 
        (
          IOobject
          (
            "Utmp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
          ),
          mesh_,
          dimVelocity
          )
      );

      volVectorField& Var(pVar());

      Info << "tmp_var init" << endl;
      float* tmp_var_u = new float[wrfInfo_.N_Time * wrfInfo_.N_bottom_top * wrfInfo_.N_south_north * wrfInfo_.N_west_east_stag];
      float* tmp_var_v = new float[wrfInfo_.N_Time * wrfInfo_.N_bottom_top * wrfInfo_.N_south_north_stag * wrfInfo_.N_west_east];
      float* tmp_var_w = new float[wrfInfo_.N_Time * wrfInfo_.N_bottom_top_stag * wrfInfo_.N_south_north * wrfInfo_.N_west_east];
      netCDF::NcVar U_wrf = dataFile_.getVar("U");
      netCDF::NcVar V_wrf = dataFile_.getVar("V");
      netCDF::NcVar W_wrf = dataFile_.getVar("W");
      Info << "load U_wrf" << endl;
      U_wrf.getVar(tmp_var_u);
      V_wrf.getVar(tmp_var_v);
      W_wrf.getVar(tmp_var_w);

      int cc = 0;
      int pcu = 0;
      int pcv = 0;
      int pcw = 0;
      for(int ibt = 0; ibt < wrfInfo_.Ncellz; ibt++)
      {
        for(int isn = 0; isn < wrfInfo_.Ncelly; isn++)
        {
          for(int iwe = 0; iwe < wrfInfo_.Ncellx; iwe++)
          {
            // The commented line produced a figure that is wrong
            // order in xy
            // This is probably due to the way data is stored 
            // in the netcdf file. 
            // cc = iwe + Ncellx*isn + Ncellx*Ncelly*ibt;
            cc = iwe + wrfInfo_.Ncellx*isn + wrfInfo_.Ncellx*wrfInfo_.Ncelly*ibt;
            pcu = nt * wrfInfo_.N_bottom_top * wrfInfo_.N_south_north *wrfInfo_.N_west_east_stag
                + ibt * wrfInfo_.N_south_north *wrfInfo_.N_west_east_stag
                + isn * wrfInfo_.N_west_east_stag + iwe;
            pcv = nt * wrfInfo_.N_bottom_top * wrfInfo_.N_south_north_stag *wrfInfo_.N_west_east
                + ibt * wrfInfo_.N_south_north_stag *wrfInfo_.N_west_east
                + isn * wrfInfo_.N_west_east + iwe;
            pcw = nt * wrfInfo_.N_bottom_top_stag * wrfInfo_.N_south_north *wrfInfo_.N_west_east
                + ibt * wrfInfo_.N_south_north *wrfInfo_.N_west_east
                + isn * wrfInfo_.N_west_east + iwe;

            Var.primitiveFieldRef()[cc][0] = 
                *(tmp_var_u + pcu);
            Var.primitiveFieldRef()[cc][1] = 
                *(tmp_var_v + pcv);
            Var.primitiveFieldRef()[cc][2] = 
                *(tmp_var_w + pcw);
          }
        }
      }
      delete [] tmp_var_u;
      tmp_var_u = NULL;
      delete [] tmp_var_v;
      tmp_var_v = NULL;
      delete [] tmp_var_w;
      tmp_var_w = NULL;
    return pVar;
}
