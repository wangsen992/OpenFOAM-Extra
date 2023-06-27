/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    testMesh

Description
    Testing extracting netcdf file for openfoam loading

    Key here is to map from the coordinate system in WRF (or lat lon based
    dataset) to 

    Conversion referenced from wrf-python utils

\*---------------------------------------------------------------------------*/
#include "LambertConverter.H"

#include "fvCFD.H"
#include "IOobject.H"
#include "DynamicList.H"
#include "Ostream.H"
#include "meshSearch.H"
#include "pointField.H"


#include <netcdf>
#include <list>
#include <map>
#include <string>
#include <math.h>

#define DATESTRLEN 19
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void printTimes(char times[][DATESTRLEN])
{
        for ( int i = 0; i < 24; i++)
        {
            std::cout << i << ": ";
            for ( int j = 0; j < 19; j++)
            {
                std::cout << times[i][j];
            }
            std::cout << std::endl;
        }
}

int main(int argc, char *argv[])
{
  
    argList::addOption("filepath", "/path/to/file");
    argList::addBoolOption("mesh");
    argList args(argc, argv);
    Time runTime(Time::controlDictName, args);

    // Loading the input netcdf file into memory (from raw wrf output format)
    const char* filename = args.optionRead<Foam::string>("filepath").c_str();
    std::cout << "input filename : " << filename << std::endl;

    netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
    std::cout << "nc name: " << dataFile.getName() << "\n"
              << std::endl;

    // Start parsing the coordinate systems of the result
      // Collect dimension sizes
      int N_Time(dataFile.getDim("Time").getSize());
      int N_west_east(dataFile.getDim("west_east").getSize());
      int N_west_east_stag(dataFile.getDim("west_east_stag").getSize());
      int N_south_north(dataFile.getDim("south_north").getSize());
      int N_south_north_stag(dataFile.getDim("south_north_stag").getSize());
      int N_bottom_top(dataFile.getDim("bottom_top").getSize());
      int N_bottom_top_stag(dataFile.getDim("bottom_top_stag").getSize());

      // Load attributes
      int MAP_PROJ;
      float CEN_LAT, CEN_LON, TRUELAT1, TRUELAT2, MOAD_CEN_LAT, STAND_LON;
      float POLE_LAT, POLE_LON, DX, DY;
      dataFile.getAtt("MAP_PROJ").getValues(&MAP_PROJ);
      dataFile.getAtt("CEN_LAT").getValues(&CEN_LAT);
      dataFile.getAtt("CEN_LON").getValues(&CEN_LON);
      dataFile.getAtt("TRUELAT1").getValues(&TRUELAT1);
      dataFile.getAtt("TRUELAT2").getValues(&TRUELAT2);
      dataFile.getAtt("MOAD_CEN_LAT").getValues(&MOAD_CEN_LAT);
      dataFile.getAtt("STAND_LON").getValues(&STAND_LON);
      dataFile.getAtt("POLE_LAT").getValues(&POLE_LAT);
      dataFile.getAtt("POLE_LON").getValues(&POLE_LON);
      dataFile.getAtt("DX").getValues(&DX);
      dataFile.getAtt("DY").getValues(&DY);

      mapProj::LambertConverter lam
      (
        STAND_LON,
        TRUELAT1,
        TRUELAT2,
        6378137,  // Earth radius
        1,        // hemi
        1,        // dx
        1,        // dy
        CEN_LAT,  // lat1
        CEN_LON,  // lon1
        0,        // x1
        0         // y1
      );

    int Nx(N_west_east_stag-2);
    int Ny(N_south_north_stag-2);
    int Nz(N_bottom_top_stag-2);
    int Ncellx(Nx-1);
    int Ncelly(Ny-1);
    int Ncellz(Nz-1);
    if (false)
    { 
      // Loading coordinates from lat/lon to xy
      float tmp_full[N_Time][N_south_north_stag][N_west_east_stag];
      // Load lat and lon (note only first time step is needed)
      float X[N_south_north][N_west_east];
      float Y[N_south_north][N_west_east];

      {
        float xlat[N_south_north][N_west_east];
        netCDF::NcVar XLAT(dataFile.getVar("XLAT"));
        XLAT.getVar(tmp_full);
        for(int isn=0; isn < N_south_north; isn++)
        {
            for(int iwe=0; iwe < N_west_east; iwe++)
            {
                xlat[isn][iwe] = tmp_full[0][isn][iwe];
            }
        }

        float xlon[N_south_north][N_west_east];
        netCDF::NcVar XLON(dataFile.getVar("XLONG"));
        XLON.getVar(tmp_full);
        for(int isn=0; isn < N_south_north; isn++)
        {
            for(int iwe=0; iwe < N_west_east; iwe++)
            {
                xlon[isn][iwe] = tmp_full[0][isn][iwe];
            }
        }

        for(int isn=0; isn < N_south_north; isn++)
        {
            for(int iwe=0; iwe < N_west_east; iwe++)
            {
                lam.get_xy
                (
                  xlat[isn][iwe],xlon[isn][iwe],
                  X[isn][iwe], Y[isn][iwe]
                );
                // std::cout << "X, Y" << ": " << X[isn][iwe] << ", " << Y[isn][iwe] << std::endl;
            }
        }
      }

      // Load lat and lon (U_stag, note only first time step is needed)
      float X_U[N_south_north][N_west_east_stag];
      float Y_U[N_south_north][N_west_east_stag];
      {
        float xlat_u[N_south_north][N_west_east_stag];
        // float tmp_full_U_stag[N_Time][N_south_north][N_west_east_stag];
        netCDF::NcVar XLAT_U(dataFile.getVar("XLAT_U"));
        XLAT_U.getVar(tmp_full);
        for(int isn=0; isn < N_south_north; isn++)
        {
            for(int iwe=0; iwe < N_west_east_stag; iwe++)
            {
                xlat_u[isn][iwe] = tmp_full[0][isn][iwe];
            }
        }

        float xlon_u[N_south_north][N_west_east_stag];
        netCDF::NcVar XLON_U(dataFile.getVar("XLONG_U"));
        XLON_U.getVar(tmp_full);
        for(int isn=0; isn < N_south_north; isn++)
        {
            for(int iwe=0; iwe < N_west_east_stag; iwe++)
            {
                xlon_u[isn][iwe] = tmp_full[0][isn][iwe];
            }
        }

        for(int isn=0; isn < N_south_north; isn++)
        {
            for(int iwe=0; iwe < N_west_east_stag; iwe++)
            {
                lam.get_xy
                (
                  xlat_u[isn][iwe],xlon_u[isn][iwe],
                  X_U[isn][iwe], Y_U[isn][iwe]
                );
                // std::cout << "X, Y" << ": " << X[isn][iwe] << ", " << Y[isn][iwe] << std::endl;
            }
        }
      }
      
      // Load lat and lon (V_stag, note only first time step is needed)
      // float tmp_full_V_stag[N_Time][N_south_north_stag][N_west_east];
      float X_V[N_south_north_stag][N_west_east];
      float Y_V[N_south_north_stag][N_west_east];
      {
        float xlat_v[N_south_north_stag][N_west_east];
        netCDF::NcVar XLAT_V(dataFile.getVar("XLAT_V"));
        XLAT_V.getVar(tmp_full);
        for(int isn=0; isn < N_south_north_stag; isn++)
        {
            for(int iwe=0; iwe < N_west_east; iwe++)
            {
                xlat_v[isn][iwe] = tmp_full[0][isn][iwe];
            }
        }

        float xlon_v[N_south_north_stag][N_west_east];
        netCDF::NcVar XLON_V(dataFile.getVar("XLONG_V"));
        XLON_V.getVar(tmp_full);
        for(int isn=0; isn < N_south_north_stag; isn++)
        {
            for(int iwe=0; iwe < N_west_east; iwe++)
            {
                xlon_v[isn][iwe] = tmp_full[0][isn][iwe];
            }
        }

        for(int isn=0; isn < N_south_north_stag; isn++)
        {
            for(int iwe=0; iwe < N_west_east; iwe++)
            {
                lam.get_xy
                (
                  xlat_v[isn][iwe],xlon_v[isn][iwe],
                  X_V[isn][iwe], Y_V[isn][iwe]
                );
            }
        }
      }
      
      // Load z (need to assign real value)
      float Z[N_bottom_top];
      float Z_stag[N_bottom_top_stag];
      {
        float Z0(0);
        float dz(1000);
        for (int i=0; i < N_bottom_top_stag; i++)
        {
            Z_stag[i] = Z0 + dz*(i-1);
        }
        for (int i=0; i < N_bottom_top; i++)
        {
            Z[i] = 0.5 * (Z_stag[i] + Z_stag[i+1]);
        }
      }
     
      // Construct OpenFOAM mesh from the given xy points
      // Note: the outer most columes and rows are ignored, and 
      //       the mesh is slightly distored due to two stagged grid are
      //       combined directly
       
        // Construct points using stagged grid points
        int N_points(Nx*Ny*Nz);
        pointField points(N_points);
        int pc=0;
        for(int i=1; i < N_west_east_stag-1; i++)
        {
          for(int j=1; j < N_south_north_stag-1; j++)
          {
            for(int k=1; k < N_bottom_top_stag-1; k++)
            {
                pc = (i-1) + Nx*(j-1) + Nx*Ny*(k-1);
                points[pc] = Foam::vector(X_U[i][j], Y_V[i][j], Z_stag[k]);
            }
          }
        }
        
        // Construct Cells
        List<FixedList<label, 8>> cells(Ncellx*Ncelly*Ncellz);
        int cc=0;
        for(int i=0; i < Ncellx; i++)
        {
          for(int j=0; j < Ncelly; j++)
          {
            for(int k=0; k < Ncellz; k++)
            {
                cc= i + Ncellx*j + Ncellx*Ncelly*k;
                cells[cc][0]= i + Nx*j + Nx*Ny*k;
                cells[cc][1]= (i+1) + Nx*j + Nx*Ny*k;
                cells[cc][2]= (i+1) + Nx*(j+1) + Nx*Ny*k;
                cells[cc][3]= i + Nx*(j+1) + Nx*Ny*k;
                cells[cc][4]= i + Nx*j + Nx*Ny*(k+1);
                cells[cc][5]= (i+1) + Nx*j + Nx*Ny*(k+1);
                cells[cc][6]= (i+1) + Nx*(j+1) + Nx*Ny*(k+1);
                cells[cc][7]= i + Nx*(j+1) + Nx*Ny*(k+1);
            }
          }
        }
        const cellModel& hex = *(cellModeller::lookup("hex"));
        cellShapeList cellShapes(cells.size());
        for (int celli = 0; celli < cells.size(); celli++)
        {
            cellShapes[celli] = cellShape(hex, labelList(cells[celli]), true);
        }

        // Now set the patches (borrwo from createBlock.C)
        FixedList<List<FixedList<label, 4>>, 6> patches;
        label patchi = 0;
        label facei = 0;
        // x-direction

        // x-min
        patches[patchi].setSize(Ncelly*Ncellz);
        for (label k=0; k<Ncellz; k++)
        {
            for (label j=0; j<Ncelly; j++)
            {
                patches[patchi][facei][0] = 0 + Nx*j + Nx*Ny*k;
                patches[patchi][facei][1] = 0 + Nx*j + Nx*Ny*(k+1);
                patches[patchi][facei][2] = 0 + Nx*(j+1) + Nx*Ny*(k+1);
                patches[patchi][facei][3] = 0 + Nx*(j+1) + Nx*Ny*k;


                facei++;
            }
        }

        // x-max
        patchi++;
        facei = 0;

        patches[patchi].setSize(Ncelly*Ncellz);

        for (label k=0; k<Ncellz; k++)
        {
            for (label j=0; j<Ncelly; j++)
            {
                patches[patchi][facei][0] = Ncellx + Nx*j + Nx*Ny*k;
                patches[patchi][facei][1] = Ncellx + Nx*(j+1) + Nx*Ny*k;
                patches[patchi][facei][2] = Ncellx + Nx*(j+1) + Nx*Ny*(k+1);
                patches[patchi][facei][3] = Ncellx + Nx*j + Nx*Ny*(k+1);
                facei++;
            }
        }

        // y-direction

        // y-min
        patchi++;
        facei = 0;

        patches[patchi].setSize(Ncellx*Ncellz);
        for (label i=0; i<Ncellx; i++)
        {
            for (label k=0; k<Ncellz; k++)
            {
                patches[patchi][facei][0] = i + 0 + Nx*Ny*k;
                patches[patchi][facei][1] = i+1 + 0 + Nx*Ny*k;
                patches[patchi][facei][2] = i+1 + 0 + Nx*Ny*(k+1);
                patches[patchi][facei][3] = i + 0 + Nx*Ny*(k+1);
                facei++;
            }
        }

        // y-max
        patchi++;
        facei = 0;

        patches[patchi].setSize(Ncellx*Ncellz);

        for (label i=0; i<Ncellx; i++)
        {
            for (label k=0; k<Ncellz; k++)
            {
                patches[patchi][facei][0] = i   + Nx*Ncelly + Nx*Ny*k;
                patches[patchi][facei][1] = i + Nx*Ncelly + Nx*Ny*(k+1);
                patches[patchi][facei][2] = i+1 + Nx*Ncelly + Nx*Ny*(k+1);
                patches[patchi][facei][3] = i+1   + Nx*Ncelly + Nx*Ny*k;
                facei++;
            }
        }

        // z-direction

        // z-min
        patchi++;
        facei = 0;

        patches[patchi].setSize(Ncellx*Ncelly);

        for (label i=0; i<Ncellx; i++)
        {
            for (label j=0; j<Ncelly; j++)
            {
                patches[patchi][facei][0] = i + Nx*j + 0;
                patches[patchi][facei][1] = i + Nx*(j+1) + 0;
                patches[patchi][facei][2] = (i+1) + Nx*(j+1) + 0;
                patches[patchi][facei][3] = (i+1) + Nx*j + 0;
                facei++;
            }
        }

        // z-max
        patchi++;
        facei = 0;

        patches[patchi].setSize(Ncellx*Ncelly);

        for (label i=0; i<Ncellx; i++)
        {
            for (label j=0; j<Ncelly; j++)
            {
                patches[patchi][facei][0] = i   + Nx*j     + Nx*Ny*Ncellz;
                patches[patchi][facei][1] = (i+1)   + Nx*j + Nx*Ny*Ncellz;
                patches[patchi][facei][2] = (i+1) + Nx*(j+1) + Nx*Ny*Ncellz;
                patches[patchi][facei][3] = i + Nx*(j+1)     + Nx*Ny*Ncellz;

                facei++;
            }
        }

        faceListList patchLists(6);
        forAll(patches, patchi)
        {
          patchLists[patchi].setSize(patches[patchi].size());
          forAll(patches[patchi], facei)
          {
            patchLists[patchi][facei] = face(labelList(patches[patchi][facei]));
          }
        }

        wordList patchNames(6);
        patchNames[0] = "west";
        patchNames[1] = "east";
        patchNames[2] = "south";
        patchNames[3] = "north";
        patchNames[4] = "bottom";
        patchNames[5] = "top";

        dictionary defaultPatchDict;
        defaultPatchDict.add("type", "zeroGradient");
        PtrList<dictionary> patchDicts(6);
        for (int i = 0; i < 6; i++)
        {
            patchDicts.set(i, &defaultPatchDict);
        }
        Info << "Constructing mesh." << endl;
      

        // Construct polyMesh
        polyMesh mesh
        (
          IOobject
          (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
          ),
          clone(points),
          cellShapes,
          patchLists,
          patchNames,
          patchDicts,
          "defaultFaces",
          emptyPolyPatch::typeName
        );

        Info << "mesh cells: " << mesh.cells().size() << endl;

        if (!mesh.write())
        {
            FatalErrorInFunction
              << "Failed writing mesh"
              << exit(FatalError);
        }
    }

    Info << "Loading fvMesh" << endl;
    fvMesh fv_mesh
    (
        IOobject
        (
          polyMesh::defaultRegion,
          runTime.constant(),
          runTime
        )
    );
      
    Info << "Retrieving var value." << endl;
    // Get a variable to look how it behaves
    Foam::volVectorField Var
    (
      IOobject
      (
        "var",
        runTime.timeName(),
        fv_mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
      ),
      fv_mesh,
      dimVelocity
    );

    {
      Info << "tmp_var init" << endl;
      float* tmp_var_u = new float[N_Time * N_bottom_top * N_south_north * N_west_east_stag];
      float* tmp_var_v = new float[N_Time * N_bottom_top * N_south_north_stag * N_west_east];
      float* tmp_var_w = new float[N_Time * N_bottom_top_stag * N_south_north * N_west_east];
      netCDF::NcVar U_wrf = dataFile.getVar("U");
      netCDF::NcVar V_wrf = dataFile.getVar("V");
      netCDF::NcVar W_wrf = dataFile.getVar("W");
      Info << "load U_wrf" << endl;
      U_wrf.getVar(tmp_var_u);
      V_wrf.getVar(tmp_var_v);
      W_wrf.getVar(tmp_var_w);

      int cc = 0;
      for(int ibt = 0; ibt < Ncellz; ibt++)
      {
        for(int isn = 0; isn < Ncelly; isn++)
        {
          for(int iwe = 0; iwe < Ncellx; iwe++)
          {
            // The commented line produced a figure that is wrong
            // order in xy
            // This is probably due to the way data is stored 
            // in the netcdf file. 
            // cc = iwe + Ncellx*isn + Ncellx*Ncelly*ibt;
            cc = isn + Ncelly*iwe + Ncellx*Ncelly*ibt;
            Var.primitiveFieldRef()[cc][0] = 
                *(tmp_var_u + ibt * N_south_north *N_west_east_stag
                    + isn * N_west_east_stag + iwe);
            Var.primitiveFieldRef()[cc][1] = 
                *(tmp_var_v + ibt * N_south_north_stag *N_west_east
                    + isn * N_west_east + iwe);
            Var.primitiveFieldRef()[cc][2] = 
                *(tmp_var_w + ibt * N_south_north *N_west_east
                    + isn * N_west_east + iwe);
          }
        }
      }
    }

    runTime++;
    runTime.write();

    Info<< "End\n" << endl; 
    return 0;
}
