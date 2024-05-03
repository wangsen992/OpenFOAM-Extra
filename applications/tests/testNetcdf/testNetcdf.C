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

#include "WrfTools.H"
// #include "WrfCaseInfo.H"
#include "WrfOfReader.H"

#define DATESTRLEN 19
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void  printTimes(char times[][DATESTRLEN])
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

    struct WrfCaseInfo wrfInfo;
    readWrfCaseInfo(&wrfInfo, dataFile);

    mapProj::LambertConverter lam
    (
      wrfInfo.STAND_LON,
      wrfInfo.TRUELAT1,
      wrfInfo.TRUELAT2,
      6378137,  // Earth radius
      1,        // hemi
      1,        // dx
      1,        // dy
      wrfInfo.CEN_LAT,  // lat1
      wrfInfo.CEN_LON,  // lon1
      0,        // x1
      0         // y1
    );

    float neiu_lat(41.980206), neiu_lon(-87.717119);
    float neiu_x, neiu_y;
    float navy_pier_lat(41.891898), navy_pier_lon(-87.598688);
    float navy_pier_x, navy_pier_y;
    lam.get_xy(41.980206, -87.717119, neiu_x, neiu_y);
    lam.get_xy(navy_pier_lat, navy_pier_lon, navy_pier_x, navy_pier_y);
    Info << "NEIU XY: " << neiu_x << "; " << neiu_y << endl;
    Info << "NAVY PIER XY: " << navy_pier_x << "; " << navy_pier_y << endl;

    // Load mesh
    Info << "Read wrf case info" << endl;

    // autoPtr<polyMesh> tmesh(meshFromNc(dataFile, runTime));
    
    // Info << "[APP] mesh nCells: " << tmesh->nCells() << endl;
    #include "createMesh.H"

    WrfOfReader wrfReader(mesh, dataFile);
    volVectorField U
    (
          IOobject
          (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
          ),
          mesh,
          dimVelocity
    );
    Info << "N_Time = " << wrfInfo.N_Time << endl;
    for (int i = 0; i < wrfInfo.N_Time; i++)
    {
      Info << "Loading at timestep " << i << endl;
      autoPtr<volVectorField> pU(wrfReader.load_U(i));
      U = pU();
      pU.clear();
      runTime++;
      runTime.write();
    }


    Info<< "End\n" << endl; 
    return 0;
}
