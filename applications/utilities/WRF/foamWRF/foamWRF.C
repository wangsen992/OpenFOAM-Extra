/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    foamWRF

Description
    A utility application for creating and transforming WRF nc files into OpenFOAM cases.

    Functions include: 
      1. Create mesh
      2. Load variables directly from the netcdf file
      3. transform the WRF mesh using proj4 string
      4. vertically transform mesh using wrf mesh (in the future with surface files)
      5. prepare variables (transformation) for simulation (from wrf model to foam model)


\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "dimensionedScalarFwd.H"
#include "fvCFD.H"
// #include "WRF.H"
#include "nesting_utils.H"
#include "treeDataFace.H"
#include "indexedOctree.H"
#include "treeBoundBox.H"
#include "fluidAtmThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption("get_wrf_proj4", "", "as name suggest");

    #include "setRootCase.H"
    #include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    IOdictionary dict
    (
      IOobject
      (
        "WrfDict",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ
      )
    );

    std::string filename(dict.lookup<string>("ncfile"));
    netCDF::NcFile ncfile(filename, netCDF::NcFile::FileMode::read);

    // Create mesh if option createMesh is specified
    autoPtr<fvMesh> pmesh;
    if(dict.lookup<bool>("createMesh"))
    {
      Info << "Create Mesh" << endl;
      pmesh.set(fvmeshFromNc(ncfile, runTime).ptr());
    }
    else
    {
      pmesh.set
      (
        new fvMesh
        (
          IOobject
          (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
          )
        )
      );
    }
    fvMesh& mesh(pmesh());

    // Read variables and write as needed
    if (dict.lookup<bool>("read_vars"))
    {
      dictionary var_dict = dict.subDict("read_varsCoeffs");
      scalar dt = var_dict.lookup<scalar>("dt");
      label tstart = var_dict.lookup<label>("tstart");
      label tend = var_dict.lookup<label>("tend");

      for(label i = tstart; i <= tend; i++)
      {
        Info << "t = " << i << endl;
        runTime.setTime(dt * i,i);
        if(var_dict.lookup<bool>("U"))
        {
          volVectorField U
          (
            load_U(pmesh(), ncfile, i)
          );
          U.write();
          U.clear();
        }

        if(var_dict.found("var3d"))
        {
          dictionary var_list = var_dict.subDict("var3d");
          for(auto varname : var_list.keys())
          {
            volScalarField var
            (
              load_var(pmesh(), ncfile, varname, var_list.subDict(varname).lookup<dimensionSet>("dim"), i)
            );
            if(var.name() == "P")
            {
              var.rename("PTURB");
            }
            var.write();
            var.clear();
          }
        }

        if(var_dict.found("var2d"))
        {
          List<word> var2d_list = var_dict.lookup<List<word>>("var2d");
          for(auto varname : var2d_list)
          {
            volScalarField var2d
            (
              load_2dvar(pmesh(), ncfile, varname, {0,0,0,1,0}, i)
            );
            var2d.write();
            var2d.clear();
          }
        }
      }
    }

    if(args.optionFound("get_wrf_proj4"))
    {
      WRF_PROJ_PARAMS params;
      getProjAtts(ncfile, params);
      std::shared_ptr<OGRSpatialReference> pcrs_wrf = getCRS(params);
      char* proj4[40];
      pcrs_wrf->exportToProj4(proj4);
      Info << *proj4 << endl;
    }

    // Trasnform Proj4
    if(dict.lookup<bool>("transform_proj4"))
    {
      dictionary proj4_dict(dict.subDict("transform_proj4Coeffs"));

      WRF_PROJ_PARAMS params;
      getProjAtts(ncfile, params);
      std::shared_ptr<OGRSpatialReference> pcrs_wrf = getCRS(params);
      std::shared_ptr<OGRSpatialReference> pcrs_foam
      (
        new OGRSpatialReference
      );

      string foam_proj4(proj4_dict.lookup<string>("proj4"));
      pcrs_foam->SetFromUserInput(foam_proj4.c_str());
      OGRCoordinateTransformation* ptransformer_;
      OGRCoordinateTransformation* pitransformer_;

      ptransformer_ = OGRCreateCoordinateTransformation(pcrs_wrf.get(), pcrs_foam.get());

      Info << "transform created" << endl;
      auto points(pmesh->points());
      std::for_each
      (
        points.begin(),
        points.end(),
        [&](point pt)
        {
          double x(pt.x()), y(pt.y()), z(pt.z());
          ptransformer_->Transform(1, &x, &y);
          return point{x,y,z};
        }
      );
      pmesh->movePoints(points);
      pmesh->setInstance(runTime.constant());
      pmesh->write();
    }

    if(dict.lookupOrDefault("transform_foam_mesh", false))
    {
      dictionary mesh_dict(dict.subDict("transform_foam_meshCoeffs"));
      Time foamTime
      (
        Time::controlDictName, 
        mesh_dict.lookup<string>("foam_case_root"),
        mesh_dict.lookup<string>("foam_case_name")
      );

      fvMesh foamMesh
      (
        IOobject
        (
          fvMesh::defaultRegion,
          foamTime.timeName(),
          foamTime,
          IOobject::MUST_READ
        )
      );

      if(mesh_dict.lookup<word>("method")=="WRF")
      {
        // Terraforming mesh
        Info << "Terraforming to WRF" << endl;
        fvMesh& wrfMesh(pmesh());
        const polyPatch& wrfGroundPatch(wrfMesh.boundaryMesh()["bottom"]);

        indexedOctree<treeDataFace> wrfTree
        (
          treeDataFace(false, wrfGroundPatch),
          treeBoundBox(boundBox(wrfGroundPatch.localPoints())),
          10,
          10,
          3
        );

        const polyPatch& foamGroundPatch(foamMesh.boundaryMesh()["bottom"]);
        indexedOctree<treeDataFace> foamTree
        (
          treeDataFace(false, foamGroundPatch),
          treeBoundBox(boundBox(foamGroundPatch.points())),
          10,
          10,
          3
        );
        // Interp the ground points from openfoam to wrf bottom patch
        Foam::vector ll{0,0,40000};
        Foam::vector zvec{0,0,1};
        auto findVec = [&](const point& pt)
        {
          auto wrfHit =  wrfTree.findLine(pt-ll, pt+ll);
          if (!wrfHit.hit())
          {
            Info << "WRF bounds " << wrfTree.bb() << endl;
            Info << "WRF not hit " << pt << endl;
          }
          point wrfPt = wrfHit.hitPoint(); 
          auto foamHit = foamTree.findLine(pt-ll, pt+ll);
          point foamPt;
          if (!foamHit.hit())
          {
            foamPt = foamTree.findNearest(pt, 1e9).hitPoint();
          }
          else
          {
            foamPt = foamHit.hitPoint(); 
          }

          return zvec * wrfPt.z();
          // compress the high altitude regions
        };

        {
          pointField foamPts = foamMesh.points();
          vectorField vec(foamPts.size());
          std::transform
          (
            foamPts.cbegin(), 
            foamPts.cend(), 
            vec.begin(),
            findVec
          );

          // assuming foam mesh is always lower than wrf mesh
          scalar zmax = gMax(foamPts.component(2));
          scalar zmin = gMin(foamPts.component(2));
          scalar vec_zmax(gMax(vec.component(2)));
          scalar vec_zmin(gMin(vec.component(2)));
          Info << vec_zmax - vec_zmin << endl;

          vec = Foam::vector{0,0,vec_zmin} 
              +(
                  (zmax - foamPts.component(2))/(zmax-zmin)
                 *(vec - Foam::vector{0,0,vec_zmin})
               );

          foamMesh.movePoints(foamPts + vec); // Some points are
                                                            // outside the wrf domain
          foamMesh.setInstance(foamMesh.time().constant());
          foamMesh.write();
          foamMesh.moving(false); // set to false so solver doesn't require V0
          Info << "Terraforming to WRF complete" << endl;
        }
      }
    }

    if(dict.lookup<bool>("prepare_variables"))
    {
      runTime.setTime(runTime.times()[0].value(), 0);
      dictionary prep_dict(dict.subDict("prepare_variablesCoeff"));
      dimensionedScalar T0 (prep_dict.lookup<dimensionedScalar>("T0"));
      dimensionedScalar P0 (prep_dict.lookup<dimensionedScalar>("P0"));
      volVectorField Uair
      (
        IOobject ( "U.air", runTime.timeName(), runTime, IOobject::NO_READ, IOobject::AUTO_WRITE),
        pmesh(),
        dimTemperature,
        "zeroGradient"
      );

      volScalarField Tair
      (
        IOobject ( "T.air", runTime.timeName(), runTime, IOobject::NO_READ, IOobject::AUTO_WRITE),
        pmesh(),
        dimensionedScalar(dimTemperature, 290),
        "zeroGradient"
      );
      Tair.write();
      volScalarField qv
      (
        IOobject ( "H2O.air", runTime.timeName(), runTime, IOobject::NO_READ, IOobject::AUTO_WRITE),
        pmesh(),
        dimless,
        "zeroGradient"
      );
      qv.write();
      volScalarField ydefault
      (
        IOobject ( "Ydefault.air", runTime.timeName(), runTime, IOobject::NO_READ, IOobject::AUTO_WRITE),
        pmesh(),
        dimensionedScalar(dimless, 1),
        "zeroGradient"
      );
      ydefault.write();
      volScalarField p
      (
        IOobject ( "p", runTime.timeName(), runTime, IOobject::NO_READ, IOobject::AUTO_WRITE),
        pmesh(),
        dimensionedScalar(dimPressure, 1e5),
        "zeroGradient"
      );
      p.write();

      volScalarField p_rgh
      (
        IOobject ( "p_rgh", runTime.timeName(), runTime, IOobject::NO_READ, IOobject::AUTO_WRITE),
        pmesh(),
        dimPressure,
        "zeroGradient"
      );
      volScalarField rho
      (
        IOobject ( "thermo:rho.air", runTime.timeName(), runTime, IOobject::NO_READ, IOobject::AUTO_WRITE),
        pmesh(),
        dimensionedScalar(dimDensity, 1),
        "zeroGradient"
      );

      #include "readpRef.H"
      #include "readGravitationalAcceleration.H"
      volVectorField ge
      (
        IOobject
        (
          "gEarth",
          runTime.constant(),
          runTime,
          IOobject::NO_READ
        ),
        pmesh(),
        dimensionedVector(dimAcceleration, vector(0,0,-9.81))
      );
      scalar re(6371000);
      scalarField h(mesh.C().component(2));
      ge.primitiveFieldRef() *=  sqr(re / (re + h));
      Info << max(ge) << endl;
      Info << min(ge) << endl;
      // #include "readhRef.H"
      // #include "gh.H"
      volScalarField gh(ge & mesh.C());

      // thermophysical models are used to compute thermodynamic properties
      Info << "Instantiating thermo model" << endl;
      autoPtr<fluidAtmThermo> pthermo
      (
        fluidAtmThermo::New(pmesh(), "air")
      );
      fluidAtmThermo& thermo(pthermo());

      for(int i = 0; i < runTime.times().size(); i++)
      {
        // Set the current time to load
        runTime.setTime(runTime.times()[i].value(), i);
        Info << "time = " << runTime.timeName() << endl;

        // Get the variables needed for this time
        volVectorField U ( IOobject ( "U", runTime.timeName(), runTime, IOobject::MUST_READ), pmesh());
        volScalarField T ( IOobject ( "T", runTime.timeName(), runTime, IOobject::MUST_READ), pmesh());
        volScalarField PTURB ( IOobject ( "PTURB", runTime.timeName(), runTime, IOobject::MUST_READ), pmesh());
        volScalarField PB ( IOobject ( "PB", runTime.timeName(), runTime, IOobject::MUST_READ), pmesh());
        volScalarField QV ( IOobject ( "QVAPOR", runTime.timeName(), runTime, IOobject::MUST_READ), pmesh());

        Info << "Updating T.air " << endl;
        Tair.primitiveFieldRef() = ((T + T0) * pow((PTURB + PB)/P0, 0.286))->primitiveField(); Tair.correctBoundaryConditions();
        Info << "Load qv" << endl;
        Uair.primitiveFieldRef() = U.primitiveField(); Uair.correctBoundaryConditions();
        qv.primitiveFieldRef() =  QV.primitiveField(); qv.correctBoundaryConditions();
        p.primitiveFieldRef() = (PTURB+PB)->primitiveField(); p.correctBoundaryConditions();
        thermo.he() = thermo.he(p, Tair);
        thermo.p() = p; 
        thermo.T() = Tair;
        thermo.correct();
        rho = thermo.rho();
        p_rgh = p - rho * gh - pRef;

        Uair.write();
        Tair.write();
        qv.write();
        p.write();
        rho.write();
        p_rgh.write();

        
        // volScalarFieldPtrTable_["thermo:rho.air"]() = volScalarFieldPtrTable_["p"]() / (volScalarFieldPtrTable_["T.air"]() * dimensionedScalar(dimEnergy/(dimMass*dimTemperature), 287.05));
        // volScalarFieldPtrTable_["thermo:rho.air"]().correctBoundaryConditions();
      }

    }

      Info << runTime.times()[0] << endl;
      Info << runTime.times()[1] << endl;
    // Conclude program with OpenFOAM boiler code
    Info << "Mesh points: " << pmesh->points().size() << endl;

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
