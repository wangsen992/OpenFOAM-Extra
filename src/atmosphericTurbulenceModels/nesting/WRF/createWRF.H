// Assuming foam mesh is already created
Info << "Creating WRF object..." << endl;
IOdictionary wrf_dict
(
  IOobject
  (
    "WRFDict",
    runTime.constant(),
    runTime,
    IOobject::MUST_READ
  )
);
Info << wrf_dict << endl;
WRF wrf
(
  wrf_dict["ncfile"], 
  wrf_dict["wrfRoot"], 
  "wrfCase",
  mesh,
  wrf_dict["proj4"],
  wrf_dict.lookupOrAddDefault("dt", 3600, false, false),
  wrf_dict.lookupOrAddDefault("T0", 290, false, false),
  wrf_dict.lookupOrAddDefault("P0", 1e5, false, false)
);

// When in parallel, terraform leads to mismatch at the proc boundaries
// by a small amount as the same point is translated by 
// slightly different amount
wrf.terraform_to_wrf(mesh);
wrf.updateVars(0);


// Load variables (besides U) to ptrTable
List<dictionary> variableList
(
  wrf_dict["variables"]
);
Info << "loaded variableList" << endl;
for(dictionary vd : variableList)
{
  Info << vd.lookupOrDefault<word>("name", "null") << endl;
}
