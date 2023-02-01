phiAirDebugc = fvc::reconstruct
               (
                  turbulence->phi()
                  *
                  mesh.magSf()
               );
divphiAirDebug = fvc::div(turbulence->phi());

Info << average(turbulence->alphaRhoPhi())<< endl;
alphaRhoPhiAirDebugc = fvc::reconstruct
               (
                  turbulence->alphaRhoPhi()
                  *
                  mesh.magSf()
               );
divalphaRhoPhiAirDebug = fvc::div(turbulence->alphaRhoPhi());
