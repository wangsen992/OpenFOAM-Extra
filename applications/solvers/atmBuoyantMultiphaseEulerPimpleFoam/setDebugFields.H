phiAirDebugc = fvc::reconstruct
               (
                  fluid.phases()[0].phi()
                  *
                  mesh.magSf()
               );
divphiAirDebug = fvc::div(fluid.phases()[0].phi());
alphaRhoPhiAirDebugc = fvc::reconstruct
               (
                  fluid.phases()[0].alphaRhoPhi()
                  *
                  mesh.magSf()
               );
divalphaRhoPhiAirDebug = fvc::div(fluid.phases()[0].alphaRhoPhi());
