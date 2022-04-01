

      Info  <<  "simple.consistent() = "
            << simple.consistent() 
            << endl;
      if (simple.consistent())
      {
          rAtU = 1.0/(1.0/rAU - UEqn.H1());
          phiHbyA +=
              fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
          HbyA -= (rAU - rAtU())*fvc::grad(p);
      }

      tUEqn.clear();

      // Update the pressure BCs to ensure flux consistency
      // The incompressible form of constrainPressure() used which basically
      // overloaded the compressible form with rho being one field
      
      constrainPressure(p, U, phiHbyA, rAtU());

      // Non-orthogonal pressure corrector loop
      while (simple.correctNonOrthogonal())
      {
          fvScalarMatrix pEqn
          (
              fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
          );

          pEqn.setReference
          (
              pressureReference.refCell(),
              pressureReference.refValue()
          );

          pEqn.solve();

          if (simple.finalNonOrthogonalIter())
          {
              phi = phiHbyA - pEqn.flux();
          }
      }

      // compute cumulativeContErr
      volScalarField contErr(fvc::div(phi));

      scalar sumLocalContErr = runTime.deltaTValue()*
          mag(contErr)().weightedAverage(mesh.V()).value();

      scalar globalContErr = runTime.deltaTValue()*
          contErr.weightedAverage(mesh.V()).value();
      cumulativeContErr += globalContErr;

      Info<< "time step continuity errors : sum local = " << sumLocalContErr
          << ", global = " << globalContErr
          << ", cumulative = " << cumulativeContErr
          << endl;

      // Explicitly relax pressure for momentum corrector
      p.relax();

      // Momentum corrector
      U = HbyA - rAtU()*fvc::grad(p);
      U.correctBoundaryConditions();
