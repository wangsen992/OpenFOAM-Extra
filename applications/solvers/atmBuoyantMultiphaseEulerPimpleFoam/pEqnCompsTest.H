PtrList<fvScalarMatrix> pEqnComps(phases.size());

{
    PtrList<volScalarField> dmdts(fluid.dmdts());
    PtrList<volScalarField> d2mdtdps(fluid.d2mdtdps());

    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;
        volScalarField& rho = phase.thermoRef().rho();

        pEqnComps.set(phasei, new fvScalarMatrix(p_rgh, dimVolume/dimTime));
        fvScalarMatrix& pEqnComp = pEqnComps[phasei];

        // Density variation
        if (!phase.isochoric() || !phase.pure())
        {
            // Adding the temporal MASS RATE CHANGE OF WATER is important to account 
            // for pressure balance in the system. Accounting for air is not
            // important. 
            if (phase.name() == "water")
            {
            pEqnComp +=
                (
                  //  fvc::div(phase.alphaRhoPhi())
                  // - fvc::Sp(fvc::div(phase.alphaPhi()), rho)
                    fvc::ddt(alpha, rho) + fvc::div(phase.alphaRhoPhi())
                  - fvc::Sp(fvc::ddt(alpha) + fvc::div(phase.alphaPhi()), rho)
                )/rho;
            }
            else
            {
            pEqnComp +=
                (
                    fvc::div(phase.alphaRhoPhi())
                  // - fvc::Sp(fvc::div(phase.alphaPhi()), rho)
                  - fvc::Sp(fvc::ddt(alpha) + fvc::div(phase.alphaPhi()), rho)
                )/rho;
            }
        }

        // Mesh dilatation correction
        if (mesh.moving())
        {
            pEqnComp += fvc::div(mesh.phi())*alpha;
        }

        // Compressibility
        // if (!phase.incompressible())
        // {
        //     if (pimple.transonic())
        //     {
        //         const surfaceScalarField phid
        //         (
        //             IOobject::groupName("phid", phase.name()),
        //             fvc::interpolate(phase.thermo().psi())*phase.phi()
        //         );

        //         pEqnComp +=
        //             correction
        //             (
        //                 (alpha/rho)*
        //                 (
        //                     phase.thermo().psi()*fvm::ddt(p_rgh)
        //                   + fvm::div(phid, p_rgh)
        //                   - fvm::Sp(fvc::div(phid), p_rgh)
        //                 )
        //             );

        //         pEqnComps[phasei].relax();
        //     }
        //     else
        //     {
        //      // This term is important 
                pEqnComp +=
                    (alpha*phase.thermo().psi()/rho)
                   *correction(fvm::ddt(p_rgh));
        //     }
        // }

        // Option sources
        if (fvModels.addsSupToField(rho.name()))
        {
            pEqnComp -= (fvModels.source(alpha, rho) & rho)/rho;
        }

        // Mass transfer
        // [Note] This term is disabled such that pressure equation reference 
        // can be set properly. 
        // [Note] Scrape the previous note. This term must be included for 
        // allowing phase fraction transfer (phase change) 
        if (dmdts.set(phasei))
        {
            pEqnComp -= dmdts[phasei]/rho;
        }
        if (d2mdtdps.set(phasei))
        {
            pEqnComp -= correction(fvm::Sp(d2mdtdps[phasei]/rho, p_rgh));
        }
    }
}
