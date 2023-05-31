{
    autoPtr<phaseSystem::specieTransferTable>
        specieTransferPtr(fluid.specieTransfer());

    phaseSystem::specieTransferTable&
        specieTransfer(specieTransferPtr());

    fluid.correctReactions();

    forAll(fluid.multiComponentPhases(), multiComponentPhasei)
    {
        phaseModel& phase = fluid.multiComponentPhases()[multiComponentPhasei];

        UPtrList<volScalarField>& Y = phase.YActiveRef();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        forAll(Y, i)
        {
            Info << "[YEqn] Yi name : " << Y[i].name() 
                 << " min: " << min(Y[i].primitiveField())
                 << " average: " << average(Y[i].primitiveField())
                 << " max: " << max(Y[i].primitiveField()) << endl;
            fvScalarMatrix YiEqn
            (
                phase.YiEqn(Y[i])
             ==
               *specieTransfer[Y[i].name()]
              + fvModels.source(alpha, rho, Y[i])
            );

            YiEqn.relax();
            YiEqn.solve("Yi");
            Info << "[YEqn] Yi name : " << Y[i].name() 
                 << " min: " << min(Y[i].primitiveField())
                 << " average: " << average(Y[i].primitiveField())
                 << " max: " << max(Y[i].primitiveField()) << endl;
        }
    }

    fluid.correctSpecies();
}
