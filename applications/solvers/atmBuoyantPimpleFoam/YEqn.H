UPtrList<volScalarField>& Y = thermo.composition().Y();
forAll(Y, i)
{
    fvScalarMatrix YiEqn
    (
        fvm::ddt(rho, Y[i])
      + fvm::div(phi, Y[i])
      + thermophysicalTransport->divj(Y[i])
    ==
      fvModels.source(rho, Y[i])
    );

    YiEqn.relax();
    YiEqn.solve("Yi."+Y[i].name());
}

thermo.composition().normalise();
