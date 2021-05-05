{
    surfaceScalarField phiC
    (
        IOobject
        (
            "phiC",
            runTime.timeName(),
            mesh
        ),
        phi + (mesh.Sf() & fvc::interpolate(Vs))
    );

    fvScalarMatrix CEqn
    (
          fvm::ddt(C)
        + fvm::div(phiC, C)
        - fvm::laplacian(turbulence->nuEff(), C)
    );


    CEqn.solve();

    Info<< "Suspended sediment volume fraction = "
        << C.weightedAverage(mesh.V()).value()
        << "  Min(C) = " << min(C).value()
        << "  Max(C) = " << max(C).value()
        << endl;

    //Clipping of C in case it goes wild
    C.min(0.6);// due to the shape of particle
    C.max(0.0);
}
