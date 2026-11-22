// Optional suspended load (constant/sedimentProperties, entries of suspendedSedimentFoam)
// Viscosity of the water (phase 1) for the wall law and the reference
// concentration of the suspended load
const scalar nuWater =
    gAverage(mixture.nuModel1().nu()().primitiveField());

// Eddy viscosity for the sediment diffusivity (zero when laminar)
volScalarField nutField
(
    IOobject("nutSediment", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
    mesh,
    dimensionedScalar(dimViscosity, Zero)
);
if (nutPtr)
{
    nutField = *nutPtr;
}

autoPtr<ibSuspendedLoad> suspendedPtr;
{
    IOobject sedIO
    (
        "sedimentProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (sedIO.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary sedimentProperties(sedIO);

        if (sedimentProperties.getOrDefault<bool>("solveSediment", true))
        {
            Info<< "Reading sedimentProperties: suspended load on" << endl;
            suspendedPtr.reset
            (
                new ibSuspendedLoad(mesh, sedimentProperties, ib, nuWater)
            );
        }
    }
}
