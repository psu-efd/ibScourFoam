// Optional suspended load (constant/sedimentProperties, entries of suspendedSedimentFoam)
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
            const scalar nu0 =
                gAverage(laminarTransport.nu()().primitiveField());
            suspendedPtr.reset
            (
                new ibSuspendedLoad(mesh, sedimentProperties, ib, nu0)
            );
        }
    }
}
