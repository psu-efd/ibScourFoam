    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    immersedBoundaryFvMesh mesh
    (
        IOobject
        (
            immersedBoundaryFvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

