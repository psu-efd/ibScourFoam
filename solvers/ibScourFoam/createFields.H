Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field gradP\n" << endl;
volVectorField gradP
(
    IOobject
    (
        "gradP",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(p)
);

Info<< "Reading field U_desired\n" << endl;
volVectorField U_desired
(
    IOobject
    (
        "U_desired",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("U_desired",U.dimensions(),vector::zero),
    U.boundaryField().types()
);

U_desired=U;

volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//sediment settling velocity field
volVectorField Vs
(
    IOobject
    (
        "Vs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("0.0", U.dimensions(), vector(0.0, 0.0, 0.0))
);


Info<< "Reading field Gamma\n" << endl;
volScalarField gamma
(
    IOobject
    (
        "Gamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);



volVectorField shearVelocity 
(
    IOobject
    (
        "shearVelocity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("0.0", U.dimensions(), vector(0.0, 0.0, 0.0))

);

volScalarField yplusIB 
(
    IOobject
    (
        "yplusIB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);

volScalarField yplus 
(
    IOobject
    (
        "yplus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);




volVectorField wallShearStress
(
    IOobject
    (
        "wallShearStress",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("wallShearStress", dimensionSet(0,2,-2,0,0,0,0),
     vector::zero)
 );
 





#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
