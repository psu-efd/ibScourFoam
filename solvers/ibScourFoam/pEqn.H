
volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
);

MRF.makeRelative(phiHbyA);

// Kill flux between dead/ib cells
mesh.ibCorrectPhi(phiHbyA);


// Adjust phiHbyA on IB faces to ensure mass conservation
mesh.immersedBoundaryAdjustPhi(phiHbyA,U);

adjustPhi(phiHbyA, U, p);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAU, MRF);



// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
    );
    pEqn.setReference(pRefCell, pRefValue);
    
    // Theorically, pEqn has to be taken care of so as to 
    // ensure Neumman BC.
    // However, it causes problems if the pressure gradient is high near IB.
    // For now, nothing has been done to pressure term.


    pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
    
    if (piso.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "immersedBoundaryContinuityErrs.H"


gradP == fvc::grad(p);

// Make gradP zero in IB/dead cells
// It is located in immersedBoundaryFvMeshPressureCorrector.C
mesh.evaluateGradP(gradP);

U = HbyA - rAU*gradP;



U.correctBoundaryConditions();

fvOptions.correct(U);
