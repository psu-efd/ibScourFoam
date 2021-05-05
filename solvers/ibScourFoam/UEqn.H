// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);


fvVectorMatrix UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
    // It is easier to modify gradP than fvc::grad(p) in pEqn.H
    fvVectorMatrix UEqn1
    (
        UEqn == -gradP
    );

    // Fix velocity in IB/dead cells
    // also apply Dirichlet BC on velocity in laminar condition
    // In turbulent condition, velocity will be modified 
    // in the last timestep wall function

    mesh.manipulateUEqn(UEqn1);

    solve(UEqn1);

    fvOptions.correct(U);
}


