/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 ibScourFoam developers (Penn State)
-------------------------------------------------------------------------------
License
    This file is part of ibScourFoam, built on OpenFOAM.

    ibScourFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ibScourFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ibScourFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "ibSuspendedLoad.H"
#include "immersedBoundary.H"
#include "ibWallModel.H"
#include "fvm.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ibSuspendedLoad::ibSuspendedLoad
(
    const fvMesh& mesh,
    const dictionary& dict,
    const immersedBoundary& ib,
    const scalar nu
)
:
    mesh_(mesh),
    ib_(ib),
    ScC_(dict.getOrDefault<scalar>("ScC", 1)),
    sigmaC_(dict.getOrDefault<scalar>("sigmaC", 1)),
    mode_
    (
        dict.subOrEmptyDict("ibBedExchange").getOrDefault<word>("mode", "HRN")
    ),
    exponential_
    (
        dict.subOrEmptyDict("ibBedExchange").getOrDefault<word>
        (
            "form", "linear"
        ) == "exponential"
    ),
    kappa_
    (
        dict.subOrEmptyDict("ibBedExchange").getOrDefault<scalar>("kappa", 0.41)
    ),
    g_(dict.getOrDefault<scalar>("g", 9.81)),
    wsMag_(0),
    ws_(Zero),
    refModel_(suspendedSediment::referenceConcentrationModel::New(dict)),
    C_
    (
        IOobject
        (
            "C",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    DC_
    (
        IOobject
        (
            "DC",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimViscosity, Zero),
        fvPatchFieldBase::calculatedType()
    ),
    E_(),
    D_(),
    Cstar_(),
    erodible_(),
    rippleHeight_(0),
    rippleLength_(0),
    Ls_(0),
    mixingLayers_(5),
    ripplePeriod_(0),
    d50_(0),
    s_(0),
    thetaWindowPtr_(),
    thetaPeakPtr_(),
    rippleWindowTime_(0),
    haveThetaPeak_(false)
{
    if
    (
        mode_ != "referenceLevel" && mode_ != "HRN" && mode_ != "logLaw"
     && mode_ != "ripple"
    )
    {
        FatalIOErrorInFunction(dict)
            << "Unknown ibBedExchange mode " << mode_
            << ": referenceLevel | HRN | logLaw | ripple" << exit(FatalIOError);
    }

    if
    (
        mode_ != "ripple"
     && refModel_->needsShear() && refModel_->skinFriction()
    )
    {
        FatalIOErrorInFunction(dict)
            << "shieldsNumber skinFriction is not available over an immersed"
            << " bed (no depth-mean velocity per bed face yet); use"
            << " shieldsNumber total" << exit(FatalIOError);
    }

    wsMag_ =
        suspendedSediment::settlingVelocityModel::New(dict)->ws(nu, g_);
    ws_ = vector(0, 0, -wsMag_);

    if (mode_ == "ripple")
    {
        readRipples(dict);
    }

    Info<< "    Suspended load: ws " << wsMag_ << " m/s, ScC " << ScC_
        << ", sigmaC " << sigmaC_ << ", bed exchange " << mode_
        << (mode_ == "HRN" ? (exponential_ ? " exponential" : " linear") : "")
        << ", delta " << refModel_->delta() << " m" << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ibSuspendedLoad::readRipples(const dictionary& dict)
{
    const dictionary& exDict = dict.subDict("ibBedExchange");
    const dictionary& rDict = exDict.subDict("ripples");

    d50_ = dict.get<scalar>("d50");
    s_ = dict.get<scalar>("s");
    ripplePeriod_ = exDict.get<scalar>("period");
    mixingLayers_ = exDict.getOrDefault<scalar>("mixingLayers", 5);

    const word model(rDict.get<word>("model"));
    if (model == "fixed")
    {
        rippleHeight_ = rDict.get<scalar>("height");
        rippleLength_ = rDict.get<scalar>("length");
    }
    else if (model == "nielsen1981")
    {
        // Regular-wave ripples (Nielsen 1981; Nielsen 1992 Eqs. 2.3.2/3):
        // eta/A = 0.275 - 0.022 psi^0.5, lambda/A = 2.2 - 0.345 psi^0.34,
        // A = Um T/(2 pi), psi = Um^2/((s - 1) g d)
        const scalar Um = rDict.get<scalar>("Um");
        const scalar A = Um*ripplePeriod_/constant::mathematical::twoPi;
        const scalar psi = sqr(Um)/((s_ - 1)*g_*d50_);
        if (psi >= 156)
        {
            FatalIOErrorInFunction(rDict)
                << "nielsen1981 ripples: mobility number " << psi
                << " >= 156 (sheet flow, no ripples)" << exit(FatalIOError);
        }
        rippleHeight_ = A*(0.275 - 0.022*sqrt(psi));
        rippleLength_ = A*(2.2 - 0.345*pow(psi, 0.34));
        Info<< "    Ripples (Nielsen 1981): A " << A << " m, psi " << psi
            << endl;
    }
    else
    {
        FatalIOErrorInFunction(rDict)
            << "Unknown ripple model " << model << ": fixed | nielsen1981"
            << exit(FatalIOError);
    }

    if
    (
        rippleHeight_ <= 0 || rippleLength_ <= 0
     || constant::mathematical::pi*rippleHeight_ >= rippleLength_
    )
    {
        FatalIOErrorInFunction(rDict)
            << "Ripple height " << rippleHeight_ << " and length "
            << rippleLength_ << " must be positive with pi eta < lambda"
            << exit(FatalIOError);
    }
    if (ripplePeriod_ <= 0)
    {
        FatalIOErrorInFunction(exDict)
            << "ripple mode: period must be positive" << exit(FatalIOError);
    }

    Ls_ = 1.4*rippleHeight_;

    for (const word& name : {word("thetaRippleWindow"), word("thetaRipplePeak")})
    {
        IOobject io
        (
            name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        );
        const bool found = io.typeHeaderOk<volScalarField>(true);
        autoPtr<volScalarField>& ptr =
            (name == "thetaRippleWindow" ? thetaWindowPtr_ : thetaPeakPtr_);
        ptr.reset
        (
            new volScalarField
            (
                io,
                mesh_,
                dimensionedScalar(dimless, Zero),
                fvPatchFieldBase::zeroGradientType()
            )
        );
        if (found && name == "thetaRipplePeak")
        {
            haveThetaPeak_ = true;
        }
    }
    // The window restarts with the run (the sum of a partial window is not
    // stored); the last completed peak is kept

    Info<< "    Ripple suspension (Nielsen 1992): eta " << rippleHeight_
        << " m, lambda " << rippleLength_ << " m, L_s " << Ls_
        << " m, eps_s " << wsMag_*Ls_ << " m2/s up to "
        << mixingLayers_*Ls_ << " m above the bed, crest factor "
        << 1/sqr(1 - constant::mathematical::pi*rippleHeight_/rippleLength_)
        << ", window " << ripplePeriod_ << " s"
        << (haveThetaPeak_ ? " (peak read back)" : "") << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ibSuspendedLoad::correctDiffusivity
(
    const volScalarField& nu,
    const volScalarField& nut
)
{
    DC_ = nu/ScC_ + nut/sigmaC_;

    if (mode_ == "ripple")
    {
        // Ripple-vortex mixing in the layer above the bed (Nielsen 1992):
        // eps_s = w_s L_s, the diffusivity of the exponential profile
        const scalar epsS = wsMag_*Ls_;
        const scalar hMix = mixingLayers_*Ls_;
        const scalarField& dist = ib_.distance().primitiveField();
        scalarField& DCi = DC_.primitiveFieldRef();
        forAll(DCi, celli)
        {
            if (dist[celli] > 0 && dist[celli] < hMix)
            {
                DCi[celli] = max(DCi[celli], epsS);
            }
        }
        DC_.correctBoundaryConditions();
    }
}


void Foam::ibSuspendedLoad::solve
(
    const surfaceScalarField& phi,
    const scalar nu
)
{
    solveImpl(phi, nu, nullptr);
}


void Foam::ibSuspendedLoad::solve
(
    const surfaceScalarField& phiWater,
    const scalar nu,
    const volScalarField& alpha
)
{
    solveImpl(phiWater, nu, &alpha);
}


void Foam::ibSuspendedLoad::solveImpl
(
    const surfaceScalarField& phi,
    const scalar nu,
    const volScalarField* alpha
)
{
    // Face water fraction (two-phase): settling and diffusion stop at the
    // free surface
    tmp<surfaceScalarField> talphaf;
    if (alpha)
    {
        talphaf = fvc::interpolate(*alpha);
    }
    // Pickup needs a wet bed (alpha >= 0.5); settling onto the bed needs only
    // water in the cell: a draining film behind a bore keeps depositing (with
    // deposition gated by the wet threshold too, the grains of a drying IB
    // cell could not leave and C = (alpha C)/alpha grew beyond 1; dam break
    // over a mobile bed, 2026-11-02)
    auto wet = [&](const label celli) -> bool
    {
        return !alpha || (*alpha)[celli] >= 0.5;
    };
    const scalar alphaDry = 1e-3;
    auto holdsWater = [&](const label celli) -> bool
    {
        return !alpha || (*alpha)[celli] >= alphaDry;
    };

    const labelList& ibc = ib_.ibCells();
    const labelList& status = ib_.status();
    const scalarField& Aw = ib_.ibWallArea();
    const ibWallModel& wm = ib_.wallModel();
    const scalarField& uStar = wm.uTau();
    const scalarField& yc = wm.yIB();
    const vectorField& nCell = ib_.normal().primitiveField();
    const scalarField& V = mesh_.V();

    // Drift flux: settling across faces between fluid cells only
    surfaceScalarField phiC("phiC", phi);
    {
        volScalarField solidMask
        (
            IOobject("solidMask", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar(dimless, Zero),
            zeroGradientFvPatchField<scalar>::typeName
        );
        forAll(status, celli)
        {
            if (status[celli] == immersedBoundary::SOLID)
            {
                solidMask[celli] = 1;
            }
        }
        solidMask.correctBoundaryConditions();

        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const surfaceVectorField& Sf = mesh_.Sf();
        scalarField& phiCi = phiC.primitiveFieldRef();

        forAll(own, facei)
        {
            if (solidMask[own[facei]] < 0.5 && solidMask[nei[facei]] < 0.5)
            {
                phiCi[facei] +=
                    (ws_ & Sf[facei])*(alpha ? talphaf()[facei] : scalar(1));
            }
        }

        // Coupled patches: settling flux where both sides are fluid; no
        // settling through the other patches (as suspendedSedimentFoam)
        forAll(phiC.boundaryField(), patchi)
        {
            const fvPatch& p = mesh_.boundary()[patchi];
            if (!p.coupled())
            {
                continue;
            }
            const scalarField maskP(solidMask.boundaryField()[patchi]);
            const scalarField maskN
            (
                solidMask.boundaryField()[patchi].patchNeighbourField()
            );
            fvsPatchScalarField& phiCp = phiC.boundaryFieldRef()[patchi];
            const vectorField& Sfp = Sf.boundaryField()[patchi];
            forAll(phiCp, i)
            {
                if (maskP[i] < 0.5 && maskN[i] < 0.5)
                {
                    phiCp[i] +=
                        (ws_ & Sfp[i])
                       *(alpha ? talphaf().boundaryField()[patchi][i] : scalar(1));
                }
            }
        }
    }

    // Bed exchange per IB cell
    const label nIb = ibc.size();
    E_.resize(nIb, Zero);
    D_.resize(nIb, Zero);
    const bool ripple = (mode_ == "ripple");
    if (ripple)
    {
        // Peak skin-friction Shields number of the wave period, per cell
        volScalarField& thW = *thetaWindowPtr_;
        volScalarField& thP = *thetaPeakPtr_;
        const scalar theta0 = (s_ - 1)*g_*d50_;
        forAll(ibc, i)
        {
            const label celli = ibc[i];
            thW[celli] = max(thW[celli], sqr(uStar[i])/theta0);
        }
        rippleWindowTime_ += mesh_.time().deltaTValue();
        if (rippleWindowTime_ >= ripplePeriod_ - SMALL)
        {
            thP.primitiveFieldRef() = thW.primitiveField();
            thW.primitiveFieldRef() = Zero;
            rippleWindowTime_ = 0;
            haveThetaPeak_ = true;
            Info<< "Suspended load: ripple window closed, peak theta' max "
                << gMax(thP.primitiveField()) << endl;
        }
        const volScalarField& th = (haveThetaPeak_ ? thP : thW);

        // C0 = 0.005 theta_r^3 at the crest-trough mean bed (Nielsen 1986)
        const scalar crest =
            1/sqr(1 - constant::mathematical::pi*rippleHeight_/rippleLength_);
        Cstar_.resize(nIb);
        forAll(ibc, i)
        {
            Cstar_[i] = 0.005*pow3(crest*th[ibc[i]]);
        }
    }
    else
    {
        Cstar_ = refModel_->Cb(uStar, nu, g_);
    }

    // No pickup from fixed (non-erodible) bed faces
    if (erodible_.size() == Cstar_.size())
    {
        Cstar_ *= erodible_;
    }

    const scalar delta = refModel_->delta();
    const scalarField& DCi = DC_.primitiveField();

    volScalarField::Internal Su
    (
        IOobject("CSu", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    );
    volScalarField::Internal Sp
    (
        IOobject("CSp", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    );

    scalarField wn(nIb, Zero);
    auto depositionFactor = [&](const label i) -> scalar
    {
        return ripple ? exp(yc[i]/Ls_) : scalar(1);
    };
    label nPeclet = 0;
    scalar maxPeclet = 0;

    forAll(ibc, i)
    {
        const label celli = ibc[i];
        wn[i] = -(ws_ & nCell[celli]);

        if (wn[i] <= 0 || Aw[i] <= 0)
        {
            // Overhang or no wall area: no exchange
            continue;
        }

        const scalar Gamma1 = max(DCi[celli], SMALL);
        scalar Eexp = 0;

        if (mode_ == "referenceLevel" || ripple)
        {
            Eexp = wn[i]*Cstar_[i];
        }
        else if (mode_ == "HRN")
        {
            const scalar P2 = wn[i]*delta/Gamma1;
            const scalar E = wn[i]*Cstar_[i]*(exponential_ ? exp(P2) : 1 + P2);
            const scalar Pe1 = wn[i]*yc[i]/Gamma1;
            maxPeclet = max(maxPeclet, Pe1);
            if (Pe1 >= 1)
            {
                ++nPeclet;
            }
            Eexp = E*(1 - Pe1);
        }
        else
        {
            const scalar Zp = sigmaC_*wn[i]/(kappa_*max(uStar[i], SMALL));
            const scalar ratio = (yc[i] > delta ? delta/yc[i] : 1);
            Eexp = wn[i]*Cstar_[i]*pow(ratio, Zp);
        }

        E_[i] = Eexp;
        if (!wet(celli))
        {
            // Dry or air-covered bed: no pickup; grains still settle out of
            // whatever water the cell holds
            E_[i] = 0;
        }
        else
        {
            Su[celli] = Eexp*Aw[i]/V[celli];
        }
        if (holdsWater(celli))
        {
            // Ripple mode: deposition of the bed value of the exponential
            // profile through the IB-cell value
            Sp[celli] = -wn[i]*depositionFactor(i)*Aw[i]/V[celli];
        }
    }

    reduce(nPeclet, sumOp<label>());
    reduce(maxPeclet, maxOp<scalar>());

    if (nPeclet > 0)
    {
        FatalErrorInFunction
            << "HRN bed exchange: first-cell Peclet number w_n y_c/Gamma_1"
            << " >= 1 in " << nPeclet << " IB cells (max " << maxPeclet
            << "); the two-point wall flux would drive C negative."
            << " Use ibBedExchange mode logLaw" << exit(FatalError);
    }

    // Face diffusivity, zero on the faces of solid cells (no diffusive
    // leak into the solid, where C is fixed to zero)
    surfaceScalarField DCf("DCf", fvc::interpolate(DC_));
    if (alpha)
    {
        DCf *= talphaf();
    }
    {
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        scalarField& DCfi = DCf.primitiveFieldRef();
        forAll(own, facei)
        {
            if
            (
                status[own[facei]] == immersedBoundary::SOLID
             || status[nei[facei]] == immersedBoundary::SOLID
            )
            {
                DCfi[facei] = 0;
            }
        }
        forAll(DCf.boundaryField(), patchi)
        {
            const fvPatch& p = mesh_.boundary()[patchi];
            fvsPatchScalarField& DCfp = DCf.boundaryFieldRef()[patchi];
            const labelUList& fc = p.faceCells();
            forAll(fc, i)
            {
                if (status[fc[i]] == immersedBoundary::SOLID)
                {
                    DCfp[i] = 0;
                }
            }
        }
    }

    // Two-phase: the conserved quantity is alpha C, transported by the
    // water flux; the small diagonal keeps the rows of air cells regular
    // (no flux reaches them, so C decays to zero there)
    tmp<fvScalarMatrix> tCEqn;
    if (alpha)
    {
        // The row of a cell keeps the margin alpha_old V/dt + the small
        // diagonal over its off-diagonals (the water flux satisfies the
        // alpha equation), so a MULES undershoot alpha_old < 0 made it lose
        // diagonal dominance and the solve diverged (wave basin, alpha_old
        // -4e-6, 2026-11-03). The negative part of alpha_old goes to the
        // diagonal as well (non-conservative by |alpha_old| C per step there)
        const scalar rDeltaT = 1/mesh_.time().deltaTValue();
        volScalarField::Internal diagGuard
        (
            IOobject("CdiagGuard", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar(dimless/dimTime, Zero)
        );
        {
            const scalarField& a0 = alpha->oldTime().primitiveField();
            forAll(diagGuard, celli)
            {
                diagGuard[celli] =
                    (1e-6 + max(-a0[celli], scalar(0)))*rDeltaT;
            }
        }
        tCEqn = tmp<fvScalarMatrix>::New
        (
            fvm::ddt(*alpha, C_)
          + fvm::div(phiC, C_)
          - fvm::laplacian(DCf, C_)
          + fvm::Sp(diagGuard, C_)
         ==
            Su
          + fvm::Sp(Sp, C_)
        );
    }
    else
    {
        tCEqn = tmp<fvScalarMatrix>::New
        (
            fvm::ddt(C_)
          + fvm::div(phiC, C_)
          - fvm::laplacian(DCf, C_)
         ==
            Su
          + fvm::Sp(Sp, C_)
        );
    }
    fvScalarMatrix& CEqn = tCEqn.ref();

    const labelList& solid = ib_.solidCells();
    if (solid.size())
    {
        CEqn.setValues(solid, scalarField(solid.size(), Zero));
    }

    const solverPerformance sp = CEqn.solve();

    // A diverged solve leaves an unbounded C that the bed exchange would
    // carry into the Exner step: stop instead of continuing
    if
    (
        sp.nIterations() > 0
     && sp.finalResidual() > max(sp.initialResidual(), scalar(1e-6))
    )
    {
        FatalErrorInFunction
            << "Concentration solve diverged: initial residual "
            << sp.initialResidual() << ", final residual "
            << sp.finalResidual() << " after " << sp.nIterations()
            << " iterations" << exit(FatalError);
    }

    // Deposition with the new concentration (any IB cell holding water)
    scalar sumE = 0;
    scalar sumD = 0;
    scalar airLoss = 0;
    const scalar dt = mesh_.time().deltaTValue();
    forAll(ibc, i)
    {
        const label celli = ibc[i];
        D_[i] =
        (
            wn[i] > 0 && Aw[i] > 0 && holdsWater(celli)
          ? wn[i]*depositionFactor(i)*C_[celli]
          : 0
        );
        if (alpha && !holdsWater(celli) && Aw[i] > 0)
        {
            // A drained IB cell: its remaining grains settle onto the bed
            // in this step instead of being discarded
            D_[i] = max((*alpha)[celli]*C_[celli], scalar(0))*V[celli]/(Aw[i]*dt);
            C_[celli] = 0;
        }
        sumE += E_[i]*Aw[i];
        sumD += D_[i]*Aw[i];
    }

    if (alpha)
    {
        // No sediment in air: the grains of a drained cell off the bed are
        // discarded (reported; a drained IB cell deposited them above)
        scalarField& Ci = C_.primitiveFieldRef();
        forAll(Ci, celli)
        {
            if ((*alpha)[celli] < alphaDry && Ci[celli] != 0)
            {
                airLoss += max((*alpha)[celli]*Ci[celli], scalar(0))*V[celli];
                Ci[celli] = 0;
            }
        }
        C_.correctBoundaryConditions();
    }
    reduce(airLoss, sumOp<scalar>());
    reduce(sumE, sumOp<scalar>());
    reduce(sumD, sumOp<scalar>());

    Info<< "Suspended load: C min/max " << gMin(C_.primitiveField()) << " / "
        << gMax(C_.primitiveField()) << ", volume "
        << (alpha ? fvc::domainIntegrate((*alpha)*C_).value()
                  : fvc::domainIntegrate(C_).value())
        << " m3, entrainment " << sumE
        << " deposition " << sumD << " m3/s, C_b* max " << gMax(Cstar_)
        << (alpha ? ", discarded in air " : "")
        << (alpha ? Foam::name(airLoss) : word())
        << (mode_ == "HRN" ? ", max Pe1 " : "")
        << (mode_ == "HRN" ? Foam::name(maxPeclet) : word())
        << endl;
}


// ************************************************************************* //
