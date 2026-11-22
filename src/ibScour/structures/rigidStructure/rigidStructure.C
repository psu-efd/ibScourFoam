/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Xiaofeng Liu, Penn State
-------------------------------------------------------------------------------
License
    This file is part of ibScourFoam, distributed under the GNU General
    Public License v3 (see the LICENSE file of OpenFOAM).

\*---------------------------------------------------------------------------*/

#include "rigidStructure.H"
#include "Time.H"
#include "IOdictionary.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::rigidStructure::rigidStructure
(
    const word& name,
    const dictionary& dict,
    const vector& defaultFree
)
:
    name_(name),
    density_(dict.getOrDefault<scalar>("density", 0)),
    mass_(dict.getOrDefault<scalar>("mass", -1)),
    volume_(dict.getOrDefault<scalar>("volume", -1)),
    freeTranslations_
    (
        dict.getOrDefault<vector>
        (
            "freeTranslations",
            dict.getOrDefault<vector>("directions", defaultFree)
        )
    ),
    freeRotations_(dict.getOrDefault<vector>("freeRotations", Zero)),
    refPoint_(dict.getOrDefault<point>("referencePoint", Zero)),
    refPointGiven_(dict.found("referencePoint")),
    centreOfMass_(dict.getOrDefault<point>("centreOfMass", Zero)),
    centreOfMassGiven_(dict.found("centreOfMass")),
    externalForce_
    (
        Function1<vector>::NewIfPresent("externalForce", dict)
    ),
    externalMoment_
    (
        Function1<vector>::NewIfPresent("externalMoment", dict)
    ),
    externalForcePoint_
    (
        dict.getOrDefault<point>("externalForcePoint", Zero)
    ),
    externalForcePointGiven_(dict.found("externalForcePoint")),
    appliedForce_(Zero),
    appliedMoment_(Zero),
    displacementFunction_
    (
        Function1<vector>::NewIfPresent("displacement", dict)
    ),
    rotationFunction_
    (
        Function1<vector>::NewIfPresent("rotation", dict)
    ),
    d_(Zero),
    theta_(Zero),
    kT_(vector::one*dict.getOrDefault<scalar>("stiffness", 1)),
    kR_(vector::one*dict.getOrDefault<scalar>("rotationalStiffness", 1)),
    K_(6, Zero),
    kT0_(kT_),
    kR0_(kR_),
    stiffnessResidualFloor_
    (
        dict.getOrDefault<scalar>("stiffnessResidualFloor", -1)
    ),
    relaxation_(dict.getOrDefault<scalar>("relaxation", 1)),
    gate_(dict.getOrDefault<scalar>("gate", 1e-3)),
    RF_(Zero),
    RM_(Zero),
    dF_(Zero),
    dM_(Zero),
    nUpdates_(0),
    motionSolver_(dict.getOrDefault<word>("motionSolver", "quasiStatic")),
    update_(dict.getOrDefault<word>("update", "inner")),
    inertia_(dict.getOrDefault<vector>("inertia", vector::one)),
    addedMass_(dict.getOrDefault<vector>("addedMass", Zero)),
    addedInertia_(dict.getOrDefault<vector>("addedInertia", Zero)),
    damping_(dict.getOrDefault<vector>("damping", Zero)),
    rotationalDamping_(dict.getOrDefault<vector>("rotationalDamping", Zero)),
    newmarkBeta_(dict.getOrDefault<scalar>("newmarkBeta", 0.25)),
    newmarkGamma_(dict.getOrDefault<scalar>("newmarkGamma", 0.5)),
    dynamicStartTime_(dict.getOrDefault<scalar>("dynamicStartTime", -GREAT)),
    dynamicActive_(false),
    nStaticSteps_(0),
    v_(Zero),
    omega_(Zero),
    a_(Zero),
    alpha_(Zero),
    d0_(Zero),
    theta0_(Zero),
    v0_(Zero),
    omega0_(Zero),
    a0_(Zero),
    alpha0_(Zero),
    dt_(0),
    dynamicStarted_(false)
{
    if (motionSolver_ != "quasiStatic" && motionSolver_ != "dynamic")
    {
        FatalIOErrorInFunction(dict)
            << "structure " << name_ << ": motionSolver " << motionSolver_
            << " is not quasiStatic or dynamic" << exit(FatalIOError);
    }
    if
    (
        dynamicSolver() && cmptMax(freeRotations_) > 0.5
     && !dict.found("inertia")
    )
    {
        FatalIOErrorInFunction(dict)
            << "structure " << name_
            << ": the dynamic solver with free rotations needs 'inertia'"
            << exit(FatalIOError);
    }
    if (dynamicSolver() && (newmarkBeta_ <= 0 || newmarkGamma_ <= 0))
    {
        FatalIOErrorInFunction(dict)
            << "structure " << name_
            << ": newmarkBeta and newmarkGamma must be positive"
            << exit(FatalIOError);
    }
    if (free())
    {
        if (density_ <= 0 && mass_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "structure " << name_
                << ": a free structure needs a positive density or a mass"
                << exit(FatalIOError);
        }
        if (cmptMin(kT_) <= 0 || cmptMin(kR_) <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "structure " << name_
                << ": stiffness and rotationalStiffness must be positive"
                << exit(FatalIOError);
        }
    }
    resetStiffness();

    if (update_ != "inner" && update_ != "outer")
    {
        FatalIOErrorInFunction(dict)
            << "structure " << name_ << ": update must be inner or outer"
            << exit(FatalIOError);
    }
    if (update_ == "outer" && dynamicSolver())
    {
        FatalIOErrorInFunction(dict)
            << "structure " << name_
            << ": the dynamic solver needs 'update inner'"
            << exit(FatalIOError);
    }
    if (stiffnessResidualFloor_ < 0)
    {
        // Outer update: the reaction is converged, and only steps made by
        // a residual above the soil accuracy carry information
        stiffnessResidualFloor_ = 1e-2;
    }
    if (mass_ >= 0 && volume_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "structure " << name_ << ": 'mass' needs 'volume'"
            << exit(FatalIOError);
    }
    if (volume_ >= 0)
    {
        Info<< "    Structure " << name_ << ": total volume " << volume_
            << " m3 (weight and buoyancy; the soil reports the embedded"
            << " volume)" << endl;
    }

    Info<< "    Structure " << name_ << ": "
        << (free() ? (dynamicSolver() ? "dynamic (Newmark)" : "quasi-static")
                   : "prescribed")
        << ", free translations " << freeTranslations_
        << ", rotations " << freeRotations_
        << (mass_ >= 0 ? ", mass " + Foam::name(mass_) + " kg"
                       : ", density " + Foam::name(density_) + " kg/m3")
        << ", stiffness estimates " << kT_.z() << " N/m, " << kR_.x()
        << " N m/rad, relaxation " << relaxation_ << endl;
    if (dynamicSolver())
    {
        Info<< "    Structure " << name_ << ": inertia " << inertia_
            << " kg m2, added mass " << addedMass_ << " kg, added inertia "
            << addedInertia_ << " kg m2, damping " << damping_
            << " N s/m, " << rotationalDamping_ << " N m s/rad, Newmark beta "
            << newmarkBeta_ << " gamma " << newmarkGamma_;
        if (dynamicStartTime_ > -GREAT/2)
        {
            Info<< ", quasi-static placement until t = " << dynamicStartTime_
                << " s, then from rest in balance";
        }
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::rigidStructure::free() const
{
    return cmptMax(freeTranslations_) > 0.5 || cmptMax(freeRotations_) > 0.5;
}


Foam::scalar Foam::rigidStructure::mass(const scalar V) const
{
    if (mass_ >= 0)
    {
        return mass_;
    }
    return density_*(volume_ >= 0 ? volume_ : V);
}


void Foam::rigidStructure::setRefPoint(const point& x)
{
    if (!refPointGiven_)
    {
        refPoint_ = x;
        refPointGiven_ = true;
    }
}


Foam::point Foam::rigidStructure::centreOfMass(const point& xV) const
{
    return (centreOfMassGiven_ ? centreOfMass_ : xV);
}


Foam::tmp<Foam::vectorField> Foam::rigidStructure::displacementAt
(
    const pointField& x
) const
{
    auto tD = tmp<vectorField>::New(x.size());
    vectorField& D = tD.ref();
    forAll(x, i)
    {
        D[i] = displacementAt(x[i]);
    }
    return tD;
}


void Foam::rigidStructure::setState(const vector& d, const vector& theta)
{
    d_ = d;
    theta_ = theta;
}


void Foam::rigidStructure::setPrescribed(const scalar t)
{
    if (displacementFunction_)
    {
        d_ = displacementFunction_->value(t);
    }
    if (rotationFunction_)
    {
        theta_ = rotationFunction_->value(t);
    }
}


Foam::vector Foam::rigidStructure::externalForce(const scalar t) const
{
    return
        appliedForce_
      + (externalForce_ ? externalForce_->value(t) : vector(Zero));
}


Foam::vector Foam::rigidStructure::externalMoment(const scalar t) const
{
    vector M(externalMoment_ ? externalMoment_->value(t) : vector(Zero));
    if (externalForce_ && externalForcePointGiven_)
    {
        M += (externalForcePoint_ - refPoint_) ^ externalForce_->value(t);
    }
    return appliedMoment_ + M;
}


void Foam::rigidStructure::gravityLoads
(
    const scalar V,
    const point& xV,
    const scalar rhoW,
    const scalar g,
    vector& F,
    vector& M
) const
{
    const scalar Vb = (volume_ >= 0 ? volume_ : V);
    const vector buoyancy(0, 0, rhoW*g*Vb);
    const vector weight(0, 0, -mass(V)*g);
    F = buoyancy + weight;
    M = ((xV - refPoint_) ^ buoyancy) + ((centreOfMass(xV) - refPoint_) ^ weight);
}


void Foam::rigidStructure::resetStiffness()
{
    K_ = Zero;
    for (direction c = 0; c < vector::nComponents; ++c)
    {
        K_(c, c) = kT_[c];
        K_(3 + c, 3 + c) = kR_[c];
    }
}


void Foam::rigidStructure::updateStiffness
(
    const vector& RF,
    const vector& RM,
    const scalar forceScale,
    const scalar momentScale
)
{
    if (nUpdates_ > 0)
    {
        // Scales per degree of freedom: the load, and the step the load
        // gives over the estimated stiffness; the step of the last update
        // and the change of the residual it caused
        FixedList<scalar, 6> Rs, qs, K0, dq, dR;
        for (direction c = 0; c < vector::nComponents; ++c)
        {
            Rs[c] = max(forceScale, SMALL);
            Rs[3 + c] = max(momentScale, SMALL);
            K0[c] = kT0_[c];
            K0[3 + c] = kR0_[c];
            qs[c] = Rs[c]/K0[c];
            qs[3 + c] = Rs[3 + c]/K0[3 + c];
            dq[c] = dF_[c];
            dq[3 + c] = dM_[c];
            dR[c] = RF[c] - RF_[c];
            dR[3 + c] = RM[c] - RM_[c];
        }

        // The residual that made the last step, in units of the load: a
        // step made by a residual within the accuracy of the reaction
        // (the soil iterate is settled to the gate only) gives a noise
        // secant
        scalar dq2 = 0;
        scalar R2Prev = 0;
        for (label i = 0; i < 6; ++i)
        {
            if (freeDof(i))
            {
                dq2 += sqr(dq[i]/qs[i]);
                R2Prev += sqr((i < 3 ? RF_[i] : RM_[i - 3])/Rs[i]);
            }
        }

        if (!outerUpdate())
        {
            // Within the correctors the reaction of the settled iterate
            // lags the step: the component secant, biased stiff by the
            // lag, is the under-relaxation the joint iteration needs (the
            // block rocking and balance tests converge with it, while
            // the coupled estimate held at the true stiffness oscillates
            // at residual 1e-4); kept when not positive, within a factor
            // 100 per update
            for (label i = 0; i < 6; ++i)
            {
                if (!freeDof(i) || mag(dq[i]) <= VSMALL) continue;
                const scalar K = -dR[i]/dq[i];
                if (K > 0)
                {
                    K_(i, i) = min(max(K, 0.01*K_(i, i)), 100*K_(i, i));
                }
            }
        }
        else if (R2Prev > sqr(stiffnessResidualFloor_) && dq2 > VSMALL)
        {
            // Broyden: the scaled stiffness Kh_ij = K_ij qs_j/Rs_i
            // satisfies Kh dqh = -dRh for the last step with the least
            // change; an update that makes a diagonal non-positive (the
            // reaction moved with the step: a reaction not settled, or a
            // changed contact state) is discarded
            const scalarSquareMatrix KPrev(K_);
            FixedList<scalar, 6> r(Zero);
            for (label i = 0; i < 6; ++i)
            {
                if (!freeDof(i)) continue;
                r[i] = -dR[i]/Rs[i];
                for (label j = 0; j < 6; ++j)
                {
                    if (freeDof(j))
                    {
                        r[i] -= K_(i, j)*dq[j]/Rs[i];
                    }
                }
            }
            for (label i = 0; i < 6; ++i)
            {
                if (!freeDof(i)) continue;
                for (label j = 0; j < 6; ++j)
                {
                    if (freeDof(j))
                    {
                        K_(i, j) += r[i]*(dq[j]/qs[j])/dq2*Rs[i]/qs[j];
                    }
                }
            }
            for (label i = 0; i < 6; ++i)
            {
                if (freeDof(i) && !(K_(i, i) > 0))
                {
                    K_ = KPrev;
                    break;
                }
            }
        }

        // Bounds: the diagonal within 1e-2..1e2 of the estimates (a
        // secant far above the estimate freezes the degree of freedom),
        // the off-diagonals within the geometric mean of the diagonals
        for (label i = 0; i < 6; ++i)
        {
            if (!freeDof(i)) continue;
            K_(i, i) = min(max(K_(i, i), 1e-2*K0[i]), 1e2*K0[i]);
        }
        for (label i = 0; i < 6; ++i)
        {
            if (!freeDof(i)) continue;
            for (label j = 0; j < 6; ++j)
            {
                if (i == j || !freeDof(j)) continue;
                const scalar lim = 0.9*sqrt(K_(i, i)*K_(j, j));
                K_(i, j) = min(max(K_(i, j), -lim), lim);
            }
        }
        for (direction c = 0; c < vector::nComponents; ++c)
        {
            kT_[c] = K_(c, c);
            kR_[c] = K_(3 + c, 3 + c);
        }
    }
    RF_ = RF;
    RM_ = RM;
}


void Foam::rigidStructure::solveStep
(
    const scalarSquareMatrix& K,
    const vector& RF,
    const vector& RM,
    vector& dF,
    vector& dM
) const
{
    dF = Zero;
    dM = Zero;

    DynamicList<label> dofs(6);
    for (label i = 0; i < 6; ++i)
    {
        if (freeDof(i)) dofs.append(i);
    }
    const label n = dofs.size();
    if (n == 0) return;

    scalarSquareMatrix A(n, Zero);
    List<scalar> b(n);
    forAll(dofs, i)
    {
        b[i] = (dofs[i] < 3 ? RF[dofs[i]] : RM[dofs[i] - 3]);
        forAll(dofs, j)
        {
            A(i, j) = K(dofs[i], dofs[j]);
        }
    }

    List<scalar> x(b);
    {
        scalarSquareMatrix LU(A);
        LUsolve(LU, x);
    }

    // A step that does not decrease the residual form (R . dq <= 0: K not
    // positive on it) falls back to the diagonal
    scalar bx = 0;
    forAll(x, i)
    {
        bx += b[i]*x[i];
    }
    if (!(bx > 0) || !std::isfinite(bx))
    {
        forAll(dofs, i)
        {
            x[i] = b[i]/A(i, i);
        }
    }

    forAll(dofs, i)
    {
        if (dofs[i] < 3)
        {
            dF[dofs[i]] = x[i];
        }
        else
        {
            dM[dofs[i] - 3] = x[i];
        }
    }
}


void Foam::rigidStructure::newmarkRates()
{
    const scalar b = newmarkBeta_;
    const scalar g = newmarkGamma_;
    const scalar dt = dt_;
    a_ = (d_ - d0_ - dt*v0_)/(b*sqr(dt)) - (1/(2*b) - 1)*a0_;
    v_ = v0_ + dt*((1 - g)*a0_ + g*a_);
    alpha_ = (theta_ - theta0_ - dt*omega0_)/(b*sqr(dt)) - (1/(2*b) - 1)*alpha0_;
    omega_ = omega0_ + dt*((1 - g)*alpha0_ + g*alpha_);
}


Foam::scalar Foam::rigidStructure::correct
(
    const vector& F,
    const vector& M,
    const scalar m,
    const scalar forceScale,
    const scalar momentScale
)
{
    if (dynamic())
    {
        return correctDynamic(F, M, m, forceScale, momentScale);
    }
    return correctQuasiStatic(F, M, forceScale, momentScale);
}


Foam::scalar Foam::rigidStructure::correctQuasiStatic
(
    const vector& F,
    const vector& M,
    const scalar forceScale,
    const scalar momentScale
)
{
    const vector RF(cmptMultiply(freeTranslations_, F));
    const vector RM(cmptMultiply(freeRotations_, M));

    updateStiffness(RF, RM, forceScale, momentScale);

    solveStep(K_, RF, RM, dF_, dM_);
    dF_ *= relaxation_;
    dM_ *= relaxation_;
    d_ += dF_;
    theta_ += dM_;
    ++nUpdates_;

    return max
    (
        cmptMax(cmptMag(RF))/max(forceScale, SMALL),
        cmptMax(cmptMag(RM))/max(momentScale, SMALL)
    );
}


Foam::scalar Foam::rigidStructure::correctDynamic
(
    const vector& F,
    const vector& M,
    const scalar m,
    const scalar forceScale,
    const scalar momentScale
)
{
    if (dt_ <= 0)
    {
        FatalErrorInFunction
            << "structure " << name_
            << ": the dynamic solver needs newStep(dt) with dt > 0"
            << exit(FatalError);
    }

    const vector RF(cmptMultiply(freeTranslations_, F));
    const vector RM(cmptMultiply(freeRotations_, M));

    // Start from rest: the initial acceleration balances the loads
    if (!dynamicStarted_)
    {
        for (direction c = 0; c < vector::nComponents; ++c)
        {
            a0_[c] = RF[c]/max(m + addedMass_[c], VSMALL);
            alpha0_[c] = RM[c]/max(inertia_[c] + addedInertia_[c], VSMALL);
        }
        a_ = a0_;
        alpha_ = alpha0_;
        dynamicStarted_ = true;
    }

    updateStiffness(RF, RM, forceScale, momentScale);
    newmarkRates();

    // Dynamic residual R = F - m a - c v at the current state and the
    // linearised balance (K + m/(beta dt^2) + c gamma/(beta dt)) dq = R
    const scalar b = newmarkBeta_;
    const scalar g = newmarkGamma_;
    const scalar dt = dt_;
    vector RFd(Zero);
    vector RMd(Zero);
    scalarSquareMatrix Keff(K_);
    for (direction c = 0; c < vector::nComponents; ++c)
    {
        const scalar mc = m + addedMass_[c];
        const scalar Ic = inertia_[c] + addedInertia_[c];
        RFd[c] = RF[c] - freeTranslations_[c]*(mc*a_[c] + damping_[c]*v_[c]);
        RMd[c] =
            RM[c]
          - freeRotations_[c]*(Ic*alpha_[c] + rotationalDamping_[c]*omega_[c]);
        Keff(c, c) += mc/(b*sqr(dt)) + damping_[c]*g/(b*dt);
        Keff(3 + c, 3 + c) += Ic/(b*sqr(dt)) + rotationalDamping_[c]*g/(b*dt);
    }
    solveStep(Keff, RFd, RMd, dF_, dM_);
    dF_ *= relaxation_;
    dM_ *= relaxation_;
    d_ += dF_;
    theta_ += dM_;
    newmarkRates();
    ++nUpdates_;

    return max
    (
        cmptMax(cmptMag(RFd))/max(forceScale, SMALL),
        cmptMax(cmptMag(RMd))/max(momentScale, SMALL)
    );
}


void Foam::rigidStructure::newStep(const scalar dt, const scalar t)
{
    // The next update has no previous residual of the same step to form a
    // secant with; the stiffness is kept. With the outer update the
    // residuals of successive solves form the secant.
    if (!outerUpdate())
    {
        nUpdates_ = 0;
    }

    if (dynamicSolver() && !dynamicActive_)
    {
        if (t > dynamicStartTime_ + SMALL)
        {
            dynamicActive_ = true;
            if (nStaticSteps_ > 0)
            {
                // Placed quasi-statically: the dynamics start from rest in
                // balance with the soil, not from the release acceleration
                // (with gamma 1/2 and a step far longer than the damped
                // settling time, the release acceleration would alternate
                // in sign from step to step and unload the contact)
                v_ = Zero;
                omega_ = Zero;
                a_ = Zero;
                alpha_ = Zero;
                dynamicStarted_ = true;
                Info<< "    Structure " << name_
                    << ": dynamic motion starts from rest at t = " << t - dt
                    << " s" << endl;
            }
        }
        else
        {
            ++nStaticSteps_;
        }
    }

    if (dynamic())
    {
        d0_ = d_;
        theta0_ = theta_;
        v0_ = v_;
        omega0_ = omega_;
        a0_ = a_;
        alpha0_ = alpha_;
        dt_ = dt;
        if (dt > 0)
        {
            // Predictor: constant acceleration over the step
            d_ = d0_ + dt*v0_ + 0.5*sqr(dt)*a0_;
            theta_ = theta0_ + dt*omega0_ + 0.5*sqr(dt)*alpha0_;
            newmarkRates();
        }
    }
}


void Foam::rigidStructure::writeHistory
(
    const fileName& caseDir,
    const scalar t
) const
{
    if (!Pstream::master())
    {
        return;
    }
    const fileName dir(caseDir/"postProcessing"/"structures");
    mkDir(dir);
    const fileName file(dir/name_ + ".dat");
    const bool exists = isFile(file);
    OFstream os(file, IOstreamOption(), IOstreamOption::APPEND_ATE);
    if (!exists)
    {
        os  << "# t [s]  d_x d_y d_z [m]  theta_x theta_y theta_z [rad]"
            << "  applied F_x F_y F_z [N]  applied M_x M_y M_z [N m]"
            << "  secant K_x K_y K_z [N/m]  Kr_x Kr_y Kr_z [N m/rad]"
            << "  v_x v_y v_z [m/s]  omega_x omega_y omega_z [rad/s]" << nl;
    }
    os  << t << ' ' << d_.x() << ' ' << d_.y() << ' ' << d_.z() << ' '
        << theta_.x() << ' ' << theta_.y() << ' ' << theta_.z() << ' '
        << appliedForce_.x() << ' ' << appliedForce_.y() << ' '
        << appliedForce_.z() << ' ' << appliedMoment_.x() << ' '
        << appliedMoment_.y() << ' ' << appliedMoment_.z() << ' '
        << kT_.x() << ' ' << kT_.y() << ' ' << kT_.z() << ' '
        << kR_.x() << ' ' << kR_.y() << ' ' << kR_.z() << ' '
        << v_.x() << ' ' << v_.y() << ' ' << v_.z() << ' '
        << omega_.x() << ' ' << omega_.y() << ' ' << omega_.z() << nl;
}


void Foam::rigidStructure::writeState(const Time& runTime) const
{
    IOdictionary dict
    (
        IOobject
        (
            name_,
            runTime.timeName(),
            "uniform"/fileName("structures"),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );
    dict.add("displacement", d_);
    dict.add("rotation", theta_);
    dict.add("velocity", v_);
    dict.add("angularVelocity", omega_);
    dict.add("acceleration", a_);
    dict.add("angularAcceleration", alpha_);
    dict.add("secantTranslation", kT_);
    dict.add("secantRotation", kR_);
    List<scalar> Kl(36);
    for (label i = 0; i < 6; ++i)
    {
        for (label j = 0; j < 6; ++j)
        {
            Kl[6*i + j] = K_(i, j);
        }
    }
    dict.add("stiffnessMatrix", Kl);
    dict.add("dynamicActive", dynamicActive_);
    dict.add("dynamicStarted", dynamicStarted_);
    dict.add("nStaticSteps", nStaticSteps_);
    dict.regIOobject::write();
}


bool Foam::rigidStructure::readState(const Time& runTime)
{
    IOobject io
    (
        name_,
        runTime.timeName(),
        "uniform"/fileName("structures"),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (!io.typeHeaderOk<IOdictionary>(true))
    {
        return false;
    }

    const IOdictionary dict(io);
    d_ = dict.get<vector>("displacement");
    theta_ = dict.get<vector>("rotation");
    v_ = dict.getOrDefault<vector>("velocity", Zero);
    omega_ = dict.getOrDefault<vector>("angularVelocity", Zero);
    a_ = dict.getOrDefault<vector>("acceleration", Zero);
    alpha_ = dict.getOrDefault<vector>("angularAcceleration", Zero);
    kT_ = dict.getOrDefault<vector>("secantTranslation", kT_);
    kR_ = dict.getOrDefault<vector>("secantRotation", kR_);
    List<scalar> Kl;
    if (dict.readIfPresent("stiffnessMatrix", Kl) && Kl.size() == 36)
    {
        for (label i = 0; i < 6; ++i)
        {
            for (label j = 0; j < 6; ++j)
            {
                K_(i, j) = Kl[6*i + j];
            }
        }
    }
    else
    {
        // Secants of an earlier version: within a decade of the estimates
        for (direction c = 0; c < vector::nComponents; ++c)
        {
            kT_[c] = min(max(kT_[c], 0.1*kT0_[c]), 10*kT0_[c]);
            kR_[c] = min(max(kR_[c], 0.1*kR0_[c]), 10*kR0_[c]);
        }
        resetStiffness();
    }
    dynamicActive_ = dict.getOrDefault<bool>("dynamicActive", false);
    dynamicStarted_ = dict.getOrDefault<bool>("dynamicStarted", false);
    nStaticSteps_ = dict.getOrDefault<label>("nStaticSteps", 0);
    d0_ = d_;
    theta0_ = theta_;
    v0_ = v_;
    omega0_ = omega_;
    a0_ = a_;
    alpha0_ = alpha_;

    Info<< "    Structure " << name_ << ": state restored from time "
        << runTime.timeName() << " (d " << d_ << " m, theta " << theta_
        << " rad)" << endl;

    return true;
}


void Foam::rigidStructure::report() const
{
    Info<< "    Structure " << name_ << ": d " << d_ << " m, theta "
        << theta_ << " rad";
    if (free())
    {
        Info<< ", secant " << kT_ << " N/m, " << kR_ << " N m/rad";
    }
    if (dynamicSolver())
    {
        Info<< ", v " << v_ << " m/s, omega " << omega_ << " rad/s";
    }
    Info<< endl;
}


// ************************************************************************* //
