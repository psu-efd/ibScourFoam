/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Xiaofeng Liu, Penn State University
-------------------------------------------------------------------------------
License
    This file is part of ibScourFoam, built on OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sedimentBedFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "turbulenceModel.H"
#include "IOdictionary.H"
#include "gravityMeshObject.H"
#include "IStringStream.H"
#include "I1Tables.H"

// * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * * //

const Foam::Enum<Foam::sedimentBedFluxFvPatchScalarField::modeType>
Foam::sedimentBedFluxFvPatchScalarField::modeTypeNames
({
    { modeType::referenceLevel, "referenceLevel" },
    { modeType::HRN, "HRN" },
    { modeType::LRN, "LRN" },
    { modeType::Robin, "Robin" },
    { modeType::nearestCell, "nearestCell" },
});


const Foam::Enum<Foam::sedimentBedFluxFvPatchScalarField::diffusivityType>
Foam::sedimentBedFluxFvPatchScalarField::diffusivityTypeNames
({
    { diffusivityType::patchValue, "patch" },
    { diffusivityType::firstCell, "firstCell" },
    { diffusivityType::logLaw, "logLaw" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * //

namespace Foam
{

//- Built-in I1(ks+, delta+) table (I1Tables.H)
static interpolation2DTable<scalar>* builtinI1Table(const word& name)
{
    const char* txt =
    (
        name == "SST" ? suspendedSediment::I1TableSST
      : name == "vanDriestCC" ? suspendedSediment::I1TableVanDriestCC
      : nullptr
    );

    if (!txt)
    {
        FatalErrorInFunction
            << "Unknown I1Model " << name
            << " (SST, vanDriestCC or file)" << exit(FatalError);
    }

    IStringStream is(txt);
    List<Tuple2<scalar, List<Tuple2<scalar, scalar>>>> data(is);

    return new interpolation2DTable<scalar>
    (
        data,
        bounds::normalBounding::CLAMP,
        fileName("I1Table" + name)
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::dictionary&
Foam::sedimentBedFluxFvPatchScalarField::sedimentDict() const
{
    const word dictName("sedimentProperties");

    if (!db().foundObject<IOdictionary>(dictName))
    {
        IOdictionary* dictPtr = new IOdictionary
        (
            IOobject
            (
                dictName,
                db().time().constant(),
                db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            )
        );
        dictPtr->store();
    }

    return db().lookupObject<IOdictionary>(dictName);
}


const Foam::suspendedSediment::referenceConcentrationModel&
Foam::sedimentBedFluxFvPatchScalarField::refModel() const
{
    if (!refModel_)
    {
        refModel_ =
            suspendedSediment::referenceConcentrationModel::New(sedimentDict());
    }

    return *refModel_;
}


Foam::tmp<Foam::scalarField>
Foam::sedimentBedFluxFvPatchScalarField::calcUStar() const
{
    const turbulenceModel& turb = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField>(UName_);
    const vectorField nf(patch().nf());
    vectorField dUdn(Up.snGrad());
    dUdn -= nf*(nf & dUdn);

    return sqrt(turb.nuEff(patch().index())*mag(dUdn));
}


Foam::tmp<Foam::scalarField>
Foam::sedimentBedFluxFvPatchScalarField::faceDiffusivity() const
{
    if (diffusivity_ == diffusivityType::patchValue)
    {
        return nullptr;
    }

    const fvPatchScalarField& DCp =
        patch().lookupPatchField<volScalarField>(DCName_);
    tmp<scalarField> tGamma(DCp.patchInternalField());

    if (diffusivity_ == diffusivityType::logLaw)
    {
        scalarField& Gamma = tGamma.ref();
        const scalar delta = refModel().delta();
        const scalar sigmaC = sedimentDict().get<scalar>("sigmaC");
        const scalarField uStar(calcUStar());
        const scalarField wn
        (
            db().lookupObject<uniformDimensionedVectorField>(wsName_).value()
          & patch().nf()
        );
        const scalarField& dc = patch().deltaCoeffs();

        forAll(Gamma, facei)
        {
            const scalar y1 = 1/dc[facei];
            if (y1 <= delta || uStar[facei] <= SMALL)
            {
                continue;       // first-cell value
            }
            const scalar L = log(y1/delta);
            const scalar x = sigmaC*wn[facei]/(kappa_*uStar[facei])*L;  // Z' ln(y1/delta)
            Gamma[facei] =
            (
                mag(x) > 1e-8
              ? wn[facei]*y1/(-expm1(-x))
              : kappa_*uStar[facei]*y1/(sigmaC*L)     // limit ws -> 0
            );
        }
    }

    return tGamma;
}


void Foam::sedimentBedFluxFvPatchScalarField::depthMean
(
    scalarField& V,
    scalarField& D
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const vectorField nf(patch().nf());
    const vectorField Cf(patch().Cf());
    const vectorField& faceCentres = mesh.faceCentres();

    // Full cell columns above the patch faces, to the boundary at the top
    // (static mesh: built once)
    if (fullColumnCells_.empty())
    {
        fullColumnCells_.setSize(size());
        fullColumnDy_.setSize(size());
        columnDepth_.setSize(size());

        forAll(*this, facei)
        {
            DynamicList<label> cells;
            DynamicList<scalar> dys;

            label facej = patch().start() + facei;
            label celli = patch().faceCells()[facei];
            scalar y0 = 0;

            while (true)
            {
                const label opp =
                    mesh.cells()[celli].opposingFaceLabel(facej, mesh.faces());

                if (opp < 0)
                {
                    FatalErrorInFunction
                        << "No opposite face above face " << facei
                        << " of patch " << patch().name()
                        << ": the skin-friction Shields number needs"
                        << " hexahedral or prismatic cell columns"
                        << " (or shieldsNumber total in sedimentProperties)"
                        << exit(FatalError);
                }

                const scalar y = -(nf[facei] & (faceCentres[opp] - Cf[facei]));

                cells.append(celli);
                dys.append(y - y0);
                y0 = y;

                if (!mesh.isInternalFace(opp))
                {
                    const polyPatch& pp =
                        mesh.boundaryMesh()[mesh.boundaryMesh().whichPatch(opp)];

                    if (pp.coupled())
                    {
                        FatalErrorInFunction
                            << "The cell column above face " << facei
                            << " of patch " << patch().name()
                            << " reaches the coupled patch " << pp.name()
                            << ": the skin-friction Shields number needs the"
                            << " whole column on one processor (decompose"
                            << " along the bed, e.g. simple in x) or"
                            << " shieldsNumber total in sedimentProperties"
                            << exit(FatalError);
                    }
                    break;
                }

                celli =
                (
                    mesh.faceOwner()[opp] == celli
                  ? mesh.faceNeighbour()[opp]
                  : mesh.faceOwner()[opp]
                );
                facej = opp;
            }

            fullColumnCells_[facei].transfer(cells);
            fullColumnDy_[facei].transfer(dys);
            columnDepth_[facei] = y0;
        }
    }

    // Height-weighted mean of the cell velocities, tangential part
    const volVectorField& U = db().lookupObject<volVectorField>(UName_);

    V.setSize(size());
    D = columnDepth_;

    forAll(V, facei)
    {
        vector Um(Zero);
        forAll(fullColumnCells_[facei], k)
        {
            Um += fullColumnDy_[facei][k]*U[fullColumnCells_[facei][k]];
        }
        Um /= max(columnDepth_[facei], VSMALL);
        Um -= nf[facei]*(nf[facei] & Um);
        V[facei] = mag(Um);
    }
}


void Foam::sedimentBedFluxFvPatchScalarField::checkFirstCellPeclet
(
    const scalarField& wn
) const
{
    const scalar sigmaC = sedimentDict().get<scalar>("sigmaC");

    scalar PeMax = 0;
    label nBad = 0;

    forAll(wn, facei)
    {
        if (uStar_[facei] > ROOTVSMALL)
        {
            const scalar Pe = sigmaC*wn[facei]/(kappa_*uStar_[facei]);
            PeMax = max(PeMax, Pe);
            if (Pe >= 1)
            {
                ++nBad;
            }
        }
    }

    reduce(PeMax, maxOp<scalar>());
    reduce(nBad, sumOp<label>());

    if (nBad)
    {
        FatalErrorInFunction
            << "mode HRN with diffusivity firstCell on patch "
            << patch().name() << ": the first-cell Peclet number"
            << " sigma_c w_s/(kappa u*) = " << PeMax << " >= 1 on " << nBad
            << " bed face(s)." << nl
            << "    With wall functions the wall-cell diffusivity is"
            << " Gamma_1 = kappa u* y1/sigma_c, so w_s y1/Gamma_1 equals this"
            << " Rouse number for every y1: the two-point wall flux then"
            << " drives C negative, and refining the mesh does not help."
            << nl
            << "    Remedy: set 'diffusivity logLaw;' in the sedimentBedFlux"
            << " entry of this patch (wall flux exact for the equilibrium"
            << " profile of a log layer)." << exit(FatalError);
    }
}


Foam::tmp<Foam::scalarField>
Foam::sedimentBedFluxFvPatchScalarField::wallIntegral
(
    const scalarField& wn,
    const scalarField& Gamma,
    const scalar delta
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const vectorField& Cc = mesh.C();
    const vectorField nf(patch().nf());
    const vectorField Cf(patch().Cf());

    // Cell columns above the patch faces (static mesh: built once)
    if (columnCells_.empty())
    {
        columnCells_.setSize(size());
        columnY_.setSize(size());

        forAll(*this, facei)
        {
            DynamicList<label> cells;
            DynamicList<scalar> ys;

            label facej = patch().start() + facei;
            label celli = patch().faceCells()[facei];

            while (true)
            {
                const scalar y = -(nf[facei] & (Cc[celli] - Cf[facei]));
                cells.append(celli);
                ys.append(y);

                if (y >= delta)
                {
                    break;
                }

                const label opp =
                    mesh.cells()[celli].opposingFaceLabel(facej, mesh.faces());

                if (opp < 0 || !mesh.isInternalFace(opp))
                {
                    FatalErrorInFunction
                        << "Cell column of face " << facei << " of patch "
                        << patch().name() << " ends below delta = " << delta
                        << " (needs hexahedral or prismatic cells)"
                        << exit(FatalError);
                }

                celli =
                (
                    mesh.faceOwner()[opp] == celli
                  ? mesh.faceNeighbour()[opp]
                  : mesh.faceOwner()[opp]
                );
                facej = opp;
            }

            columnCells_[facei].transfer(cells);
            columnY_[facei].transfer(ys);
        }
    }

    // Trapezoidal rule for w_n/eps from the wall (eps = Gamma_f) to delta,
    // through the cell centres; the last interval is cut at delta
    const volScalarField& DC = db().lookupObject<volScalarField>(DCName_);

    auto tI = tmp<scalarField>::New(size(), Zero);
    scalarField& I = tI.ref();

    forAll(I, facei)
    {
        scalar y0 = 0;
        scalar f0 = 1/Gamma[facei];

        forAll(columnCells_[facei], k)
        {
            const scalar y1 = columnY_[facei][k];
            const scalar f1 = 1/DC[columnCells_[facei][k]];

            if (y1 >= delta)
            {
                const scalar fd = f0 + (f1 - f0)*(delta - y0)/(y1 - y0);
                I[facei] += 0.5*(f0 + fd)*(delta - y0);
                break;
            }

            I[facei] += 0.5*(f0 + f1)*(y1 - y0);
            y0 = y1;
            f0 = f1;
        }

        I[facei] *= wn[facei];
    }

    return tI;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sedimentBedFluxFvPatchScalarField::sedimentBedFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    mode_(modeType::HRN),
    diffusivity_(diffusivityType::firstCell),
    kappa_(0.41),
    exponential_(false),
    useTable_(false),
    I1Model_("vanDriestCC"),
    ks_(0),
    beta_(0),
    DCName_("DC"),
    wsName_("ws"),
    UName_("U"),
    E_(p.size(), Zero),
    beta_f_(p.size(), Zero),
    Cstar_(p.size(), Zero),
    uStar_(p.size(), Zero),
    uStarSkin_(p.size(), Zero)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::sedimentBedFluxFvPatchScalarField::sedimentBedFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    mode_(modeTypeNames.get("mode", dict)),
    diffusivity_
    (
        // default: firstCell for HRN (reproduces Liu 2014), patch
        // for the other modes (LRN: the wall value nu/ScC is the eps_s at
        // the wall; referenceLevel and Robin: the face value)
        diffusivityTypeNames.getOrDefault
        (
            "diffusivity",
            dict,
            mode_ == modeType::HRN
          ? diffusivityType::firstCell
          : diffusivityType::patchValue
        )
    ),
    kappa_(dict.getOrDefault<scalar>("kappa", 0.41)),
    exponential_(dict.getOrDefault<word>("form", "linear") == "exponential"),
    useTable_(dict.getOrDefault<word>("integral", "computed") == "table"),
    I1Model_
    (
        dict.getOrDefault<word>("I1Model", dict.found("I1Table") ? "file" : "vanDriestCC")
    ),
    ks_(dict.getOrDefault<scalar>("ks", 0)),
    beta_(dict.getOrDefault<scalar>("beta", 0)),
    DCName_(dict.getOrDefault<word>("DC", "DC")),
    wsName_(dict.getOrDefault<word>("ws", "ws")),
    UName_(dict.getOrDefault<word>("U", "U")),
    E_(p.size(), Zero),
    beta_f_(p.size(), Zero),
    Cstar_(p.size(), Zero),
    uStar_(p.size(), Zero),
    uStarSkin_(p.size(), Zero)
{
    fvPatchFieldBase::readDict(dict);

    if (mode_ == modeType::LRN && useTable_)
    {
        if (I1Model_ == "file")
        {
            I1Table_.reset
            (
                new interpolation2DTable<scalar>(dict.subDict("I1Table"))
            );
        }
        else
        {
            I1Table_.reset(builtinI1Table(I1Model_));
        }
    }

    if (mode_ == modeType::Robin && !dict.found("beta"))
    {
        FatalIOErrorInFunction(dict)
            << "mode Robin needs the deposition velocity beta [m/s]"
            << exit(FatalIOError);
    }

    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;

    if (!this->readValueEntry(dict))
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }
}


Foam::sedimentBedFluxFvPatchScalarField::sedimentBedFluxFvPatchScalarField
(
    const sedimentBedFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    mode_(ptf.mode_),
    diffusivity_(ptf.diffusivity_),
    kappa_(ptf.kappa_),
    exponential_(ptf.exponential_),
    useTable_(ptf.useTable_),
    I1Model_(ptf.I1Model_),
    I1Table_
    (
        ptf.I1Table_ ? new interpolation2DTable<scalar>(*ptf.I1Table_) : nullptr
    ),
    ks_(ptf.ks_),
    beta_(ptf.beta_),
    DCName_(ptf.DCName_),
    wsName_(ptf.wsName_),
    UName_(ptf.UName_),
    E_(p.size(), Zero),
    beta_f_(p.size(), Zero),
    Cstar_(p.size(), Zero),
    uStar_(p.size(), Zero),
    uStarSkin_(p.size(), Zero)
{}


Foam::sedimentBedFluxFvPatchScalarField::sedimentBedFluxFvPatchScalarField
(
    const sedimentBedFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    mode_(ptf.mode_),
    diffusivity_(ptf.diffusivity_),
    kappa_(ptf.kappa_),
    exponential_(ptf.exponential_),
    useTable_(ptf.useTable_),
    I1Model_(ptf.I1Model_),
    I1Table_
    (
        ptf.I1Table_ ? new interpolation2DTable<scalar>(*ptf.I1Table_) : nullptr
    ),
    ks_(ptf.ks_),
    beta_(ptf.beta_),
    DCName_(ptf.DCName_),
    wsName_(ptf.wsName_),
    UName_(ptf.UName_),
    E_(ptf.E_),
    beta_f_(ptf.beta_f_),
    Cstar_(ptf.Cstar_),
    uStar_(ptf.uStar_),
    uStarSkin_(ptf.uStarSkin_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sedimentBedFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Diffusivity of the C equation: patch and first-cell values (with
    // diffusivity firstCell or logLaw the solver has already replaced the
    // patch value of DC by faceDiffusivity())
    const fvPatchScalarField& DCp =
        patch().lookupPatchField<volScalarField>(DCName_);
    const scalarField Gamma(DCp);
    const scalarField GammaP(DCp.patchInternalField());

    // Settling velocity through the patch (> 0: out of the fluid)
    const uniformDimensionedVectorField& ws =
        db().lookupObject<uniformDimensionedVectorField>(wsName_);
    const scalarField wn(ws.value() & patch().nf());

    const scalar g = mag(meshObjects::gravity::New(db().time()).value());

    const turbulenceModel& turb = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalar nu = gAverage(turb.nu(patchi)());

    const bool hrnFirstCell =
    (
        mode_ == modeType::HRN && diffusivity_ == diffusivityType::firstCell
    );

    if
    (
        refModel().needsShear()
     || mode_ == modeType::LRN
     || diffusivity_ == diffusivityType::logLaw
     || hrnFirstCell
    )
    {
        uStar_ = calcUStar();
    }

    if (hrnFirstCell)
    {
        checkFirstCellPeclet(wn);
    }

    // Shear velocity of the Shields number: skin friction (from the
    // depth-mean velocity and the depth above each face) or total
    if (refModel().needsShear() && refModel().skinFriction())
    {
        scalarField V, D;
        depthMean(V, D);
        uStarSkin_ = refModel().uStarSkin(uStar_, V, D, g);
    }
    else
    {
        uStarSkin_ = uStar_;
    }

    Cstar_ = refModel().Cb(uStarSkin_, nu, g);

    // No pickup from fixed (non-erodible) bed faces
    if (erodible_.size() == Cstar_.size())
    {
        Cstar_ *= erodible_;
    }

    if (debug)
    {
        // DebugSwitches { sedimentBedFlux 1; } in system/controlDict
        Info<< "sedimentBedFlux " << patch().name()
            << ": u* min/max " << gMin(uStar_) << " " << gMax(uStar_)
            << ", Shields u* min/max " << gMin(uStarSkin_) << " "
            << gMax(uStarSkin_)
            << ", C_b* min/max " << gMin(Cstar_) << " " << gMax(Cstar_)
            << endl;
    }
    const scalar delta = refModel().delta();

    beta_f_ = wn;

    switch (mode_)
    {
        case modeType::referenceLevel:
        {
            E_ = wn*Cstar_;
            break;
        }

        case modeType::HRN:
        {
            if (diffusivity_ == diffusivityType::logLaw)
            {
                // the face value stands for C at delta
                E_ = wn*Cstar_;
                break;
            }
            const scalarField epsS
            (
                diffusivity_ == diffusivityType::patchValue ? Gamma : GammaP
            );
            const scalarField P2(wn*delta/epsS);

            E_ = wn*Cstar_*(exponential_ ? exp(P2) : 1 + P2);

            break;
        }

        case modeType::LRN:
        {
            scalarField I;

            if (useTable_)
            {
                I.setSize(size());
                forAll(I, facei)
                {
                    const scalar u = max(uStar_[facei], SMALL);
                    I[facei] = wn[facei]/u
                       *(*I1Table_)(ks_*u/nu, delta*u/nu);
                }
            }
            else
            {
                I = wallIntegral(wn, Gamma, delta);
            }

            E_ = wn*Cstar_*exp(I);
            break;
        }

        case modeType::Robin:
        {
            beta_f_ = beta_;
            E_ = beta_*Cstar_;
            break;
        }

        case modeType::nearestCell:
        {
            // E - w_n C_P enters the C equation as a source of the wall
            // cell (CEqn.H); the face itself carries no flux
            E_ = wn*Cstar_;
            refValue() = Zero;
            refGrad() = Zero;
            valueFraction() = Zero;
            mixedFvPatchScalarField::updateCoeffs();
            return;
        }
    }

    const scalarField grad(E_/Gamma);

    // Mixed coefficients: Gamma snGrad(C) + (beta - w_n) C_f = E
    const scalarField K(beta_f_ - wn);
    const scalarField GammaDelta(Gamma*patch().deltaCoeffs());

    if (gMin(K + GammaDelta) <= 0)
    {
        FatalErrorInFunction
            << "beta - w_n + Gamma deltaCoeffs <= 0 on patch "
            << patch().name() << ": deposition velocity too small for this"
            << " mesh (Robin with beta < w_s)"
            << exit(FatalError);
    }

    refValue() = Zero;
    refGrad() = grad;
    valueFraction() = K/(K + GammaDelta);

    mixedFvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::scalarField> Foam::sedimentBedFluxFvPatchScalarField::D() const
{
    if (mode_ == modeType::nearestCell)
    {
        return beta_f_*patchInternalField();
    }

    return beta_f_*(*this);
}


void Foam::sedimentBedFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("mode", modeTypeNames[mode_]);

    os.writeEntry("diffusivity", diffusivityTypeNames[diffusivity_]);
    os.writeEntryIfDifferent<scalar>("kappa", 0.41, kappa_);
    if (mode_ == modeType::HRN)
    {
        os.writeEntry("form", word(exponential_ ? "exponential" : "linear"));
    }
    if (mode_ == modeType::LRN)
    {
        os.writeEntry("integral", word(useTable_ ? "table" : "computed"));
        if (useTable_)
        {
            os.writeEntry("I1Model", I1Model_);
            os.writeEntry("ks", ks_);
            if (I1Model_ == "file")
            {
                os.beginBlock("I1Table");
                I1Table_->write(os);
                os.endBlock();
            }
        }
    }
    if (mode_ == modeType::Robin)
    {
        os.writeEntry("beta", beta_);
    }

    os.writeEntryIfDifferent<word>("DC", "DC", DCName_);
    os.writeEntryIfDifferent<word>("ws", "ws", wsName_);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        sedimentBedFluxFvPatchScalarField
    );
}

// ************************************************************************* //
