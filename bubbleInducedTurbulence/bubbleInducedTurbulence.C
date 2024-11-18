/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |

\*---------------------------------------------------------------------------*/

#include "bubbleInducedTurbulence.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSST::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*this->nu()/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSST::
F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*this->nu()/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSST::
F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*this->nu()/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSST::
F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
void kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
void kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::
correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), F23());
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    return betaStar_*omega_();
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::
omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSST
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    MomentumTransportModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    ci_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ci",
            this->coeffDict_,
            0.4
        )
    ),
    ce_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce",
            this->coeffDict_,
            0.17
        )
    ),
    cd_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cd",
            this->coeffDict_,
            0.71
        )
    ),
    dp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "dp",
            this->coeffDict_,
            0.01
        )
    ),
    cOmegaD_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cOmegaD",
            this->coeffDict_,
            0.0828
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
bool kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::read()
{
    if (MomentumTransportModel::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

// dh-add
template<class MomentumTransportModel, class BasicMomentumTransportModel>
const phaseCompressibleMomentumTransportModel& kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::gasTurbulence() const
{
    return this->gasTurbulence();
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
const phaseModel& kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::liquidPhase() const
{
    return refCast<const phaseModel>(this->properties());
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField> kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::relativeVelocityMagnitude() const
{
    // Get the gas and liquid phase velocities
    const phaseCompressibleMomentumTransportModel& gasTurbulence = this->gasTurbulence();
    const volVectorField& gasVelocity = gasTurbulence.U();
    const volVectorField& liquidVelocity = this->U_;

    // Compute and store the magnitude of the relative velocity as ur
    tmp<volScalarField> ur = mag(liquidVelocity - gasVelocity);

    
    return ur ;
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
void kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    MomentumTransportModel::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu(dev(twoSymm(tgradU()())) && tgradU()());
    volScalarField::Internal G(this->GName(), nut()*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());
    volScalarField ur_magnitude = mag(relativeVelocityMagnitude());
    tmp<volScalarField> fd = (3.0/(4.0*dp_))*cd_*rho*alpha*relativeVelocityMagnitude()*ur_magnitude;

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma
           *min
            (
                GbyNu,
                (c1_/a1_)*betaStar_*omega_()
               *max(a1_*omega_(), b1_*F23()*sqrt(S2()))
            )
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*cOmegaD_*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvModels.source(alpha, rho, omega_)
          + (1.0/(betaStar_*k_))*(ce_*ci_*relativeVelocityMagnitude()*fd*relativeVelocityMagnitude()/dp_)
          - (omega_/k_)*ci_*relativeVelocityMagnitude()*fd
        );

        omegaEqn.ref().relax();
        fvConstraints.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvConstraints.constrain(omega_);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha()*rho()*Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*betaStar_*rho()*omega_*k_, k_) //dissipation term
      + kSource()
      + ci_*relativeVelocityMagnitude()*fd
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

    correctNut(S2, F23);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
