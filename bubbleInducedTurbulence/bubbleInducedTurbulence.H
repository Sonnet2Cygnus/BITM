/*---------------------------------------------------------------------------*\  
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef bubbleInducedTurbulence_H
#define bubbleInducedTurbulence_H

#include "dimensionedScalar.H"
#include "volFields.H"
#include "Switch.H"
#include "fvc.H"
#include "fvm.H"
#include "viscosityModel.H"
#include "physicalProperties.H"
#include "phaseModel.H"
#include "phaseSystem.H"
#include "phaseCompressibleMomentumTransportModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class Declaration
\*---------------------------------------------------------------------------*/
template<class MomentumTransportModel, class BasicMomentumTransportModel>
class kOmegaSST : public MomentumTransportModel
{
protected:

    // ========================================================================
    // Model Coefficients
    // ========================================================================
    dimensionedScalar alphaK1_, alphaK2_;
    dimensionedScalar alphaOmega1_, alphaOmega2_;
    dimensionedScalar gamma1_, gamma2_;
    dimensionedScalar beta1_, beta2_, betaStar_;
    dimensionedScalar a1_, b1_, c1_;
    Switch F3_;
    dimensionedScalar ci_, ce_, cd_, dp_, cOmegaD_;

    // ========================================================================
    // Model Fields
    // ========================================================================
    const volScalarField& y_;           // Wall distance
    volScalarField k_, omega_;          // Turbulence variables

    // ========================================================================
    // Protected Member Functions
    // ========================================================================

    // Control Coefficients
    virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
    virtual tmp<volScalarField> F2() const;
    virtual tmp<volScalarField> F3() const;
    virtual tmp<volScalarField> F23() const;

    // General Blending Function
    template<typename FieldType>
    tmp<FieldType> blend
    (
        const FieldType& F1,
        const dimensionedScalar& psi1,
        const dimensionedScalar& psi2
    ) const
    {
        return F1 * (psi1 - psi2) + psi2;
    }

    tmp<volScalarField> alphaK(const volScalarField& F1) const;
    tmp<volScalarField> alphaOmega(const volScalarField& F1) const;
    tmp<volScalarField::Internal> beta(const volScalarField::Internal& F1) const;
    tmp<volScalarField::Internal> gamma(const volScalarField::Internal& F1) const;

    // Nut Correction
    virtual void correctNut(const volScalarField& S2, const volScalarField& F2);
    virtual void correctNut();

    // Turbulence Production and Dissipation
    virtual tmp<volScalarField::Internal> Pk(const volScalarField::Internal& G) const;
    virtual tmp<volScalarField::Internal> epsilonByk
    (
        const volScalarField::Internal& F1,
        const volScalarField::Internal& F2
    ) const;

    // Source Terms
    virtual tmp<fvScalarMatrix> kSource() const;
    virtual tmp<fvScalarMatrix> omegaSource() const;
    virtual tmp<fvScalarMatrix> Qsas
    (
        const volScalarField::Internal& S2,
        const volScalarField::Internal& gamma,
        const volScalarField::Internal& beta
    ) const;

    // Phase Properties
    const phaseCompressibleMomentumTransportModel& gasTurbulence() const;
    const phaseModel& liquidPhase() const;

public:

    // ========================================================================
    // Type Aliases
    // ========================================================================
    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;

    // ========================================================================
    // Constructors and Destructor
    // ========================================================================
    kOmegaSST
    (
        const word& type,
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const viscosity& viscosity
    );

    kOmegaSST(const kOmegaSST&) = delete;

    virtual ~kOmegaSST() {}

    // ========================================================================
    // Public Member Functions
    // ========================================================================
    virtual bool read();                         // Update model coefficients
    tmp<volScalarField> DkEff(const volScalarField& F1) const;  // Diffusivity for k
    tmp<volScalarField> DomegaEff(const volScalarField& F1) const;  // Diffusivity for omega
    tmp<volScalarField> relativeVelocityMagnitude() const;      // Relative velocity
    tmp<volScalarField> calculateFd() const;    // Calculate additional drag

    // Turbulence Properties
    virtual tmp<volScalarField> k() const;      // Return turbulence kinetic energy
    virtual tmp<volScalarField> epsilon() const; // Return dissipation rate
    virtual tmp<volScalarField> omega() const;  // Return turbulence frequency
    virtual void correct();                     // Correct turbulence properties

    // Operators
    void operator=(const kOmegaSST&) = delete;
};

} 

#endif 

// ************************************************************************* //
