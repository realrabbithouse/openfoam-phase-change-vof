/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Class
    Foam::compressibleInterPhaseTransportModel

Description
    Transport model selection class for the compressibleInterFoam family of
    solvers.

    By default the standard mixture transport modelling approach is used in
    which a single momentum stress model (laminar, non-Newtonian, LES or RAS) is
    constructed for the mixture.  However if the \c simulationType in
    constant/turbulenceProperties is set to \c twoPhaseTransport the alternative
    Euler-Euler two-phase transport modelling approach is used in which separate
    stress models (laminar, non-Newtonian, LES or RAS) are instantiated for each
    of the two phases allowing for different modeling for the phases.
    

SourceFiles
    compressibleInterPhaseTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleInterPhaseTransportModel_H
#define compressibleInterPhaseTransportModel_H

#include "twoPhaseMixture.H"
#include "twoPhaseMixtureThermo.H"
#include "turbulentFluidThermoModel.H"
#include "VoFphaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class compressibleInterPhaseTransportModel Declaration
\*---------------------------------------------------------------------------*/

class compressibleInterPhaseTransportModel
{
    // Private data

        //- Switch to select two-phase or mixture transport modelling
        Switch twoPhaseTransport_;

        //- Two-phase mixture
        const twoPhaseMixtureThermo& mixture_;

        //- Mixture volumetric flux
        const surfaceScalarField& phi_;

        //- Phase volumetric flux
        const surfaceScalarField& alphaPhi10_;

        //- Phase-1 mass-flux (constructed for two-phase transport)
        tmp<surfaceScalarField> alphaRhoPhi1_;

        //- Phase-2 mass-flux (constructed for two-phase transport)
        tmp<surfaceScalarField> alphaRhoPhi2_;

        //- Mixture transport model (constructed for mixture transport)
        autoPtr<compressible::turbulenceModel> turbulence_;

        //- Phase-1 transport model (constructed for two-phase transport)
        autoPtr
        <
            ThermalDiffusivity<PhaseCompressibleTurbulenceModel<fluidThermo>>
        > turbulence1_;

        //- Phase-2 transport model (constructed for two-phase transport)
        autoPtr
        <
            ThermalDiffusivity<PhaseCompressibleTurbulenceModel<fluidThermo>>
        > turbulence2_;


    // Private Member Functions

        //- No copy construct
        compressibleInterPhaseTransportModel
        (
            const compressibleInterPhaseTransportModel&
        ) = delete;

        //- No copy assignment
        void operator=(const compressibleInterPhaseTransportModel&) = delete;


public:

    // Constructors

        //- Construct from components
        compressibleInterPhaseTransportModel
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const surfaceScalarField& rhoPhi,
            const surfaceScalarField& alphaPhi10,
            const twoPhaseMixtureThermo& mixture
        );


    // Member Functions

        //- Return the effective temperature transport coefficient
        tmp<volScalarField> alphaEff() const;

        //- Return the effective thermal conductivity coefficient
        tmp<volScalarField> kappaEff() const;

        //- Return the effective momentum stress divergence
        tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Correct the phase mass-fluxes
        // (required for the two-phase transport option)
        void correctPhasePhi();

        //- Correct the phase or mixture transport models
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
