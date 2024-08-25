/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::phaseChangeModels::noPhaseChange

Description

SourceFiles
    noPhaseChange.C

\*---------------------------------------------------------------------------*/

#ifndef noPhaseChange_H
#define noPhaseChange_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace phaseChangeModels
{

/*---------------------------------------------------------------------------*\
                         Class noPhaseChange Declaration
\*---------------------------------------------------------------------------*/

class noPhaseChange
:
    public phaseChangeModel
{

public:

    //- Runtime type information
    TypeName("noPhaseChange");


    // Constructors

        //- Construct from twoPhaseMixtureThermo
        noPhaseChange(const twoPhaseMixtureThermo& mixture);

        //- Disallow null constructor
        noPhaseChange() = delete;

        //- Disallow default bitwise copy construction
        noPhaseChange(const noPhaseChange&) = delete;


    //- Destructor
    virtual ~noPhaseChange();


    // Member Functions

        //- Correct phase change model
        void correct()
        {
            saturationTemperaturePtr_->correct();
        }

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate [kg/(s.m^3)]
        virtual Pair<tmp<volScalarField>> mDotAlphal() const;

        //- Return the mass condensation and vaporisation rates as coefficients [kg/(s.m^3)]
        virtual Pair<tmp<volScalarField>> mDot() const;

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (Tsat - T) for the condensation rate
        //  and a coefficient to multiply  (T - Tsat) for the vaporisation rate
        //  mDot = mDotDeltaT()[0]*(TSat - T) + mDotDeltaT()[1]*(TSat - T)
        //  mDotDeltaT()[0]*(TSat - T) >= 0; mDotDeltaT()[1]*(TSat - T) <= 0
        virtual Pair<tmp<volScalarField>> mDotDeltaT() const;

        //- Source for T equation, implicit energy source
        virtual tmp<fvScalarMatrix> TSource() const;

        //- Read phase change dictionary and update
        virtual bool read(const dictionary& thermoDict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const noPhaseChange&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
