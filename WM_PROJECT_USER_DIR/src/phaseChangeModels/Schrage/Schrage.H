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
    Foam::phaseChangeModels::Schrage

Description
    This assumes liquid and vapour are in equilibrium, whereby the
    accommodation coefficients are equivalent for the interface.
    This relation is known as the Hertz-Knudsen-Schrage.

    References:
    \verbatim
        Van P. Carey, Liquid-Vapor Phase Change Phenomena,
        ISBN 0-89116836, 1992, pp. 112-121.
    \endverbatim

SourceFiles
    Schrage.C

\*---------------------------------------------------------------------------*/

#ifndef Schrage_H
#define Schrage_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace phaseChangeModels
{

/*---------------------------------------------------------------------------*\
                         Class Schrage Declaration
\*---------------------------------------------------------------------------*/

class Schrage
:
    public phaseChangeModel
{

public:

    // Public data types

        //- Enumeration defining the method for calculating interface area density
        //  枚举类型可以输出, 但不能输入
        enum interfaceAreaMethods
        {
            gradAlpha,
            isoAlpha,
            isoAdvection
        };


private:

    // Private Data

        //- Phase change accommodation coefficient between 0 and 1
        scalar gamma_;

        //- Mathematical constant pi
        const scalar pi_;

        // "gradAlpha", "isoAlpha", or "isoAdvection"
        word interfaceAreaMethod_;

        //- Interface area per unit volume
        volScalarField interfaceArea_;

        //- Phase change related parameter
        volScalarField phaseChangeCoeff_;


    // Private Member Functions

        //- Phase Change coefficients
        tmp<volScalarField> phaseChangeCoeff() const;


public:
        
    //- Runtime type information
    TypeName("Schrage");

    //- Ideal gas constant with dimensions
    static const dimensionedScalar R0;


    // Constructors

        //- Construct from twoPhaseMixtureThermo
        Schrage
        (
            const twoPhaseMixtureThermo& mixture
        );

        //- Disallow default null constructor
        Schrage() = delete;

        //- Disallow default bitwise copy construction
        Schrage(const Schrage&) = delete;


    //- Destructor
    virtual ~Schrage();


    // Member Functions

        //- Update interface area density
        virtual void updateInterface(const volVectorField& normal);


        //- 通过Pair实现分别计算蒸发的冷凝的相变

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


        //- Correct the phase change properties
        virtual void correct();

        //- Read the phaseChnageModel dictionary and update
        virtual bool read(const dictionary& thermoDict);

        virtual bool isExplicit() const
        {
            return false;
        }

        //- Is phase change model based on the interface area?
        virtual bool isCalcInterface() const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Schrage&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "SchrageI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
