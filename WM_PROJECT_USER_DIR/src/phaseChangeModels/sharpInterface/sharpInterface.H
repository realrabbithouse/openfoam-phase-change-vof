/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 AUTHOR,AFFILIATION
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
    Foam::sharpInterface

Description
    grad(alphal) & grad(T) > 0 ---> evaporation
    grad(alphal) & grad(T) < 0 ---> condensation
    事实上只要满足相界面处温度梯度非零，则必然会产生相变传质，与饱和温度没有关系，该单元的
    相体积分数显然介于0到1之间
    
SourceFiles
    sharpInterface.C

\*---------------------------------------------------------------------------*/

#ifndef sharpInterface_H
#define sharpInterface_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace phaseChangeModels
{

/*---------------------------------------------------------------------------*\
                         Class sharpInterface Declaration
\*---------------------------------------------------------------------------*/

class sharpInterface
:
    public phaseChangeModel
{
    // Private Data

        //- Interface area density
        volVectorField interfaceArea_;

        // "gradAlpha", "isoAlpha", or "isoAdvection"
        word interfaceAreaMethod_;

        //- The flux field determined from "grad(T) & interfaceArea_"
        volScalarField flux_;

        //- Explicit or implicit evaluation
        Switch isExplicit_;


    // Private Member Functions

        //- Disallow default null constructor
        sharpInterface() = delete;

        //- No copy construct
        sharpInterface(const sharpInterface&) = delete;

        //- No copy assignment
        void operator=(const sharpInterface&) = delete;


public:

    //- Runtime type information
    TypeName("sharpInterface");


    // Constructors

        //- Construct from twoPhaseMixtureThermo
        sharpInterface(const twoPhaseMixtureThermo& mixture);


    //- Destructor
    virtual ~sharpInterface();


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

        //- Explicit energy source
        virtual tmp<volScalarField> Sq() const;


        //- Correct the phase change properties
        virtual void correct();

        //- Read the phaseChnageModel dictionary and update
        virtual bool read(const dictionary& thermoDict);

        //- Is phase change model based on the interface area?
        virtual bool isCalcInterface() const
        {
            return true;
        }

        virtual bool isExplicit() const
        {
            if (isExplicit_)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "sharpInterfaceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
