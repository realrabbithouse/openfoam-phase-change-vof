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
    Foam::staTempModels::Clapeyron

Description

SourceFiles
    Clapeyron.C

\*---------------------------------------------------------------------------*/

#ifndef Clapeyron_H
#define Clapeyron_H

#include "saturationTemperatureModel.H"
#include "thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace saturationTemperatureModels
{

/*---------------------------------------------------------------------------*\
                         Class Clapeyron Declaration
\*---------------------------------------------------------------------------*/

class Clapeyron
:
    public saturationTemperatureModel
{

protected:

    // protected Data

        //- Reference saturation prssure
        dimensionedScalar p0_;

        //- Reference saturation temperature
        dimensionedScalar T0_;


public:

    //- Runtime type information
    TypeName("Clapeyron");

    //- Ideal gas constant with dimensions
    static const dimensionedScalar R0;


    // Constructors

        //- Construct from components
        Clapeyron(const twoPhaseMixtureThermo& mixture);

        //- Disallow null constructor
        Clapeyron() = delete;

        //- Disallow default bitwise copy construction
        Clapeyron(const Clapeyron&) = delete;

        //- Disallow default bitwise assignment
        Clapeyron& operator=(const Clapeyron&) = delete;


    //- Destructor
    virtual ~Clapeyron();


    // Member Functions

        //- Correct saturation temperature
        virtual void correct();

        //- Read saturation properties from dictionary
        virtual bool read(const dictionary& thermoDict);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationTemperatureModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "ClapeyronI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
