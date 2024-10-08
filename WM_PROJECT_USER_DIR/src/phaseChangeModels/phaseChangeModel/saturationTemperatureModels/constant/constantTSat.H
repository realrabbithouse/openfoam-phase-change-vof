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
    Foam::saturationTemperatureModels::constantTSat

Description

SourceFiles
    constantTSat.C

\*---------------------------------------------------------------------------*/

#ifndef constantTSat_H
#define constantTSat_H

#include "saturationTemperatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace saturationTemperatureModels
{

/*---------------------------------------------------------------------------*\
                         Class constantTSat Declaration
\*---------------------------------------------------------------------------*/

class constantTSat
:
    public saturationTemperatureModel
{

protected:

    // Protected Data

        //- Saturation temperature, units [T]
        dimensionedScalar constTSat_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from components
        constantTSat(const twoPhaseMixtureThermo& mixture);

        //- Disallow null constructor
        constantTSat() = delete;

        //- Disallow default bitwise copy construction
        constantTSat(const constantTSat&) = delete;

        //- Disallow default bitwise assignment
        constantTSat& operator=(const constantTSat&) = delete;


    //- Destructor
    virtual ~constantTSat();


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

// #include "constantTSatI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
