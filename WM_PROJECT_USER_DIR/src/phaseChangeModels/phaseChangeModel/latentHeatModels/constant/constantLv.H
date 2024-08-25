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
    Foam::latentHeatModels::constantLv

Description

SourceFiles
    constantLv.C

\*---------------------------------------------------------------------------*/

#ifndef constantLv_H
#define constantLv_H

#include "latentHeatModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace latentHeatModels
{

/*---------------------------------------------------------------------------*\
                         Class constantLv Declaration
\*---------------------------------------------------------------------------*/

class constantLv
:
    public latentHeatModel
{

protected:

    // Protected Data

        //- Const specific latent heat of vaporization, units [J/kg]
        dimensionedScalar constLv_;


public:

    //- Runtime type information
    TypeName("constant");

    // Constructors

        //- No null construct
        constantLv() = delete;

        //- No copy construct
        constantLv(const constantLv&) = delete;

        //- No copy assignment
        constantLv operator=(const constantLv&) = delete;

        //- Construct from components
        constantLv(const twoPhaseMixtureThermo& mixture);


    //- Destructor
    virtual ~constantLv();


    // Member Functions

        // virtual tmp<volScalarField> Lv() const;

        virtual void correct();

        virtual bool read(const dictionary& thermoDict);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace latentHeatModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "constantLvI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
