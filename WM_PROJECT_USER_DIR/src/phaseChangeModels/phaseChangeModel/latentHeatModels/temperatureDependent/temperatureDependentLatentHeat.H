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
    Foam::latentHeatModels::temperatureDependentLatentHeat

Description

SourceFiles
    temperatureDependentLatentHeat.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureDependentLatentHeat_H
#define temperatureDependentLatentHeat_H

#include "latentHeatModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace latentHeatModels
{

/*---------------------------------------------------------------------------*\
                         Class temperatureDependentLatentHeat Declaration
\*---------------------------------------------------------------------------*/

class temperatureDependentLatentHeat
:
    public latentHeatModel
{

protected:

    // Protected Data

        //- Latent-heat function, units [J/kg]
        autoPtr<Function1<scalar>> LvFunc_;


public:

    //- Runtime type information
    TypeName("temperatureDependent");

    // Constructors

        //- No null construct
        temperatureDependentLatentHeat() = delete;

        //- No copy construct
        temperatureDependentLatentHeat(const temperatureDependentLatentHeat&) = delete;

        //- No copy assignment
        temperatureDependentLatentHeat operator=(const temperatureDependentLatentHeat&) = delete;

        //- Construct from components
        temperatureDependentLatentHeat(const twoPhaseMixtureThermo& mixture);


    //- Destructor
    virtual ~temperatureDependentLatentHeat();


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

// #include "temperatureDependentLatentHeatI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
