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

\*---------------------------------------------------------------------------*/

#include "latentHeatModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(latentHeatModel, 0);
    defineRunTimeSelectionTable(latentHeatModel, dictionary);
}

const Foam::word Foam::latentHeatModel::latentHeatDictName("latentHeatProperties");

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::latentHeatModel::latentHeatModel(const twoPhaseMixtureThermo& mixture)
:
    dictionary(mixture.subDict(latentHeatDictName)),
    mixture_(mixture),
    mesh_(mixture.T().mesh()),
    Lv_
    (
        IOobject
        (
            "Lv",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0,2,-2,0,0,0,0), Zero)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::latentHeatModel::~latentHeatModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::latentHeatModel::read(const dictionary& thermoDict)
{
    return true;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
