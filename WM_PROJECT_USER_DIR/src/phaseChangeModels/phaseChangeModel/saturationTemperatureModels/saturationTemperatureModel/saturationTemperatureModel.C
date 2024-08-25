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

\*---------------------------------------------------------------------------*/

#include "saturationTemperatureModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(saturationTemperatureModel, 0);
    defineRunTimeSelectionTable(saturationTemperatureModel, dictionary);
}

const Foam::word Foam::saturationTemperatureModel::satTempDictName("saturationTemperatureProperties");

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationTemperatureModel::saturationTemperatureModel
(
    const twoPhaseMixtureThermo& mixture
)
:
    dictionary(mixture.subDict(satTempDictName)),
    mixture_(mixture),
    mesh_(mixture.T().mesh()),
    p_(mixture.p()),
    latentHeatPtr_
    (
        latentHeatModel::New(mixture)
    ),
    TSat_
    (
        IOobject
        (
            "TSat",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0,0,0,1,0,0,0), Zero)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationTemperatureModel::~saturationTemperatureModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::latentHeatModel& Foam::saturationTemperatureModel::latentHeat() const
{
    return latentHeatPtr_();
}

bool Foam::saturationTemperatureModel::read(const dictionary& thermoDict)
{
    if (latentHeatPtr_->read(thermoDict))
    {   
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
