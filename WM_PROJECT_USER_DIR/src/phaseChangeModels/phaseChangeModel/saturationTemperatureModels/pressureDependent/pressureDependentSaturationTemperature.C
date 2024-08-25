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

#include "pressureDependentSaturationTemperature.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace saturationTemperatureModels
{
    defineTypeNameAndDebug(pressureDependentSaturationTemperature, 0);
    addToRunTimeSelectionTable
    (
        saturationTemperatureModel,
        pressureDependentSaturationTemperature,
        dictionary
    );
} // End namespace Foam

} // End namespace saturationTemperatureModels


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationTemperatureModels::pressureDependentSaturationTemperature::pressureDependentSaturationTemperature
(
    const twoPhaseMixtureThermo& mixture
)
:
    saturationTemperatureModel(mixture),
    TSatFunc_(Function1<scalar>::New("TSat", *this))
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationTemperatureModels::pressureDependentSaturationTemperature::~pressureDependentSaturationTemperature()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::saturationTemperatureModels::pressureDependentSaturationTemperature::correct()
{
    latentHeatPtr_->correct();

    TSat_.field() = TSatFunc_->value(p_.field());

    volScalarField::Boundary& TSatBf = TSat_.boundaryFieldRef();
    const volScalarField::Boundary& pBf = p_.boundaryField();

    forAll(TSatBf, patchi)
    {
        TSatBf[patchi] = TSatFunc_->value(pBf[patchi]);
    }
}

bool Foam::saturationTemperatureModels::pressureDependentSaturationTemperature::read(const dictionary& thermoDict)
{
    if (saturationTemperatureModel::read(thermoDict))
    {
        const dictionary& satTempDict = thermoDict.subDict
        (
            saturationTemperatureModel::satTempDictName
        );

        TSatFunc_ = Function1<scalar>::New("TSat", satTempDict);

        return true;
    }
    
    return false;
}


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
