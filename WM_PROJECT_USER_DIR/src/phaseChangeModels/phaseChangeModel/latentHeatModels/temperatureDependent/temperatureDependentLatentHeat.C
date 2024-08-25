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

#include "temperatureDependentLatentHeat.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace latentHeatModels
{
    defineTypeNameAndDebug(temperatureDependentLatentHeat, 0);
    addToRunTimeSelectionTable
    (
        latentHeatModel,
        temperatureDependentLatentHeat,
        dictionary
    );
} // End namespace Foam

} // End namespace latentHeatModels


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::latentHeatModels::temperatureDependentLatentHeat::temperatureDependentLatentHeat
(
    const twoPhaseMixtureThermo& mixture
)
:
    latentHeatModel(mixture),
    LvFunc_(Function1<scalar>::New("Lv", *this))
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::latentHeatModels::temperatureDependentLatentHeat::~temperatureDependentLatentHeat()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::latentHeatModels::temperatureDependentLatentHeat::correct()
{
    const volScalarField& T = mixture_.T();

    Lv_.field() = LvFunc_->value(T.field());

    volScalarField::Boundary& LvBf = Lv_.boundaryFieldRef();
    const volScalarField::Boundary& TBf = T.boundaryField();

    forAll(LvBf, patchi)
    {
        LvBf[patchi] = LvFunc_->value(TBf[patchi]);
    }
}

bool Foam::latentHeatModels::temperatureDependentLatentHeat::read(const dictionary& thermoDict)
{
    if (latentHeatModel::read(thermoDict))
    {
        const dictionary& latentHeatDict = thermoDict.subDict
        (
            latentHeatModel::latentHeatDictName
        );

        LvFunc_ = Function1<scalar>::New("Lv", latentHeatDict);

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
