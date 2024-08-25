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

#include "constantTSat.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace saturationTemperatureModels
{
    defineTypeNameAndDebug(constantTSat, 0);
    addToRunTimeSelectionTable
    (
        saturationTemperatureModel,
        constantTSat,
        dictionary
    );
} // End namespace Foam

} // End namespace saturationTemperatureModels


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationTemperatureModels::constantTSat::constantTSat
(
    const twoPhaseMixtureThermo& mixture
)
:
    saturationTemperatureModel(mixture),
    constTSat_("TSat", dimensionSet(0,0,0,1,0,0,0), *this)
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationTemperatureModels::constantTSat::~constantTSat()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::saturationTemperatureModels::constantTSat::correct()
{
    latentHeatPtr_->correct();
    
    TSat_ = constTSat_;
}

bool Foam::saturationTemperatureModels::constantTSat::read(const dictionary& thermoDict)
{
    if (saturationTemperatureModel::read(thermoDict))
    {
        const dictionary& satTempDict = thermoDict.subDict
        (
            saturationTemperatureModel::satTempDictName
        );

        satTempDict.lookup("TSat") >> constTSat_;

        return true;
    }
    
    return false;
}


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
