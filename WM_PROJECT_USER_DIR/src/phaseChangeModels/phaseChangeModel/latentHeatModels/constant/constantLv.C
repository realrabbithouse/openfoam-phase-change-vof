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

#include "constantLv.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace latentHeatModels
{
    defineTypeNameAndDebug(constantLv, 0);
    addToRunTimeSelectionTable
    (
        latentHeatModel,
        constantLv,
        dictionary
    );
} // End namespace Foam

} // End namespace latentHeatModels


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::latentHeatModels::constantLv::constantLv
(
    const twoPhaseMixtureThermo& mixture
)
:
    latentHeatModel(mixture),
    constLv_("Lv", dimensionSet(0,2,-2,0,0,0,0), *this) // *this is a dictionary
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::latentHeatModels::constantLv::~constantLv()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::latentHeatModels::constantLv::correct()
{
    Lv_ = constLv_;
}

bool Foam::latentHeatModels::constantLv::read(const dictionary& thermoDict)
{
    if (latentHeatModel::read(thermoDict))
    {
        const dictionary& latentHeatDict = thermoDict.subDict
        (
            latentHeatModel::latentHeatDictName
        );

        latentHeatDict.lookup("Lv") >> constLv_;

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
