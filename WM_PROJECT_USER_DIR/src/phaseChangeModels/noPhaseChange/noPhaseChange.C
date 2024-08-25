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

#include "noPhaseChange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace phaseChangeModels
{
    defineTypeNameAndDebug(noPhaseChange, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeModel,
        noPhaseChange,
        dictionary
    );
}

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::noPhaseChange::noPhaseChange
(
    const twoPhaseMixtureThermo& mixture
)
:
    phaseChangeModel(mixture)
{
    // empty
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseChangeModels::noPhaseChange::~noPhaseChange()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::noPhaseChange::mDotAlphal() const
{
    return Pair<tmp<volScalarField>>
    (
        volScalarField::New
        (
            "mDotcAlphal",
            mesh_,
            dimensionedScalar(dimensionSet(1,-3,-1,0,0,0,0), Zero)
        ),
        volScalarField::New
        (
            "mDotvAlphal",
            mesh_,
            dimensionedScalar(dimensionSet(1,-3,-1,0,0,0,0), Zero)
        )
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::noPhaseChange::mDot() const
{
    return Pair<tmp<volScalarField>>
    (
        volScalarField::New
        (
            "mDotc",
            mesh_,
            dimensionedScalar(dimensionSet(1,-3,-1,0,0,0,0), Zero)
        ),
        volScalarField::New
        (
            "mDotv",
            mesh_,
            dimensionedScalar(dimensionSet(1,-3,-1,0,0,0,0), Zero)
        )
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::noPhaseChange::mDotDeltaT() const
{
    return Pair<tmp<volScalarField>>
    (
        volScalarField::New
        (
            "mDotcDeltaT",
            mesh_,
            dimensionedScalar(dimensionSet(1,-3,-1,-1,0,0,0), Zero)
        ),
        volScalarField::New
        (
            "mDotvDeltaT",
            mesh_,
            dimensionedScalar(dimensionSet(1,-3,-1,-1,0,0,0), Zero)
        )
    );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::phaseChangeModels::noPhaseChange::TSource() const
{
    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();

    tmp<fvScalarMatrix> tTSource
    (
        new fvScalarMatrix
        (
            T,
            dimEnergy/dimTime // 焦耳/秒
        )
    );

    fvScalarMatrix& TSource = tTSource.ref();

    const volScalarField zeroCoeff
    (
        volScalarField::New
        (
            "zeroCoeff",
            mesh_,
            dimensionedScalar(dimensionSet(1,-1,-3,-1,0,0,0), Zero)
        )
    );

    TSource =
        fvm::Sp(zeroCoeff, T) - zeroCoeff*TSat
      + fvm::Sp(zeroCoeff, T) - zeroCoeff*TSat;

    return tTSource;

    /*
    const volScalarField& T = mixture_.T();

    tmp<fvScalarMatrix> tTSource
    (
        new fvScalarMatrix
        (
            T,
            dimEnergy/dimTime // 焦耳/秒
        )
    );

    return tTSource;
    */
}


bool Foam::phaseChangeModels::noPhaseChange::read
(
    const dictionary& thermoDict
)
{
    // basic Properties must be provided, although they have never been used
    if (phaseChangeModel::read(thermoDict))
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
