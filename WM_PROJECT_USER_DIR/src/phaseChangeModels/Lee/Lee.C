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

#include "Lee.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace phaseChangeModels
{
    defineTypeNameAndDebug(Lee, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeModel,
        Lee,
        dictionary
    );
}

}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::Lee::Lee
(
    const twoPhaseMixtureThermo& mixture
)
:
    phaseChangeModel(mixture),
    rc_("rc", dimensionSet(0,0,-1,0,0,0,0), *this), // condensation
    rv_("rv", dimensionSet(0,0,-1,0,0,0,0), *this) // evaporation
{
    // empty
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseChangeModels::Lee::~Lee()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::Lee::mDotAlphal() const
{
    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();

    const dimensionedScalar T0(dimTemperature, Zero);

    return Pair<tmp<volScalarField>>
    (
        rc_*mixture_.thermo2().rho()*(max(TSat - T.oldTime(), T0)/TSat), // > 0, condensation
        -rv_*mixture_.thermo1().rho()*(max(T.oldTime() - TSat, T0)/TSat) // < 0, evaporation
    );

    // if (phaseChangeSwitch_)
    // {
    //    /.../
    // }
    // else
    // {
    //     return Pair<tmp<volScalarField>>
    //     (
    //         volScalarField::New
    //         (
    //             "mDotcAlphal",
    //             mesh_,
    //             dimensionedScalar(dimensionSet(1,-3,-1,0,0,0,0), Zero)
    //         ),
    //         volScalarField::New
    //         (
    //             "mDotvAlphal",
    //             mesh_,
    //             dimensionedScalar(dimensionSet(1,-3,-1,0,0,0,0), Zero)
    //         )
    //     );
    // }
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::Lee::mDot() const
{
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();

    const dimensionedScalar T0(dimTemperature, Zero);

    return Pair<tmp<volScalarField>>
    (
        rc_*mixture_.thermo2().rho()*(scalar(1) - limitedAlpha1)
        *(max(TSat - T.oldTime(), T0)/TSat), // > 0, condensation
        -rv_*mixture_.thermo1().rho()*limitedAlpha1
        *(max(T.oldTime() - TSat, T0)/TSat) // < 0, evaporation
    );
}


void Foam::phaseChangeModels::Lee::correct()
{
    saturationTemperaturePtr_->correct();

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();

    const dimensionedScalar T0(dimTemperature, Zero);

    mDot_ =
        rc_*mixture_.thermo2().rho()*(scalar(1) - limitedAlpha1)
        *( max(TSat - T.oldTime(), T0)/TSat ) // condensation, > 0
      + rv_*mixture_.thermo1().rho()*limitedAlpha1
        *( min(TSat - T.oldTime(), T0)/TSat ); // evaporation, < 0
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::Lee::mDotDeltaT() const
{
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();

    return Pair<tmp<volScalarField>>
    (
        rc_*mixture_.thermo2().rho()*(scalar(1) - limitedAlpha1)*pos(TSat - T.oldTime())/TSat,
        rv_*mixture_.thermo1().rho()*limitedAlpha1*pos(T.oldTime() - TSat)/TSat
    );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::phaseChangeModels::Lee::TSource() const
{
    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();
    const volScalarField& Lv = this-> Lv();

    tmp<fvScalarMatrix> tTSource
    (
        new fvScalarMatrix
        (
            T,
            dimEnergy/dimTime // 焦耳/秒
        )
    );

    fvScalarMatrix& TSource = tTSource.ref();

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    const volScalarField Ccoeff
    (
        rc_*mixture_.thermo2().rho()*(scalar(1) - limitedAlpha1)*Lv*pos(TSat - T)/TSat
    );

    const volScalarField Vcoeff
    (
        rv_*mixture_.thermo1().rho()*limitedAlpha1*Lv*pos(T - TSat)/TSat
    );

    // 注: coeff*TSat的单位是J/m^3/s, 而不是J/s
    // condensation <= 0, evaporation >= 0, 所以当liquid作为主相时
    // TSource应当放在能量方程的左边, 放在右边的话需要添加负号
    TSource =
        fvm::Sp(Vcoeff, T) - Vcoeff*TSat
      + fvm::Sp(Ccoeff, T) - Ccoeff*TSat;

    return tTSource;
}


bool Foam::phaseChangeModels::Lee::read
(
    const dictionary& thermoDict
)
{
    if (phaseChangeModel::read(thermoDict))
    {
        const dictionary& phaseChangeDict = thermoDict.subDict
        (
            phaseChangeModel::phaseChangeDictName
        );

        phaseChangeDict.readEntry("rc", rc_);
        phaseChangeDict.readEntry("rv", rv_);

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
