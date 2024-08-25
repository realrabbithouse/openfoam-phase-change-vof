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

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeModel, 0);
    defineRunTimeSelectionTable(phaseChangeModel, dictionary);
}

const Foam::word Foam::phaseChangeModel::phaseChangeDictName("phaseChangeProperties");

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModel::phaseChangeModel
(
    const twoPhaseMixtureThermo& mixture
)
:
    dictionary(mixture.subDict(phaseChangeDictName)),
    mixture_(mixture),
    mesh_(mixture.T().mesh()),
    mDot_
    (
        IOobject
        (
            "mDot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),
    saturationTemperaturePtr_
    (
        saturationTemperatureModel::New(mixture)
    )
{
    // empty
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseChangeModel::~phaseChangeModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::saturationTemperatureModel&
Foam::phaseChangeModel::saturationTemperature() const
{
    return saturationTemperaturePtr_();
}


const Foam::latentHeatModel&
Foam::phaseChangeModel::latentHeat() const
{
    return saturationTemperaturePtr_->latentHeat();
}


void Foam::phaseChangeModel::updateInterface
(
    const volVectorField& normal
)
{
    NotImplemented;

    return;
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModel::vDotAlphal() const
{
    // temperaturePhaseChangeTwoPhaseMixture并没有limitedAlpha1, 而是
    // 直接使用mixture_.alpha1
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField alphalCoeff
    (
        scalar(1)/mixture_.thermo1().rho() - limitedAlpha1
        *(scalar(1)/mixture_.thermo1().rho() - scalar(1)/mixture_.thermo2().rho())
    );

    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*mDotAlphal[0], // condensation, >= 0
        alphalCoeff*mDotAlphal[1] // evaporation, <= 0
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModel::vDot() const
{
    volScalarField pCoeff // < 0
    (
        scalar(1)/mixture_.thermo1().rho() - scalar(1)/mixture_.thermo2().rho()
    );

    Pair<tmp<volScalarField>> mDot = this->mDot();

    return Pair<tmp<volScalarField>>
    (
        pCoeff*mDot[0], // condensation <= 0
        pCoeff*mDot[1] // evaporation >= 0
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>> // [W/m^3]
Foam::phaseChangeModel::qDot() const
{
    const volScalarField& Lv = this->Lv();
    Pair<tmp<volScalarField>> mDot = this->mDot();

    return Pair<tmp<volScalarField>>
    (
        Lv*mDot[0],
        Lv*mDot[1]
    );
}


// explicit phase change source
Foam::tmp<Foam::volScalarField> Foam::phaseChangeModel::Salphal() const
{
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField alphalCoeff
    (
        scalar(1)/mixture_.thermo1().rho() - limitedAlpha1
        *(scalar(1)/mixture_.thermo1().rho() - scalar(1)/mixture_.thermo2().rho())
    );

    return alphalCoeff*mDot_;
}


Foam::Pair<Foam::tmp<Foam::volScalarField>> Foam::phaseChangeModel::alphalSuSp() const
{
    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField SuCoeff // > 0
    (
        scalar(1)/mixture_.thermo1().rho()
    );

    volScalarField SpCoeff // > 0
    (
        scalar(1)/mixture_.thermo2().rho() - scalar(1)/mixture_.thermo1().rho()
    );

    volScalarField alphalCoeff // > 0
    (
        scalar(1)/mixture_.thermo1().rho() - limitedAlpha1
        *(scalar(1)/mixture_.thermo1().rho() - scalar(1)/mixture_.thermo2().rho())
    );

    tmp<volScalarField> talphalSu
    (
        volScalarField::New
        (
            "alphalSu",
            mesh_,
            dimensionedScalar(dimless/dimTime, Zero)
        )
    );

    tmp<volScalarField> talphalSp
    (
        volScalarField::New
        (
            "alphalSp",
            mesh_,
            dimensionedScalar(dimless/dimTime, Zero)
        )
    );

    volScalarField& alphalSu = talphalSu.ref();
    volScalarField& alphalSp = talphalSp.ref();

    forAll(mDot_, celli)
    {
        if (mDot_[celli] > 0.0)
        {
            alphalSu[celli] = alphalCoeff[celli]*mDot_[celli];
            alphalSp[celli] = 0.0;
        }
        else if (mDot_[celli] < 0.0)
        {
            alphalSu[celli] = SuCoeff[celli]*mDot_[celli];
            alphalSp[celli] = SpCoeff[celli]*mDot_[celli];
        }
    }

    return Pair<tmp<volScalarField>>
    (
        talphalSu,
        talphalSp
    );
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeModel::Sv() const
{
    volScalarField pCoeff // < 0
    (
        scalar(1)/mixture_.thermo1().rho() - scalar(1)/mixture_.thermo2().rho()
    );

    // generally, negetive for condensation, positive for evaporation
    return pCoeff*mDot_;
}


Foam::tmp<Foam::volScalarField> Foam::phaseChangeModel::Sq() const
{
    const volScalarField& Lv = this->Lv();

    return Lv*mDot_;
}


bool Foam::phaseChangeModel::read(const dictionary& thermoDict)
{
    if (saturationTemperaturePtr_->read(thermoDict))
    {
        // const dictionary& phaseChangeDict = thermoDict.subDict
        // (
        //     phaseChangeDictName
        // );

        // phaseChangeDict.readIfPresent<Switch>
        // (
        //     "isExplicit", isExplicit_
        // );

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
