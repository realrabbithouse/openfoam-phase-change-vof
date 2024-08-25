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

#include "Clapeyron.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace saturationTemperatureModels
{
    defineTypeNameAndDebug(Clapeyron, 0);
    addToRunTimeSelectionTable
    (
        saturationTemperatureModel,
        Clapeyron,
        dictionary
    );
} // End namespace Foam

} // End namespace saturationTemperatureModels

const Foam::dimensionedScalar Foam::saturationTemperatureModels::Clapeyron::R0
(
    "R0", dimensionSet(1,2,-2,-1,-1,0,0), RR
);


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationTemperatureModels::Clapeyron::Clapeyron
(
    const twoPhaseMixtureThermo& mixture
)
:
    saturationTemperatureModel(mixture),
    p0_("p0", dimensionSet(1,-1,-2,0,0,0,0), *this),
    T0_("T0", dimensionSet(0,0,0,1,0,0,0), *this)
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationTemperatureModels::Clapeyron::~Clapeyron()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::saturationTemperatureModels::Clapeyron::correct()
{
    latentHeatPtr_->correct();
    
    // Info<< "Correcting saturation temperature by Clausius-Clapeyron equation" << nl;

    if (min(p_).value() > 0)
    {
        const volScalarField& Lv = latentHeatPtr_->Lv();

        // molecular weight of the gas phase
        tmp<volScalarField> tW = mixture_.thermo2().W();
        const volScalarField& W = tW();

        // Info<< "max(Lv): " << max(Lv).value() << " "
        //     << "min(Lv): " << min(Lv).value() << nl << endl;
        // Info<< "p0: " << p0_.value() << " " << "T0: " << T0_.value() << nl << endl;
        // Info<< "max(p): " << max(p_).value() << " "
        //     << "min(p): " << min(p_).value() << nl << endl;
        // Info<< "max(W): " << max(W).value() << " "
        //     << "min(W): " << min(W).value() << nl << endl;

        // 这是因为对返回类型为tmp<volScalarField>的接口的理解错误
        // 很奇怪，当采用以下方式访问时运行时都会报错
        // const volScalarField& W = mixture_.thermo2().W();
        // const volScalarField& W(mixture_.thermo2().W());
        // const volScalarField& W(mixture_.thermo2().W()());

        TSat_ = scalar(1.0)/( scalar(1.0)/T0_ - (R0/W)*log(p_/p0_)/Lv );

        // Info<< "max(TSat): " << max(TSat_).value() << " "
        //     << "min(TSat): " << min(TSat_).value() << nl << endl;
    }
    else
    {
        FatalErrorIn
        (
            "Incorrect Pressure Field"
        )   << "Minimum pressure = " << nl
            << min(p_) << nl << nl;
    }
}

bool Foam::saturationTemperatureModels::Clapeyron::read(const dictionary& thermoDict)
{
    if (saturationTemperatureModel::read(thermoDict))
    {
        const dictionary& satTempDict = thermoDict.subDict
        (
            saturationTemperatureModel::satTempDictName
        );

        satTempDict.lookup("p0") >> p0_;
        satTempDict.lookup("T0") >> T0_;

        return true;
    }

    return false;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
