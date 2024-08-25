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

#include "interfaceEquilibrium.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace phaseChangeModels
{
    defineTypeNameAndDebug(interfaceEquilibrium, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeModel,
        interfaceEquilibrium,
        dictionary
    );
}

}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::interfaceEquilibrium::interfaceEquilibrium
(
    const twoPhaseMixtureThermo& mixture
)
:
    phaseChangeModel(mixture),
    global_(lookupOrDefault<Switch>("global", "no")),
    oldFashion_(lookupOrDefault<Switch>("oldFashion", "no")),
    relaxationFactor_(lookupOrDefault<scalar>("relaxationFactor", 1.0)),
    interfaceScan_(mesh_, mixture.alpha1()),
    condThresh_(lookupOrDefault<scalar>("condThresh", 0.5)),
    evapThresh_(lookupOrDefault<scalar>("evapThresh", 0.5)),
    interfaceField_
    (
        IOobject
        (
            "interfaceField",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    wallField_
    (
        IOobject
        (
            "wallField",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseChangeModels::interfaceEquilibrium::~interfaceEquilibrium()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::interfaceEquilibrium::mDotAlphal() const
{
    NotImplemented;

    return Pair<tmp<volScalarField>>
    (
        volScalarField::New
        (
            "mDotcAlphal",
            mesh_,
            dimensionedScalar(dimDensity/dimTime, Zero)
        ),
        volScalarField::New
        (
            "mDotvAlphal",
            mesh_,
            dimensionedScalar(dimDensity/dimTime, Zero)
        )
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::interfaceEquilibrium::mDot() const
{
    NotImplemented;

    return Pair<tmp<volScalarField>>
    (
        volScalarField::New
        (
            "mDotc",
            mesh_,
            dimensionedScalar(dimDensity/dimTime, Zero)
        ),
        volScalarField::New
        (
            "mDotv",
            mesh_,
            dimensionedScalar(dimDensity/dimTime, Zero)
        )
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::interfaceEquilibrium::mDotDeltaT() const
{
    NotImplemented;

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
Foam::phaseChangeModels::interfaceEquilibrium::TSource() const
{
    NotImplemented;

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
}


void Foam::phaseChangeModels::interfaceEquilibrium::correct()
{
    saturationTemperaturePtr_->correct();
    
    const dimensionedScalar deltaT = mesh_.time().deltaT();
    const dimensionedScalar rDeltaT = scalar(1)/mesh_.time().deltaT();

    const volScalarField& T = mixture_.T();
    const volScalarField& rho1 = mixture_.thermo1().rho();
    const volScalarField& rho2 = mixture_.thermo2().rho();
    const volScalarField& Lv = this->Lv();
    const volScalarField& TSat = this->TSat();

    // specific heat capacity of the mixture, not for single component!
    tmp<volScalarField> tCp = mixture_.Cp();
    const volScalarField& Cp = tCp();

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    if (!mesh_.objectRegistry::foundObject<volScalarField>("rho"))
    {
        FatalErrorIn
        (
            "interfaceEquilibrium::correct"
        )   << "mixture density rho not found in the objectRegistry"
            << exit(FatalError);
    }
    const volScalarField& rho = 
        mesh_.objectRegistry::lookupObjectRef<volScalarField>("rho");


    if (global_)
    {
        volScalarField S1 // unit [W/m^3]
        (
            rho*Cp*(TSat - T.oldTime())*rDeltaT
        );
        volScalarField S2C
        (
            (scalar(1) - limitedAlpha1)*rho2*Lv*rDeltaT
        );
        volScalarField S2V
        (
            -limitedAlpha1*rho1*Lv*rDeltaT
        );
        volScalarField S3C
        (
            Lv*rDeltaT/(scalar(1)/rho2 - scalar(1)/rho1)
        );
        volScalarField S3V
        (
            -Lv*rDeltaT/(scalar(1)/rho2 - scalar(1)/rho1)
        );

        mDot_ = pos(S1)*( min(min(S1, S2C), S3C)/Lv )
                + neg(S1)*( max(max(S1, S2V), S3V)/Lv );
        mDot_ *= relaxationFactor_;
    }
    else
    {
        List<interfaceMarker::cellPairInfo> condCellPairs;
        List<interfaceMarker::cellPairInfo> evapCellPairs;

        interfaceScan_.reset();
        interfaceScan_.getInterfaceCellPairs(condCellPairs, condThresh_);

        interfaceScan_.reset();
        interfaceScan_.getInterfaceCellPairs(evapCellPairs, evapThresh_);

        interfaceField_ = 0.0;
        wallField_ = 0.0;

        surfaceScalarField Tf = fvc::interpolate(T.oldTime());
        surfaceScalarField TSatf = fvc::interpolate(TSat);

        forAll(condCellPairs, i)
        {
            const label faceI = condCellPairs[i].faceIndex;

            // 如果 T <= TSat, 这标记为相界面
            if (Tf[faceI] <= TSatf[faceI])
            {
                interfaceField_[ condCellPairs[i].cell1 ] = 1.0;
                interfaceField_[ condCellPairs[i].cell2 ] = 1.0;
            }
        }

        forAll(evapCellPairs, i)
        {
            const label faceI = evapCellPairs[i].faceIndex;

            if (Tf[faceI] >= TSatf[faceI])
            {
                interfaceField_[ evapCellPairs[i].cell1 ] = 1.0;
                interfaceField_[ evapCellPairs[i].cell2 ] = 1.0;
            }
        }

        labelList wallCells;
        forAll(mesh_.boundaryMesh(), patchI)
        {
            if (isA<wallFvPatch>(mesh_.boundary()[patchI])
                    || isA<wallPolyPatch>(mesh_.boundaryMesh()[patchI]))
            {
                wallCells.append(mesh_.boundaryMesh()[patchI].faceCells());
            }
        }

        forAll(wallCells, i)
        {
            wallField_[wallCells[i]] = 1.0;

            // 同时把wall cells纳入interfaceField_中
            interfaceField_[wallCells[i]] = 1.0;
        }


        /* ---------------------------------------------------------------------- *\
                                    Old Fashion Road Begin
        \* ---------------------------------------------------------------------- */
        if (oldFashion_)
        {
            // Unlimited phase change mass transfer
            // > 0 for evaporation, < 0 for condensation
            volScalarField qDot =
                 interfaceField_*rho*Cp*(T.oldTime() - TSat)*rDeltaT;

            // Fluid availability limiting
            volScalarField limCond = (1.0 - limitedAlpha1)*rho2*Lv*rDeltaT; // > 0
            // No evaporation on wall
            // volScalarField limEvap = (1.0 - wallField_)*limitedAlpha1*rho1*Lv*rDeltaT; // > 0
            // Allowed evaporation on wall
            volScalarField limEvap = limitedAlpha1*rho1*Lv*rDeltaT;

            // Apply fluid limiting
            volScalarField qDot_fluid =
                neg(qDot)*max(qDot, -limCond) + pos(qDot)*min(qDot, limEvap);

            // Volume-based limiting (i.e. relative phase change rate can't exceed |1| per time step)
            volScalarField volumeFac = 
                (qDot/Lv)*( (1.0/rho2) - (1.0/rho1) )*deltaT;

            // Volume generation/sink based limiting
            // Again, regular evaporation on wall is not allowed
            // volScalarField qDot_volume =
            //     qDot*mag( min(max(1.0/(volumeFac + SMALL), -1.0), (1.0 - wallField_)) );
            // Allowed evaporation on wall
            volScalarField qDot_volume =
                qDot*mag( min(max(1.0/(volumeFac + SMALL), -1.0), 1.0) );

            // Final limiting
            qDot =
                   neg(qDot)*max(max(qDot, qDot_fluid), qDot_volume) 
                 + pos(qDot)*min(min(qDot, qDot_fluid), qDot_volume);

            // Under relax phase change rate per user specification
            mDot_ = -qDot/Lv;
            mDot_ *= relaxationFactor_;
        }
        /* ---------------------------------------------------------------------- *\
                                    Old Fashion Road End
        \* ---------------------------------------------------------------------- */
        else
        {
            // 不限制壁面处的蒸发和冷凝
            volScalarField S1 // 冷凝时大于0, 蒸发时小于0
            (
                interfaceField_*rho*Cp*(TSat - T.oldTime())*rDeltaT
            );

            volScalarField S2C // 大于0
            (
                (scalar(1) - limitedAlpha1)*rho2*Lv*rDeltaT
            );
            volScalarField S2V // 小于0
            (
                -limitedAlpha1*rho1*Lv*rDeltaT
            );

            volScalarField S3C // 大于0
            (
                Lv*rDeltaT/(scalar(1)/rho2 - scalar(1)/rho1)
            );
            volScalarField S3V // 小于0
            (
                -Lv*rDeltaT/(scalar(1)/rho2 - scalar(1)/rho1)
            );

            mDot_ = pos(S1)*( min(min(S1, S2C), S3C)/Lv )
                    + neg(S1)*( max(max(S1, S2V), S3V)/Lv );
                    
            mDot_ *= relaxationFactor_;
        }
    }
}


bool Foam::phaseChangeModels::interfaceEquilibrium::read
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

        phaseChangeDict.readIfPresent<Switch>
        (
            "global", global_
        );
        phaseChangeDict.readIfPresent<Switch>
        (
            "oldFashion", oldFashion_
        );
        phaseChangeDict.readIfPresent<scalar>
        (
            "relaxationFactor", relaxationFactor_
        );
        phaseChangeDict.readIfPresent<scalar>
        (
            "condThresh", condThresh_
        );
        phaseChangeDict.readIfPresent<scalar>
        (
            "evapThresh", evapThresh_
        );

        return true;
    }

    return false;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
