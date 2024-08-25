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

#include "limiterSharpInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "cutCellIso.H"
#include "volPointInterpolation.H"
#include "wallPolyPatch.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace phaseChangeModels
{
    defineTypeNameAndDebug(limiterSharpInterface, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeModel,
        limiterSharpInterface,
        dictionary
    );
}

}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::phaseChangeModels::limiterSharpInterface::restrainInter()
{

    volScalarField& T = const_cast<volScalarField&>(mixture_.T());
    const volScalarField& TSat = this->TSat();

    List<label> interCells(0);
    interfaceScan_.reset();

    // interfaceScan_.getDoubleLayerInterfaceCells(interCells, interVOF_);
    interfaceScan_.getInterfaceCells(interCells, interVOF_);

    forAll(interCells, i)
    {
        label celli = interCells[i];
        T[celli] = TSat[celli];
    }

    // forAll(interfaceArea_, celli)
    // {
    //     if (mag(interfaceArea_[celli]) > 1e-3)
    //     {
    //         T[celli] = TSat[celli];
    //     }
    // }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::limiterSharpInterface::limiterSharpInterface
(
    const twoPhaseMixtureThermo& mixture
)
:
    phaseChangeModel(mixture),
    interfaceArea_
    (
        IOobject
        (
            "interfaceArea",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector(0, 0, 0))
    ),
    interfaceAreaMethod_
    (
        lookupOrDefault<word>("interfaceAreaMethod", "gradAlpha")
    ),
    mDotc_
    (
        IOobject
        (
            "mDotc",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),
    mDotv_
    (
        IOobject
        (
            "mDotv",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),
    flux_
    (
        IOobject
        (
            "fluxCmpt",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature/dimArea, Zero)
    ),
    limiter_(lookupOrDefault<Switch>("limiter", "yes")),
    isExplicit_(lookupOrDefault<Switch>("isExplicit", "yes")),
    interfaceScan_(mesh_, mixture.alpha1()),
    interVOF_(lookupOrDefault<scalar>("interVOF", 0.5))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseChangeModels::limiterSharpInterface::~limiterSharpInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseChangeModels::limiterSharpInterface::updateInterface
(
    const volVectorField& normal
)
{
    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();

    if (interfaceAreaMethod_ == "gradAlpha")
    {
        volScalarField limitedAlpha1
        (
            min(max(mixture_.alpha1(), scalar(0)), scalar(1))
        );

        interfaceArea_ = fvc::grad(limitedAlpha1);
    }
    else if (interfaceAreaMethod_ == "isoAlpha")
    {
        scalar isoAlphaVal_( lookupOrDefault<scalar>("isoAlphaVal", 0.5) );

        // Interpolating alpha1 cell centre values to mesh points (vertices)
        scalarField ap
        (
            volPointInterpolation::New(mesh_).interpolate(mixture_.alpha1())
        );

        cutCellIso cutCell(mesh_, ap);

        if ((isoAlphaVal_ >= 1) || (isoAlphaVal_ <= 0))
        {
            FatalErrorIn
            (
                "limiterSharpInterface::updateInterface"
            )   << "Inavailable isoAlpha value " << isoAlphaVal_ << nl
                << "Valid value should be range in (0, 1)"
                << exit(FatalError);
        }

        forAll(interfaceArea_, celli)
        {
            label status = cutCell.calcSubCell(celli, isoAlphaVal_);
            interfaceArea_[celli] = vector(0, 0, 0);

            if (status == 0) // cell is cut
            {
                interfaceArea_[celli] =
                    cutCell.faceArea()/mesh_.V()[celli];
            }
        }
    }
    else if (interfaceAreaMethod_ == "isoAdvection")
    {
        forAll(interfaceArea_, celli)
        {
            interfaceArea_[celli] = 
                normal[celli]/mesh_.V()[celli];
        }
    }
    else
    {
        FatalErrorIn
        (
            "limiterSharpInterface::updateInterface"
        )   << "Unknown interface area method, available methods are:\n"
            << "gradAlpha, isoAlpha, isoAdvection"
            << exit(FatalError);
    }

    // 把近壁面单元格计入interfaceArea_中
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchi)
    {
        if (isA<wallPolyPatch>(pbm[patchi]))
        {
            const polyPatch& pp = pbm[patchi];
            forAll(pp.faceCells(), i)
            {
                const label pCelli = pp.faceCells()[i];

                // 仅当冷凝时计入, 蒸发时不计入
                if
                (
                    (TSat[pCelli] - T[pCelli]) > 0
                  && mixture_.alpha1()[pCelli] < 0.9
                )
                {
                    interfaceArea_[pCelli] =
                        pp.faceAreas()[i]/mesh_.V()[pCelli];
                }

                // // 不考虑饱和温度的影响
                // if (mixture_.alpha1()[pCelli] < 0.9)
                // {
                //     interfaceArea_[pCelli] =
                //         pp.faceAreas()[i]/mesh_.V()[pCelli];
                // }
            }
        }
    }
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::limiterSharpInterface::mDotAlphal() const
{
    // tmp<volScalarField> tkappa1 = mixture_.thermo1().kappa();
    // tmp<volScalarField> tkappa2 = mixture_.thermo2().kappa();
    // const volScalarField& kappa1 = tkappa1();
    // const volScalarField& kappa2 = tkappa2();

    // const volScalarField& T = mixture_.T();
    // const volScalarField& Lv = this->Lv();

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(mixture_.alpha2(), scalar(0)), scalar(1))
    );

    return Pair<tmp<volScalarField>>
    (
        mDotc_/(limitedAlpha2 + SMALL),
        mDotv_/(limitedAlpha1 + SMALL)
    );

    // NotImplemented;

    // return Pair<tmp<volScalarField>>
    // (
    //     volScalarField::New
    //     (
    //         "mDotcAlphal",
    //         mesh_,
    //         dimensionedScalar(dimDensity/dimTime, Zero)
    //     ),
    //     volScalarField::New
    //     (
    //         "mDotvAlphal",
    //         mesh_,
    //         dimensionedScalar(dimDensity/dimTime, Zero)
    //     )
    // );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::limiterSharpInterface::mDot() const
{
    return Pair<tmp<volScalarField>>
    (
        tmp<volScalarField>(mDotc_),
        tmp<volScalarField>(mDotv_)
    );

    // NotImplemented;

    // return Pair<tmp<volScalarField>>
    // (
    //     volScalarField::New
    //     (
    //         "mDotc",
    //         mesh_,
    //         dimensionedScalar(dimDensity/dimTime, Zero)
    //     ),
    //     volScalarField::New
    //     (
    //         "mDotv",
    //         mesh_,
    //         dimensionedScalar(dimDensity/dimTime, Zero)
    //     )
    // );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::limiterSharpInterface::mDotDeltaT() const
{
    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();

    dimensionedScalar SMALL_T("SMALL_T", dimTemperature, SMALL);

    return Pair<tmp<volScalarField>> // 有很大的除零风险
    (
        // mDotc_*pos(TSat - T)/max(1e-6*cmptAv(TSat), TSat - T),
        // mDotv_*pos(T - TSat)/min(-1e-6*cmptAv(TSat), TSat - T)
        // 参考interfaceHeatResistance
        mDotc_*pos(TSat - T)/(TSat - T + SMALL_T),
        mDotv_*pos(T - TSat)/(TSat - T - SMALL_T)
    );

    // NotImplemented;

    // return Pair<tmp<volScalarField>>
    // (
    //     volScalarField::New
    //     (
    //         "mDotcDeltaT",
    //         mesh_,
    //         dimensionedScalar(dimensionSet(1,-3,-1,-1,0,0,0), Zero)
    //     ),
    //     volScalarField::New
    //     (
    //         "mDotvDeltaT",
    //         mesh_,
    //         dimensionedScalar(dimensionSet(1,-3,-1,-1,0,0,0), Zero)
    //     )
    // );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::phaseChangeModels::limiterSharpInterface::TSource() const
{
    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();
    const volScalarField& Lv = this->Lv();

    tmp<fvScalarMatrix> tTSource
    (
        new fvScalarMatrix
        (
            T,
            dimEnergy/dimTime
        )
    );

    fvScalarMatrix& TSource = tTSource.ref();

    // Pair<tmp<volScalarField>> mDotDeltaT = 
    //     this->mDotDeltaT();

    // const volScalarField& mDotcDeltaT = mDotDeltaT[0]();
    // const volScalarField& mDotvDeltaT = mDotDeltaT[1]();

    // // interface heat resistance
    // const volScalarField IHRCoeff = 
    //     (mDotcDeltaT + mDotvDeltaT)*Lv;

    // // condensation negetive, evaporation positive
    // TSource = fvm::Sp(IHRCoeff, T) - IHRCoeff*TSat;

    return tTSource;

    // NotImplemented;

    // const volScalarField& T = mixture_.T();

    // tmp<fvScalarMatrix> tTSource
    // (
    //     new fvScalarMatrix
    //     (
    //         T,
    //         dimEnergy/dimTime // 焦耳/秒
    //     )
    // );

    // return tTSource;
}


void Foam::phaseChangeModels::limiterSharpInterface::correct()
{
    saturationTemperaturePtr_->correct();
    restrainInter();

    tmp<volScalarField> tkappa1 = mixture_.thermo1().kappa();
    tmp<volScalarField> tkappa2 = mixture_.thermo2().kappa();
    const volScalarField& kappa1 = tkappa1();
    const volScalarField& kappa2 = tkappa2();
    const volScalarField& rho1 = mixture_.thermo1().rho();
    const volScalarField& rho2 = mixture_.thermo2().rho();

    const volScalarField& T = mixture_.T();
    const volScalarField& Lv = this->Lv();

    // flux_ < 0 for evaporation, flux_ > 0 for condensation
    // vapor is unsaturated phase for evaporation
    // liquid is unsaturated phase for condensation
    flux_ = fvc::grad(T) & interfaceArea_;

    mDotc_ = pos(flux_)*flux_*(2.0*kappa1)/Lv; // >= 0
    mDotv_ = neg(flux_)*flux_*(2.0*kappa2)/Lv; // <= 0

    // mDotc_ = pos(TSat - T)*pos(flux_)*flux_*(2.0*kappa1)/Lv; // >= 0
    // mDotv_ = pos(T - TSat)*neg(flux_)*flux_*(2.0*kappa2)/Lv; // <= 0

    if (limiter_)
    {
        forAll(mDotc_, celli) // 限制器
        {
            scalar maxCond = 
                mixture_.alpha2()[celli]*rho2[celli]
                /mesh_.time().deltaTValue();
            scalar maxEvap = 
                mixture_.alpha1()[celli]*rho1[celli]
                /mesh_.time().deltaTValue();
            
            mDotc_[celli] = min(mDotc_[celli], maxCond);
            mDotv_[celli] = max(mDotv_[celli], -maxEvap);
        }
    }

    mDot_ = mDotc_ + mDotv_;
}


Foam::tmp<Foam::volScalarField>
Foam::phaseChangeModels::limiterSharpInterface::Sq() const
{
    // tmp<volScalarField> tCp1 = mixture_.thermo1().Cp();
    // tmp<volScalarField> tCp2 = mixture_.thermo2().Cp();
    // const volScalarField& Cp1 = tCp1();
    // const volScalarField& Cp2 = tCp2();

    // const volScalarField& TSat = this->TSat();
    // const volScalarField& T  = mixture_.T();
    // const volScalarField& Lv = this->Lv();

    // // explicit evaluation
    // return (Lv + (Cp1 - Cp2)*(TSat - T))*mDot_;

    return tmp<volScalarField>
    (
        volScalarField::New
        (
            "Sq",
            mesh_,
            dimensionedScalar(dimPower/dimVolume, Zero)
        )
    );
}


bool Foam::phaseChangeModels::limiterSharpInterface::read
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

        phaseChangeDict.readIfPresent<word>
        (
            "interfaceAreaMethod", interfaceAreaMethod_
        );

        phaseChangeDict.readIfPresent<Switch>
        (
            "limiter", limiter_
        );

        phaseChangeDict.readIfPresent<Switch>
        (
            "isExplicit", isExplicit_
        );

        phaseChangeDict.readIfPresent<scalar>
        (
            "interVOF", interVOF_
        );

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
