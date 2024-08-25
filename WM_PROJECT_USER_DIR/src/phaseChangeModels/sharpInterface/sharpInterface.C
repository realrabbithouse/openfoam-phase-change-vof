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

#include "sharpInterface.H"
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
    defineTypeNameAndDebug(sharpInterface, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeModel,
        sharpInterface,
        dictionary
    );
}

}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::sharpInterface::sharpInterface(const twoPhaseMixtureThermo& mixture)
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
    isExplicit_(lookupOrDefault<Switch>("isExplicit", "yes"))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseChangeModels::sharpInterface::~sharpInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseChangeModels::sharpInterface::updateInterface
(
    const volVectorField& normal
)
{
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

        // interface heat resistance
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
                "sharpInterface::updateInterface"
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
            "sharpInterface::updateInterface"
        )   << "Unknown interface area method, available methods are:\n"
            << "gradAlpha, isoAlpha, isoAdvection"
            << exit(FatalError);
    }

    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();

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
                    (TSat[pCelli] - T.oldTime()[pCelli]) > 0
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
Foam::phaseChangeModels::sharpInterface::mDotAlphal() const
{
    tmp<volScalarField> tkappa1 = mixture_.thermo1().kappa();
    tmp<volScalarField> tkappa2 = mixture_.thermo2().kappa();
    const volScalarField& kappa1 = tkappa1();
    const volScalarField& kappa2 = tkappa2();

    const volScalarField& Lv = this->Lv();

    return Pair<tmp<volScalarField>>
    (
        pos(flux_)*flux_*(2.0*kappa1)/Lv,
        neg(flux_)*flux_*(2.0*kappa2)/Lv
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
Foam::phaseChangeModels::sharpInterface::mDot() const
{
    tmp<volScalarField> tkappa1 = mixture_.thermo1().kappa();
    tmp<volScalarField> tkappa2 = mixture_.thermo2().kappa();
    const volScalarField& kappa1 = tkappa1();
    const volScalarField& kappa2 = tkappa2();

    const volScalarField& Lv = this->Lv();

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    return Pair<tmp<volScalarField>>
    (
        pos(flux_)*(scalar(1) - limitedAlpha1)*flux_*(2.0*kappa1)/Lv,
        neg(flux_)*limitedAlpha1*flux_*(2.0*kappa2)/Lv
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
Foam::phaseChangeModels::sharpInterface::mDotDeltaT() const
{
    tmp<volScalarField> tkappa1 = mixture_.thermo1().kappa();
    tmp<volScalarField> tkappa2 = mixture_.thermo2().kappa();
    const volScalarField& kappa1 = tkappa1();
    const volScalarField& kappa2 = tkappa2();

    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();
    const volScalarField& Lv = this->Lv();

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    dimensionedScalar SMALL_T("SMALL_T", dimTemperature, SMALL);

    return Pair<tmp<volScalarField>> // 有很大的除零风险
    (
        // pos(flux_)*(scalar(1) - limitedAlpha1)*flux_*(2.0*kappa1)
        // /( Lv*max(1e-6*cmptAv(TSat), TSat - T.oldTime()) ),
        // neg(flux_)*limitedAlpha1*flux_*(2.0*kappa2)
        // /( Lv*min(-1e-6*cmptAv(TSat), TSat - T.oldTime()) )
        pos(flux_)*flux_*(scalar(1) - limitedAlpha1)*(2.0*kappa1)
        *pos(TSat - T.oldTime())/( Lv*(TSat - T.oldTime() + SMALL_T) ),
        neg(flux_)*flux_*limitedAlpha1*(2.0*kappa2)
        *pos(T.oldTime() - TSat)/( Lv*(TSat - T.oldTime() - SMALL_T) )
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
Foam::phaseChangeModels::sharpInterface::TSource() const
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

    Pair<tmp<volScalarField>> mDotDeltaT = 
        this->mDotDeltaT();

    const volScalarField& mDotcDeltaT = mDotDeltaT[0]();
    const volScalarField& mDotvDeltaT = mDotDeltaT[1]();

    // interface heat resistance
    const volScalarField IHRCoeff = 
        (mDotcDeltaT + mDotvDeltaT)*Lv;

    // condensation negetive, evaporation positive
    TSource = fvm::Sp(IHRCoeff, T) - IHRCoeff*TSat;

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


void Foam::phaseChangeModels::sharpInterface::correct()
{
    saturationTemperaturePtr_->correct();

    tmp<volScalarField> tkappa1 = mixture_.thermo1().kappa();
    tmp<volScalarField> tkappa2 = mixture_.thermo2().kappa();
    const volScalarField& kappa1 = tkappa1();
    const volScalarField& kappa2 = tkappa2();

    const volScalarField& T = mixture_.T();
    const volScalarField& Lv = this->Lv();

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    // flux_ < 0 for evaporation, flux_ > 0 for condensation
    // vapor is unsaturated phase for evaporation
    // liquid is unsaturated phase for condensation
    flux_ = fvc::grad(T.oldTime()) & interfaceArea_;

    // correct() 之前必须先 updateInterface(normal)
    mDot_ = pos(flux_)*(scalar(1) - limitedAlpha1)*flux_*(2.0*kappa1)/Lv
            + neg(flux_)*limitedAlpha1*flux_*(2.0*kappa2)/Lv;
}


Foam::tmp<Foam::volScalarField>
Foam::phaseChangeModels::sharpInterface::Sq() const
{
    tmp<volScalarField> tCp1 = mixture_.thermo1().Cp();
    tmp<volScalarField> tCp2 = mixture_.thermo2().Cp();
    const volScalarField& Cp1 = tCp1();
    const volScalarField& Cp2 = tCp2();

    const volScalarField& TSat = this->TSat();
    const volScalarField& T  = mixture_.T();
    const volScalarField& Lv = this->Lv();

    // explicit evaluation
    return (Lv + (Cp1 - Cp2)*(TSat - T.oldTime()))*mDot_;
}


bool Foam::phaseChangeModels::sharpInterface::read
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
            "isExplicit", isExplicit_
        );

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
