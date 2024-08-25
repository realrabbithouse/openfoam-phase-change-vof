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

#include "Schrage.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "thermodynamicConstants.H"
#include "cutCellIso.H"
#include "volPointInterpolation.H"
#include "wallPolyPatch.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace phaseChangeModels
{
    defineTypeNameAndDebug(Schrage, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeModel,
        Schrage,
        dictionary
    );
}

}

const Foam::dimensionedScalar Foam::phaseChangeModels::Schrage::R0
(
    "R0", dimensionSet(1,2,-2,-1,-1,0,0), Foam::constant::thermodynamic::RR
);


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeModels::Schrage::phaseChangeCoeff() const
{
    // molecular weight of the gas phase
    tmp<volScalarField> tW = mixture_.thermo2().W();
    const volScalarField& W = tW();

    const volScalarField& Lv = this->Lv();
    const volScalarField& TSat = this->TSat();
    const volScalarField& rho2 = mixture_.thermo2().rho();

    // dimensions = [1, -2, -1, -1, 0]
    return (2.0*gamma_/(2.0 - gamma_))*sqrt(W/(2.0*pi_*R0))
                *rho2*Lv/sqrt(pow3(TSat));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModels::Schrage::Schrage(const twoPhaseMixtureThermo& mixture)
:
    phaseChangeModel(mixture),
    gamma_(lookupOrDefault<scalar>("gamma", 1.0)),
    pi_(Foam::constant::mathematical::pi),
    interfaceAreaMethod_
    (
        lookupOrDefault<word>("interfaceAreaMethod", "gradAlpha")
    ),
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
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    phaseChangeCoeff_
    (
        IOobject
        (
            "phaseChangeCoeff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimensionSet(1,-2,-1,-1,0,0,0), Zero)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseChangeModels::Schrage::~Schrage()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseChangeModels::Schrage::updateInterface
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

        interfaceArea_ = mag(fvc::grad(limitedAlpha1));
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
                "Schrage::updateInterface"
            )   << "Inavailable isoAlpha value " << isoAlphaVal_ << nl
                << "Valid value should be range in (0, 1)"
                << exit(FatalError);
        }

        forAll(interfaceArea_, celli)
        {
            label status = cutCell.calcSubCell(celli, isoAlphaVal_);
            interfaceArea_[celli] = 0;
            if (status == 0) // cell is cut
            {
                interfaceArea_[celli] =
                    mag(cutCell.faceArea())/mesh_.V()[celli];
            }
        }
    }
    else if (interfaceAreaMethod_ == "isoAdvection")
    {
        forAll(interfaceArea_, celli)
        {
            interfaceArea_[celli] = 
                mag(normal[celli])/mesh_.V()[celli];
        }
    }
    else
    {
        FatalErrorIn
        (
            "Schrage::updateInterface"
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
                        mag(pp.faceAreas()[i])/mesh_.V()[pCelli];
                }
            }
        }
    }
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::Schrage::mDotAlphal() const
{
    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();
    const dimensionedScalar T0(dimTemperature, Zero);

    // volScalarField mDotcAlphal
    // (
    //     "mDotcAlphal",
    //     phaseChangeCoeff_*interfaceArea_*max(TSat - T.oldTime(), T0)
    // );

    // volScalarField mDotvAlphal
    // (
    //     "mDotvAlphal",
    //     phaseChangeCoeff_*interfaceArea_*min(TSat - T.oldTime(), T0)
    // );

    return Pair<tmp<volScalarField>>
    (
        phaseChangeCoeff_*interfaceArea_*max(TSat - T.oldTime(), T0),
        phaseChangeCoeff_*interfaceArea_*min(TSat - T.oldTime(), T0)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::Schrage::mDot() const
{
    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();
    const dimensionedScalar T0(dimTemperature, Zero);

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    // volScalarField mDotc
    // (
    //     "mDotc",
    //     phaseChangeCoeff_*interfaceArea_*max(TSat - T.oldTime(), T0)
    //     *(scalar(1) - limitedAlpha1)
    // );

    // volScalarField mDotv
    // (
    //     "mDotv",
    //     phaseChangeCoeff_*interfaceArea_*min(TSat - T.oldTime(), T0)
    //     *limitedAlpha1
    // );

    return Pair<tmp<volScalarField>>
    (
        phaseChangeCoeff_*interfaceArea_
        *max(TSat - T.oldTime(), T0)*(scalar(1) - limitedAlpha1),
        phaseChangeCoeff_*interfaceArea_
        *min(TSat - T.oldTime(), T0)*limitedAlpha1
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeModels::Schrage::mDotDeltaT() const
{
    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    // volScalarField mDotcDeltaT
    // (
    //     "mDotcDeltaT",
    //     phaseChangeCoeff_*interfaceArea_*pos(TSat - T.oldTime())
    //     *(scalar(1) - limitedAlpha1)
    // );

    // volScalarField mDotvDeltaT
    // (
    //     "mDotvDeltaT",
    //     phaseChangeCoeff_*interfaceArea_*pos(T.oldTime() - TSat)
    //     *limitedAlpha1
    // );

    return Pair<tmp<volScalarField>>
    (
        phaseChangeCoeff_*interfaceArea_
        *(scalar(1) - limitedAlpha1)*pos(TSat - T.oldTime()),
        phaseChangeCoeff_*interfaceArea_
        *limitedAlpha1*pos(T.oldTime() - TSat)
    );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::phaseChangeModels::Schrage::TSource() const
{
    const volScalarField& T = mixture_.T();
    const volScalarField& TSat = this->TSat();
    const volScalarField& Lv = this-> Lv();

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
}


void Foam::phaseChangeModels::Schrage::correct()
{
    saturationTemperaturePtr_->correct();

    // molecular weight of the gas phase
    tmp<volScalarField> tW = mixture_.thermo2().W();
    const volScalarField& W = tW();

    const volScalarField& Lv = this->Lv();
    const volScalarField& TSat = this->TSat();
    const volScalarField& rho1 = mixture_.thermo1().rho();
    const volScalarField& rho2 = mixture_.thermo2().rho();

    // dimensions = [1, -2, -1, -1, 0]
    phaseChangeCoeff_ =
        (2.0*gamma_/(2.0 - gamma_))*sqrt(W/(2.0*pi_*R0))
        *(rho1*rho2/(rho1 - rho2))*Lv/sqrt(pow3(TSat));

    const volScalarField& T = mixture_.T();
    const dimensionedScalar T0(dimTemperature, Zero);

    volScalarField limitedAlpha1
    (
        min(max(mixture_.alpha1(), scalar(0)), scalar(1))
    );

    // correct() 之前必须先 updateInterface(normal)
    mDot_ = phaseChangeCoeff_*interfaceArea_*max(TSat - T.oldTime(), T0)
            *(scalar(1) - limitedAlpha1)
          + phaseChangeCoeff_*interfaceArea_*min(TSat - T.oldTime(), T0)
            *limitedAlpha1;
}


bool Foam::phaseChangeModels::Schrage::read
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

        phaseChangeDict.readIfPresent<scalar>
        (
            "gamma", gamma_
        );

        phaseChangeDict.readIfPresent<word>
        (
            "interfaceAreaMethod", interfaceAreaMethod_
        );

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
