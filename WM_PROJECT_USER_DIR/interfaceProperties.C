/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "interfaceProperties.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad() * acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::interfaceProperties::calculateK()
{
    if (smoothedCurvature_)
    {
        /*
        ###需要回答的问题###
        1. smoothed alpha1是否会影响alpha1梯度和nHatf_的计算?
        2. smothed alpha1是否会影响correctContactAngle(...)?
        */
        const fvMesh& mesh = alpha1_.mesh();
        const surfaceVectorField& Sf = mesh.Sf();
        const surfaceScalarField& magSf = mesh.magSf();

        smoothedAlpha1_ = alpha1_;
        for (label i = 0; i < smoothedCycles_; ++i)
        {
            smoothedAlpha1_ = fvc::surfaceSum(fvc::interpolate(smoothedAlpha1_)*magSf)
                /fvc::surfaceSum(magSf);
        }

        // const surfaceVectorField gradSAlphaf
        // (
        //     IOobject
        //     (
        //         "gradSAlphaf",
        //         alpha1_.time().timeName(),
        //         alpha1_.mesh()
        //     ),
        //     fvc::interpolate(fvc::grad(smoothedAlpha1_))
        // );

        // Cell gradient of alpha
        const volVectorField gradSAlpha(fvc::grad(smoothedAlpha1_, "nHat"));

        // Interpolated face-gradient of alpha
        surfaceVectorField gradSAlphaf(fvc::interpolate(gradSAlpha));

        surfaceVectorField nHatSfv(gradSAlphaf/(mag(gradSAlphaf) + deltaN_));

        // correct dynamic contact angle based on the smoothed alpha1
        correctContactAngle(nHatSfv.boundaryFieldRef(), gradSAlphaf.boundaryField());

        nHatf_ = nHatSfv & Sf;
        K_ = -fvc::div(nHatf_);
    }
    else
    {
        const fvMesh& mesh = alpha1_.mesh();
        const surfaceVectorField& Sf = mesh.Sf();

        // Cell gradient of alpha
        const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));

        // Interpolated face-gradient of alpha
        surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

        // gradAlphaf -=
        //    (mesh.Sf()/mesh.magSf())
        //   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

        // Face unit interface normal
        surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
        // surfaceVectorField nHatfv
        // (
        //     (gradAlphaf + deltaN_*vector(0, 0, 1)
        //    *sign(gradAlphaf.component(vector::Z)))/(mag(gradAlphaf) + deltaN_)
        // );
        correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField());

        // Face unit interface normal flux
        nHatf_ = nHatfv & Sf;

        // Simple expression for curvature
        K_ = -fvc::div(nHatf_);

        // Complex expression for curvature.
        // Correction is formally zero but numerically non-zero.
        /*
        volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
        forAll(nHat.boundaryField(), patchi)
        {
            nHat.boundaryFieldRef()[patchi] = nHatfv.boundaryField()[patchi];
        }

        K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
        */
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        alpha1.mesh().solverDict(alpha1.name()).get<scalar>("cAlpha")
    ),

    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),

    deltaN_
    (
        "deltaN",
        1e-8/cbrt(average(alpha1.mesh().V()))
    ),

    alpha1_(alpha1),
    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimArea, Zero)
    ),

    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    ),

    smoothedCurvature_
    (
        alpha1.mesh().solverDict(alpha1.name()).getOrDefault<Switch>("smoothedCurvature", "no")
    ),
    smoothedCycles_
    (
        alpha1.mesh().solverDict(alpha1.name()).getOrDefault<label>("smoothedCycles", 2)
    ),
    smoothedAlpha1_
    (
        IOobject
        (
            "smoothedAlpha1",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_
    )
{
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::sigmaK() const
{
    return sigmaPtr_->sigma()*K_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::surfaceTensionForce() const
{
    tmp<surfaceScalarField> sigmaKf = fvc::interpolate(sigmaK());
    return  sigmaKf*fvc::snGrad(smoothedAlpha1_);
    
    // return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_);
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}


void Foam::interfaceProperties::correct()
{
    calculateK();
}


bool Foam::interfaceProperties::read()
{
    alpha1_.mesh().solverDict(alpha1_.name()).readEntry("cAlpha", cAlpha_);
    sigmaPtr_->readDict(transportPropertiesDict_);

    return true;
}


// ************************************************************************* //
