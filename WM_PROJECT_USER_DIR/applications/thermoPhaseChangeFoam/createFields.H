#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


// Info<< "Reading transportProperties\n" << endl;
// immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);
Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
twoPhaseMixtureThermo mixture(U, phi);

Info<< "Constructing phaseChangeModel\n" << endl;
autoPtr<phaseChangeModel> phaseChange
(
    phaseChangeModel::New(mixture)
);

if (!mixture.incompressible())
{
    FatalErrorIn("twoPhaseMixtureThermo mixture(U, phi)")
        << "When constructing twoPhaseMixtureThermo, one of the components "
        << "is not incompressible" << nl
        << exit(FatalError);
}

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<< "Reading thermophysical properties\n" << endl;
const volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& rho2 = mixture.thermo2().rho();

Info<< "Reading field p\n" << endl;
volScalarField& p = mixture.p();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Need to store rho for ddt(rhoCp, T)
volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*rho1*mixture.thermo1().Cp() // 初始化不重要！
    + alpha2*rho2*mixture.thermo2().Cp()
);
rhoCp.oldTime();


#include "createAlphaFluxes.H"
// Construct compressible turbulence model
compressibleInterPhaseTransportModel turbulence
(
    rho,
    U,
    phi,
    rhoPhi,
    alphaPhi10, // 湍流模型需要alphaPhi10构造！！！
    mixture
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

volVectorField normal
(
    IOobject
    (
        "normal",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector(dimless, vector(0, 0, 0))
);

#include "createMRF.H"
#include "createFvOptions.H"
