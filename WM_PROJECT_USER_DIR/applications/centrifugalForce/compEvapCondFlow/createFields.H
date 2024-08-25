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

Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
twoPhaseMixtureThermo mixture(U, phi);

Info<< "Constructing phaseChangeModel\n" << endl;
autoPtr<phaseChangeModel> phaseChange
(
    phaseChangeModel::New(mixture)
);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<< "Reading thermophysical properties\n" << endl;

const volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& rho2 = mixture.thermo2().rho();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
// rho.oldTime();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    mixture
);

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
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

volScalarField dgdt
(
    IOobject
    (
        "dgdt",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*fvc::div(phi)
);

#include "createAlphaFluxes.H"

isoAdvection advector(alpha1, phi, U);

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

#include "createK.H"

#include "createMRF.H"
#include "createFvOptions.H"

// Need to store rhoCv for ddt(rhoCv,T)
volScalarField rhoCv
(
    IOobject
    (
        "rhoCv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*rho1*mixture.thermo1().Cv() // 初始化不重要！
    + alpha2*rho2*mixture.thermo2().Cv()
);
// rhoCv.oldTime();

dimensionedVector omega("omega", dimless/dimTime, mixture);
dimensionedVector origin("origin", dimLength, mixture);

// // 测试相变热源
// volScalarField qDot_
// (
//     IOobject
//     (
//         "qDot_",
//         mesh.time().timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar(dimPower/dimVolume, Zero)
// );

// // 另一种方式测试相变热源
// volScalarField Sq_
// (
//     IOobject
//     (
//         "Sq_",
//         mesh.time().timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar(dimPower/dimVolume, Zero)
// );