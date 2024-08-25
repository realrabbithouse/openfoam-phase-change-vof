IOobject alphaPhi10Header
(
    IOobject::groupName("alphaPhi0", alpha1.group()),
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

const bool alphaRestart =
    alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true);

if (alphaRestart)
{
    Info << "Restarting alpha" << endl;
}

// MULES flux from previous time-step
surfaceScalarField alphaPhi10
(
    alphaPhi10Header,
    phi*fvc::interpolate(alpha1)
);
/*
alphaPhi10的作用:
1) 构造turbulence model;
2) 当turbulence model的twoPhaseTransport_为true时, 构造各相的mass flux.
*/

// MULES Correction
tmp<surfaceScalarField> talphaPhi1Corr0;
