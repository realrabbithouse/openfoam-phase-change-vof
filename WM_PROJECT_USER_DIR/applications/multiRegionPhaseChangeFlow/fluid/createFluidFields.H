// Initialise fluid field pointer lists
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceVectorField> UfFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<twoPhaseMixtureThermo> mixtureFluid(fluidRegions.size());
PtrList<phaseChangeModel> phaseChangeFluid(fluidRegions.size());

// alpha1, alpha2, rho1, rho2

PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<dimensionedScalar> pMinFluid(fluidRegions.size());

// set flux required for p_rgh and alpha1

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());

PtrList<surfaceScalarField> rhoPhiFluid(fluidRegions.size());
PtrList<volScalarField> dgdtFluid(fluidRegions.size());
PtrList<volScalarField> rhoCvFluid(fluidRegions.size());

// createAlphaFluxes.H
PtrList<surfaceScalarField> alphaPhi10Fluid(fluidRegions.size());

PtrList<isoAdvection> advectorFluid(fluidRegions.size());
PtrList<compressibleInterPhaseTransportModel> turbulenceFluid(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

List<bool> frozenFlowFluid(fluidRegions.size(), false);


forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to UfFluid\n" << endl;
    UfFluid.set
    (
        i,
        new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(UFluid[i])
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::flux(UFluid[i])
        )
    );

    Info<< "    Adding to mixtureFluid\n" << endl;
    mixtureFluid.set
    (
        i,
        new twoPhaseMixtureThermo(UFluid[i], phiFluid[i])
    );

    Info<< "    Adding to phaseChangeFluid\n" << endl;
    phaseChangeFluid.set
    (
        i,
        phaseChangeModel::New(mixtureFluid[i]).ptr()
    );

    // 当离开作用域时, 这几个名字会消失
    volScalarField& alpha1 = mixtureFluid[i].alpha1();
    volScalarField& alpha2 = mixtureFluid[i].alpha2();
    const volScalarField& rho1 = mixtureFluid[i].thermo1().rho();
    const volScalarField& rho2 = mixtureFluid[i].thermo2().rho();


    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            alpha1*rho1 + alpha2*rho2
        )
    );
    // rhoFluid[i].oldTime();

    Info<< "    Adding to pMinFluid\n" << endl;
    pMinFluid.set
    (
        i,
        new dimensionedScalar
        (
            "pMin",
            dimPressure,
            mixtureFluid[i]
        )
    );

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());
    fluidRegions[i].setFluxRequired(mixtureFluid[i].alpha1().name());


    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, Zero)
        )
    );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, Zero)
    );

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );


    Info<< "    Adding to rhoPhiFluid\n" << endl;
    rhoPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "rhoPhi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            // initialization doesn't matter
            fvc::interpolate(rhoFluid[i])*phiFluid[i]
        )
    );

    Info<< "    Adding to dgdtFluid\n" << endl;
    dgdtFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dgdt",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            alpha1*fvc::div(phiFluid[i])
        )
    );

    Info<< "    Adding to rhoCvFluid\n" << endl;
    rhoCvFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rhoCv",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            alpha1*rho1*mixtureFluid[i].thermo1().Cv()
            + alpha2*rho2*mixtureFluid[i].thermo2().Cv()
        )
    );
    // rhoCvFluid[i].oldTime();


    Info<< "    Adding to alphaPhi10Fluid\n" << endl;
    alphaPhi10Fluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("alphaPhi0", alpha1.group()),
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(alpha1)*phiFluid[i]
        )
    );


    Info<< "    Adding to advectorFluid\n" << endl;
    advectorFluid.set
    (
        i,
        new isoAdvection
        (
            alpha1,
            phiFluid[i],
            UFluid[i]
        )
    );

    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        new compressibleInterPhaseTransportModel
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            rhoPhiFluid[i],
            alphaPhi10Fluid[i],
            mixtureFluid[i]
        )
    );

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UFluid[i])
        )
    );

    Info<< "    Adding to MRFFluid\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );
}

