{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    volScalarField::Internal divU
    (
        mesh.moving()
      ? fvc::div(phi + mesh.phi())
      : fvc::div(phi)
    );

    zeroField Su;
    zeroField Sp;

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    {
        tmp<fv::ddtScheme<scalar>> tddtAlpha
        (
            fv::ddtScheme<scalar>::New
            (
                mesh,
                mesh.ddtScheme("ddt(alpha)")
            )
        );
        const fv::ddtScheme<scalar>& ddtAlpha = tddtAlpha();

        if
        (
            isType<fv::EulerDdtScheme<scalar>>(ddtAlpha)
         || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha)
        )
        {
            ocCoeff = 0;
        }
        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha))
        {
            if (nAlphaSubCycles > 1)
            {
                FatalErrorInFunction
                    << "Sub-cycling is not supported "
                       "with the CrankNicolson ddt scheme"
                    << exit(FatalError);
            }

            if
            (
                alphaRestart
             || mesh.time().timeIndex() > mesh.time().startTimeIndex() + 1
            )
            {
                ocCoeff =
                    refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha)
                   .ocCoeff();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Only Euler and CrankNicolson ddt schemes are supported"
                << exit(FatalError);
        }
    }

    // Set the time blending factor, 1 for Euler
    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    // Standard face-flux compression coefficient
    surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    // Add the optional shear compression contribution
    if (scAlpha > 0)
    {
        phic +=
            scAlpha*mag(mesh.delta() & fvc::interpolate(symm(fvc::grad(U))));
    }


    surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime();
    }


    // Phase change sources
    volScalarField alphalDotSu
    (
        IOobject
        (
            "alphalDotSu",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless/dimTime, Zero)
    );
    volScalarField alphalDotSp
    (
        IOobject
        (
            "alphalDotSp",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless/dimTime, Zero)
    );

    if (phaseChange->isExplicit())
    {   
        // Implicit
        alphalDotSu = phaseChange->alphalSuSp()[0];
        alphalDotSp = phaseChange->alphalSuSp()[1];

        // Explicit
        // alphalDotSu = phaseChange->Salphal();
    }
    else
    {
        Pair<tmp<volScalarField>> vDotAlphal = 
            phaseChange->vDotAlphal();
        alphalDotSu = vDotAlphal[0](); // Su, > 0
        // const volScalarField& vDotvAlphal = vDotAlphal[1](); // < 0
        alphalDotSp = vDotAlphal[1]() - vDotAlphal[0](); // Sp
    }


    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
              : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiCN,
                upwind<scalar>(mesh, phiCN)
            ).fvmDiv(phiCN, alpha1)
       // - fvm::Sp(fvc::ddt(dimensionedScalar("1", dimless, 1), mesh)
       //           + fvc::div(phiCN), alpha1)
         ==
            alphalDotSu() + Su
            + fvm::Sp(divU + alphalDotSp() + Sp, alpha1)
        );

        alpha1Eqn.solve(); // Unbounded alpha1

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        tmp<surfaceScalarField> talphaPhi1UD(alpha1Eqn.flux());
        alphaPhi10 = talphaPhi1UD(); // 低阶通量

        if (alphaApplyPrevCorr && talphaPhi1Corr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct
            (
                geometricOneField(),
                alpha1,
                alphaPhi10,
                talphaPhi1Corr0.ref(),
                oneField(),
                zeroField()
            );

            alphaPhi10 += talphaPhi1Corr0();
        }

        // Cache the upwind-flux
        talphaPhi1Corr0 = talphaPhi1UD;

        // alpha1 = min(max(alpha1, scalar(0)), scalar(1)); // 使alpha1介与0到1之间

        alpha2 = 1.0 - alpha1; // Unbounded alpha2

        mixture.correct();
    }

    volScalarField alpha10("alpha10", alpha1);

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phir(phic*mixture.nHatf());

        tmp<surfaceScalarField> talphaPhi1Un
        (
            fvc::flux
            (
                phiCN(),
                cnCoeff*alpha1 + (1.0 - cnCoeff)*alpha1.oldTime(),
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhi1Corr(talphaPhi1Un() - alphaPhi10);
            volScalarField alpha10("alpha10", alpha1);

            // alpha100是alpha10上一次迭代的结果, 起初, alpha100 = alpha10 = alpha1 (Unbounded)
            volScalarField alpha100("alpha100", alpha10);

            alpha10 = alpha1; // alpha10是上一次MULESCorr更新的alpha1

            volScalarField::Internal SpCorr =
            alphalDotSp() + Sp;

            volScalarField::Internal SuCorr = 
            divU*(alpha10() - alpha100()) - (alphalDotSp() + Sp)*alpha10();
            // volScalarField::Internal SuCorr = -(alphalDotSp() + Sp)*alpha10();

            MULES::correct
            (
                geometricOneField(),
                alpha1, // non-const, initial unbounded alpha1 from solving alpha1Eqn
                talphaPhi1Un(), // uncorrected higher-order flux
                talphaPhi1Corr.ref(), // non-const, return corrected anti-diffusive flux
                SpCorr,
                SuCorr,
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi10 += talphaPhi1Corr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                alphaPhi10 += 0.5*talphaPhi1Corr();
            }
        }
        else
        {
            alphaPhi10 = talphaPhi1Un; // uncorrected higher-order flux

            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1, // non-const
                phiCN, // 用于构造低阶通量
                alphaPhi10, // non-const, uncorrected higher-order flux -> bounded higher-order flux
                alphalDotSp(), // (Sp + vDotvmcAlphal())()
                (alphalDotSu() + divU*min(alpha1(), scalar(1)))(),
                oneField(),
                zeroField()
            );
        }

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        // So-called anti-diffusive flux from previous time step
        talphaPhi1Corr0 = alphaPhi10 - talphaPhi1Corr0;
        talphaPhi1Corr0.ref().rename("alphaPhi1Corr0");
    }
    else
    {
        talphaPhi1Corr0.clear();
    }

    surfaceScalarField rho1f(fvc::interpolate(rho1));
    surfaceScalarField rho2f(fvc::interpolate(rho2));

    if
    (
        word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::EulerDdtScheme<vector>::typeName
     || word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::localEulerDdtScheme<vector>::typeName
    )
    {
        rhoPhi = alphaPhi10*(rho1f - rho2f) + phiCN*rho2f;
    }
    else
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step alpha flux
            alphaPhi10 =
                (alphaPhi10 - (1.0 - cnCoeff)*alphaPhi10.oldTime())/cnCoeff;
        }

        // Calculate the end-of-time-step mass flux
        rhoPhi = alphaPhi10*(rho1f - rho2f) + phi*rho2f;
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
