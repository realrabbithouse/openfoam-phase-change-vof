{
    rhoCv == alpha1*rho1*mixture.thermo1().Cv()
    + alpha2*rho2*mixture.thermo2().Cv();

    // 在pEqn.H之前计算的话
    const surfaceScalarField rhoCvPhi
    (
        "rhoCvPhi",
        fvc::interpolate(rhoCv/rho)*rhoPhi
    );

    // // 在pEqn.H之后计算的话
    // const surfaceScalarField rhoCvPhi
    // (
    //     "rhoCvPhi",
    //     fvc::interpolate(rhoCv)*phi
    // );

    // // rhoPhi = (fvc::interpolate(alpha1)*fvc::interpolate(rho1 - rho2)
    // //           + fvc::interpolate(rho2))*phi;
    // rhoPhi = fvc::interpolate(rho)*phi;

    // volScalarField::Internal contErrAft
    // (
    //     (
    //         fvc::ddt(rho) + fvc::div(rhoPhi)
    //       - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
    //       - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
    //     )()
    // );

    // scalar sumLocalContErr = runTime.deltaTValue()*
    //         mag(contErrAft)().weightedAverage(mesh.V()).value();
    // scalar globalContErr = runTime.deltaTValue()*
    //         contErrAft.weightedAverage(mesh.V()).value();

    // Info<< "Min/max contErr: " << min(mag(contErrAft)).value() << ' '
    //         << max(mag(contErrAft)).value() << endl;
    // Info<< "sumLocalContErr = " << sumLocalContErr
    //     << " globalContErr = " << globalContErr << endl;

    if (phaseChange->isExplicit())
    {
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCv, T) + fvm::div(rhoCvPhi, T)
          + fvc::ddt(rho, K) + fvc::div(rhoPhi, K) - contErr*K
          - fvm::Sp(fvc::ddt(rhoCv) + fvc::div(rhoCvPhi), T)
          // - fvm::Sp(contErr*rhoCv/rho, T)
          - fvm::laplacian(turbulence.kappaEff(), T)
          + divUp - contErr/rho*p
          - phaseChange->Sq()
         ==
            fvOptions(rhoCv, T)
        );

        TEqn.relax();

        fvOptions.constrain(TEqn);

        TEqn.solve();

        fvOptions.correct(T);
    }
    else
    {
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCv, T) + fvm::div(rhoCvPhi, T)
          + fvc::ddt(rho, K) + fvc::div(rhoPhi, K) - contErr*K
          - fvm::Sp(fvc::ddt(rhoCv) + fvc::div(rhoCvPhi), T)
          // - fvm::Sp(contErr*rhoCv/rho, T)
          - fvm::laplacian(turbulence.kappaEff(), T)
          + divUp - contErr/rho*p
          + phaseChange->TSource()
         ==
            fvOptions(rhoCv, T)
        );

        TEqn.relax();

        fvOptions.constrain(TEqn);

        TEqn.solve();

        fvOptions.correct(T);
    }

    // 根据压力p和温度T计算内能he, 然后根据压力p和内能he更新每一相的T, psi, rho, mu和alpha
    mixture.correctThermo();

    // 更新mixture的psi, mu和alpha (thermal diffusivity), 以及更新interfaceProperties (PS: 并没有更新mixture的密度)
    // 实际上更新interfaceProperties什么都不会做, 因为表面张力仅与相体积分数有关
    // 只是sigma可能与温度有关系
    mixture.correct();

    Info<< "Min/max T: " << min(T).value() << ' '
        << max(T).value() << endl;

    // // 测试相变能量源项
    // Pair<tmp<volScalarField>> mDotDeltaT = phaseChange->mDotDeltaT();
    // const volScalarField& mDotcDeltaT = mDotDeltaT[0]();
    // const volScalarField& mDotvDeltaT = mDotDeltaT[1]();

    // const volScalarField& Lv = phaseChange->Lv();
    // const volScalarField& TSat = phaseChange->TSat();

    // qDot_ = Lv*(mDotcDeltaT + mDotvDeltaT)*(TSat - T.oldTime());
    // Sq_ = phaseChange->Sq();
}
