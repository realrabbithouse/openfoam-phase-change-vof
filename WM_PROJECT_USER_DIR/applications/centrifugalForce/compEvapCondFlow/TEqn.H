{
    const surfaceScalarField rhoCvPhi
    (
        "rhoCvPhi",
        fvc::interpolate(rhoCv/rho)*rhoPhi
    );

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

