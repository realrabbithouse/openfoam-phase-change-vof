{
    const surfaceScalarField rhoCvPhi
    (
        "rhoCvPhi",
        fvc::interpolate(rhoCv/rho)*rhoPhi
    );

    if (phaseChange.isExplicit())
    {
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCv, T) + fvm::div(rhoCvPhi, T)
          + fvc::ddt(rho, K) + fvc::div(rhoPhi, K) - contErr*K
          - fvm::Sp(fvc::ddt(rhoCv) + fvc::div(rhoCvPhi), T)
          - fvm::laplacian(turbulence.kappaEff(), T)
          + divUp - contErr/rho*p
          - phaseChange.Sq()
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
          - fvm::laplacian(turbulence.kappaEff(), T)
          + divUp - contErr/rho*p
          + phaseChange.TSource()
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

    Info<< "Region: " << mesh.name() << " Min/max T:" << min(T).value()
        << ' ' << max(T).value() << endl;
}

