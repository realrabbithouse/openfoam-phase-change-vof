{
    // rhoCp = rho*cp;
    // kappaEff = thermo->kappa() + rho*cp*turbulence->nut()/Prt;
    // const surfaceScalarField rhoCpPhi
    // (
    //     "rhoCpPhi",
    //     rhoPhi*(Cp1 - Cp2) + phi*rho2*Cp2
    // );
    
    rhoCp == alpha1*rho1*mixture.thermo1().Cp()
        + alpha2*rho2*mixture.thermo2().Cp();

    const surfaceScalarField rhoCpPhi
    (
        "rhoCpPhi",
        fvc::interpolate(rhoCp/rho)*rhoPhi
    );

    if (phaseChange->isExplicit())
    {
        fvScalarMatrix TEqn
        (
             fvm::ddt(rhoCp, T)
           + fvm::div(rhoCpPhi, T)
           - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T) // 连续性误差
           - fvm::laplacian(turbulence.kappaEff(), T)
           // + phaseChange->TSource()
           - phaseChange->Sq()
           ==
             fvOptions(rhoCp, T)
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
             fvm::ddt(rhoCp, T)
           + fvm::div(rhoCpPhi, T)
           - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T) // 连续性误差
           - fvm::laplacian(turbulence.kappaEff(), T)
           + phaseChange->TSource()
           ==
             fvOptions(rhoCp, T)
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

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
