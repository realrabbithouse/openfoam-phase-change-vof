{
    // Temporarily making U relative to mesh motion
    if (mesh.moving())
    {
        U -= fvc::reconstruct(mesh.phi());
    }

    volScalarField::Internal divU
    (
        mesh.moving()
      ? fvc::div(phi + mesh.phi())
      : fvc::div(phi)
    );

    zeroField Su;
    zeroField Sp;

    if (phaseChange->isExplicit())
    {
        // Evaluate implicitly, in another way
        Pair<tmp<volScalarField>> alphalSuSp = 
            phaseChange->alphalSuSp();
        const volScalarField& alphalSu = alphalSuSp[0]();
        const volScalarField& alphalSp = alphalSuSp[1]();

        // Updating alpha1
        advector.advect
        (
            alphalSp(),
            (alphalSu() + divU*min(alpha1(), scalar(1)))()
        );

        // // Evaluate explicitly, for real
        // tmp<volScalarField> Salphal = phaseChange->Salphal();
        // const volScalarField& alphalSu = Salphal();

        // // Updating alpha1
        // advector.advect
        // (
        //     Sp,
        //     (Su + alphalSu() + divU*min(alpha1(), scalar(1)))()
        // );
    }
    else
    {
        Pair<tmp<volScalarField>> vDotAlphal = 
            phaseChange->vDotAlphal();
        const volScalarField& vDotcAlphal = vDotAlphal[0](); // Su
        const volScalarField& vDotvAlphal = vDotAlphal[1]();
        const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal); // Sp

        advector.advect
        (
            vDotvmcAlphal(),
            (divU*min(alpha1(), scalar(1)) + vDotcAlphal())()
        );
    }
    
    // Making U absolute again after advection step
    if (mesh.moving())
    {
        U += fvc::reconstruct(mesh.phi());
    }

    rhoPhi = advector.getRhoPhi(rho1, rho2);

    alphaPhi10 = advector.alphaPhi();
    /*
    In isoAdvection algorithm, alphaPhi_ is updated after isoAdvection::advect(Sp, Su)
    which satisfies alphaPhi_ = dVf_/mesh.time.deltaT()
    */

    alpha2 = 1.0 - alpha1;
}

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;
