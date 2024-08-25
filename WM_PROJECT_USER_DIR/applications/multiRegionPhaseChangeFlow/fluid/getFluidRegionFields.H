    fvMesh& mesh = fluidRegions[i];

    volScalarField& p_rgh = p_rghFluid[i];
    volVectorField& U = UFluid[i];
    surfaceVectorField& Uf = UfFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    twoPhaseMixtureThermo& mixture = mixtureFluid[i];
    phaseChangeModel& phaseChange = phaseChangeFluid[i];

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    const volScalarField& rho1 = mixture.thermo1().rho();
    const volScalarField& rho2 = mixture.thermo2().rho();

    volScalarField& rho = rhoFluid[i];
    dimensionedScalar& pMin = pMinFluid[i];

    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    surfaceScalarField& rhoPhi = rhoPhiFluid[i];
    volScalarField& dgdt = dgdtFluid[i];
    volScalarField& rhoCv = rhoCvFluid[i];

    surfaceScalarField& alphaPhi10 = alphaPhi10Fluid[i];

    isoAdvection& advector = advectorFluid[i];
    compressibleInterPhaseTransportModel& turbulence = turbulenceFluid[i];
    volScalarField& K = KFluid[i];
    IOMRFZoneList& MRF = MRFfluid[i];
    fv::options& fvOptions = fluidFvOptions[i];

    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();
    const volScalarField& psi1 = mixture.thermo1().psi();
    const volScalarField& psi2 = mixture.thermo2().psi();

    bool frozenFlow = frozenFlowFluid[i];

    