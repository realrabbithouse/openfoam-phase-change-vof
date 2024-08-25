// Update alpha1
#include "alphaSuSp.H"

// Pair<tmp<volScalarField>> vDotAlphal = 
// 	phaseChange.vDotAlphal();
// const volScalarField& vDotcAlphal = vDotAlphal[0](); // Su
// const volScalarField& vDotvAlphal = vDotAlphal[1]();
// const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal); // Sp

// advector.advect
// (
// 	(Sp + vDotvmcAlphal())(),
// 	(Su + divU*min(alpha1(), scalar(1)) + vDotcAlphal())()
// );

if (phaseChange.isExplicit())
{
	Pair<tmp<volScalarField>> alphalSuSp = 
	phaseChange.alphalSuSp();
	const volScalarField& alphalSu = alphalSuSp[0]();
	const volScalarField& alphalSp = alphalSuSp[1]();

    // Updating alpha1
	advector.advect
	(
		(Sp + alphalSp())(),
		(Su + alphalSu() + divU*min(alpha1(), scalar(1)))()
	);
}
else
{
	Pair<tmp<volScalarField>> vDotAlphal = 
		phaseChange.vDotAlphal();
	const volScalarField& vDotcAlphal = vDotAlphal[0](); // Su
	const volScalarField& vDotvAlphal = vDotAlphal[1]();
	const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal); // Sp

	advector.advect
	(
		(Sp + vDotvmcAlphal())(),
		(Su + divU*min(alpha1(), scalar(1)) + vDotcAlphal())()
	);
}

// Update rhoPhi
rhoPhi = advector.getRhoPhi(rho1, rho2);
alphaPhi10 = advector.alphaPhi();

alpha2 = 1.0 - alpha1;

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") - 1 = " << max(alpha1).value() - 1
    << endl;
