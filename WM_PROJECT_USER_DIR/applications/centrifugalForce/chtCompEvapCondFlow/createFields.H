#include "createFluidFields.H"
#include "createSolidFields.H"

IOdictionary centrifugalForceDict
(
	IOobject
	(
		"centrifugalForceProperties",
		runTime.time().constant(),
		runTime.db(),
		IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
	)
);

dimensionedVector omega("omega", dimless/dimTime, centrifugalForceDict);
dimensionedVector origin("origin", dimLength, centrifugalForceDict);
