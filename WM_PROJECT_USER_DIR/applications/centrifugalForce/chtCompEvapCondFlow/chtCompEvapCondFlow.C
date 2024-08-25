/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    multiRegionPhaseChangeFlow

Group
    grpMultiPhaseSolvers

Description
    Transient solver for buoyant, turbulent fluid flow and solid heat
    conduction with conjugate heat transfer between solid and fluid regions.

    It handles secondary fluid or solid circuits which can be coupled
    thermally with the main fluid region. i.e radiators, etc.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "subCycle.H"
#include "compressibleInterPhaseTransportModel.H"
// #include "dynamicRefineFvMesh.H"
#include "isoAdvection.H"
#include "twoPhaseMixtureThermo.H"
#include "phaseChangeModel.H"
#include "alphaCourantNo.H"
#include "fluidCourantNo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H" // adjustTimeStep(false), maxCo(1), maxDeltaT(GREAT)
    #include "readSolidTimeControls.H" // maxDi(10)
    #include "fluidRegionCourantNo.H" // 计算fluid region的Courant No.
    #include "solidRegionDiffusionNo.H" // 计算solid region的Fourier No.
    #include "setInitialMultiRegionDeltaT.H"

    while (runTime.run())
    {
        #include "readTimeControls.H" // adjustTimeStep, maxCo, maxDeltaT
        #include "readSolidTimeControls.H" // maxDi(10)

        // 一共有三处fvSolution, 全局的fvSloution只需定义nOuterCorrectors, 默认为1
        #include "readPIMPLEControls.H"

        #include "fluidRegionCourantNo.H"
        #include "fluidRegionAlphaCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // 不确定有没有必要
        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                p_rghFluid[i].storePrevIter();
                rhoFluid[i].storePrevIter();
            }
        }

        // 更新相变模型
        forAll(fluidRegions, i)
        {
        	if (phaseChangeFluid[i].isCalcInterface())
        	{
        		advectorFluid[i].surf().reconstruct();
        		phaseChangeFluid[i].updateInterface
        		(
        			advectorFluid[i].surf().normal()
        		);
        		phaseChangeFluid[i].correct();
        	}
        	else
        	{
        		phaseChangeFluid[i].correct();
        	}
        }

        // --- PIMPLE loop
        for (int oCorr = 0; oCorr < nOuterCorr; ++oCorr)
        {
            const bool finalIter = (oCorr == nOuterCorr - 1);

            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "getFluidRegionFields.H"
                #include "readFluidMultiRegionPIMPLEControls.H"
                /*
                    const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

                    const int nCorr =
                        pimple.getOrDefault<int>("nCorrectors", 1);

                    const int nNonOrthCorr =
                        pimple.getOrDefault<int>("nNonOrthogonalCorrectors", 0);

                    const bool momentumPredictor =
                        pimple.getOrDefault("momentumPredictor", true);
                */
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "getSolidRegionFields.H"
                #include "readSolidMultiRegionPIMPLEControls.H"
                /*
                const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

                int nNonOrthCorr =
                    pimple.getOrDefault<int>("nNonOrthogonalCorrectors", 0);
                */
                #include "solveSolid.H"
            }

            // Additional loops for energy solution only
            // Only for the first PIMPLE Outer Corrector
            if (!oCorr && nOuterCorr > 1)
            {
                loopControl looping(runTime, pimple, "energyCoupling");

                while (looping.loop())
                {
                    Info<< nl << looping << nl;

                    forAll(fluidRegions, i)
                    {
                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                        #include "getFluidRegionFields.H"
                        #include "readFluidMultiRegionPIMPLEControls.H"

                        frozenFlow = true;
                        #include "solveFluid.H"
                    }

                    forAll(solidRegions, i)
                    {
                        Info<< "\nSolving for solid region "
                            << solidRegions[i].name() << endl;
                        #include "getSolidRegionFields.H"
                        #include "readSolidMultiRegionPIMPLEControls.H"
                        #include "solveSolid.H"
                    }
                }
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
