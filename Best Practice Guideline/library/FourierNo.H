/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

Global
    FourierNo

Description
    Calculate Fourier number for conduction

\*---------------------------------------------------------------------------*/

scalar maxFo
(
    runTime.controlDict().lookupOrDefault<scalar>("maxFo", 100.0)
);

scalar FoNum = 0.0;
scalar meanFoNum = 0.0;

{
    scalarField alphaEff // 热扩散系数
    (
        turbulence.alphaEff()().primitiveField()
    );

    FoNum = gMax
    (
        alphaEff*runTime.deltaTValue()/pow(minCellDims.field(), 2)
    );

    meanFoNum = gAverage
    (
        alphaEff*runTime.deltaTValue()/pow(minCellDims.field(), 2)
    );
}

Info<< "Fourier number mean: " << meanFoNum
    << " max: " << FoNum << endl;

// ************************************************************************* //
