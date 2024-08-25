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
    cellDims

Description
    Find the characteristic length (minimum/maximum edge length) for mesh cells.

\*---------------------------------------------------------------------------*/

volScalarField::Internal minCellDims // 返回每个cell的最短边
(
    IOobject
    (
        "cellMinDims",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar(dimLength, Zero)
);

volScalarField::Internal maxCellDims // 返回每个cell的最长边
(
    IOobject
    (
        "MaxCellDims",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar(dimLength, Zero)
);

const cellList& cells = mesh.cells();

forAll(cells, cellI)
{
    scalar minEdge = great;
    scalar maxEdge = small;

    const labelList& edges = mesh.cellEdges()[cellI]; // edge list of a cell

    forAll(edges, edgeI)
    {
        scalar edgeLength = mesh.edges()[ edges[edgeI] ].mag( mesh.points() );

        if (edgeLength < minEdge)
        {  
            minEdge = edgeLength;
        }
        if (edgeLength > maxEdge)
        {
            maxEdge = edgeLength;
        }
    }

    minCellDims[cellI] = minEdge;
    maxCellDims[cellI] = maxEdge;
}

Info<< "Minimum cell dimension = " <<  gMin(minCellDims) << endl;
