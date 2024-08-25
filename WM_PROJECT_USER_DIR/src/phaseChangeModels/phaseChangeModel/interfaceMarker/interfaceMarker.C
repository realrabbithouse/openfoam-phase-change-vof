/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "interfaceMarker.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceMarker::interfaceMarker
(
    const fvMesh& mesh, const volScalarField& alpha1
)
:
    mesh_(mesh),
    alpha1_(alpha1),
    cellsInfo_(mesh.nCells()),
    facesInfo_(mesh.owner().size())
{
    forAll(mesh.cells(), cellI)
    {
        cellsInfo_[cellI].center = mesh.C()[cellI];
    }

    forAll(mesh.owner(), faceI)
    {
        facesInfo_[faceI].owner = mesh.owner()[faceI];
        facesInfo_[faceI].neighbour = mesh.neighbour()[faceI];
        facesInfo_[faceI].normal = mesh.Sf()[faceI]/mesh.magSf()[faceI];
        facesInfo_[faceI].center = mesh.Cf()[faceI];
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceMarker::~interfaceMarker()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interfaceMarker::reset()
{
    if (mesh_.dynamic())
    {
        cellsInfo_.clear();
        facesInfo_.clear();

        cellsInfo_.setSize(mesh_.nCells());
        facesInfo_.setSize(mesh_.owner().size());

        forAll(mesh_.cells(), cellI)
        {
            cellsInfo_[cellI].center = mesh_.C()[cellI];
            cellsInfo_[cellI].alpha1Val = alpha1_[cellI];
        }

        forAll(mesh_.owner(), faceI)
        {
            facesInfo_[faceI].owner = mesh_.owner()[faceI];
            facesInfo_[faceI].neighbour = mesh_.neighbour()[faceI];
            facesInfo_[faceI].normal = mesh_.Sf()[faceI]/mesh_.magSf()[faceI];
            facesInfo_[faceI].center = mesh_.Cf()[faceI];
        }
    }
    else
    {
        forAll(mesh_.cells(), cellI)
        {
            cellsInfo_[cellI].alpha1Val = alpha1_[cellI];
        }
    }
}


void Foam::interfaceMarker::getInterfaceCells(labelList& interCells, const scalar interVal)
{
    label nFaces = facesInfo_.size();

    for (label faceI = 0; faceI < nFaces; faceI++)
    {
        // Get two connecting cells' cellInfo
        cellInfo& Own = cellsInfo_[facesInfo_[faceI].owner];
        cellInfo& Nei = cellsInfo_[facesInfo_[faceI].neighbour];

        // Get values and difference
        scalar OwnVal = Own.alpha1Val;
        scalar NeiVal = Nei.alpha1Val;
        scalar dVal = (NeiVal - OwnVal) + SMALL;

        // 非相界面
        if (dVal == 0.0)
        {
            continue;
        }

        // interVal在OwnVal和NeiVal的相对位置
        scalar relativeLoc = (interVal - OwnVal)/dVal;

        // interVal不在OwnVal和NeiVal之间
        if ( (relativeLoc < 0) || (relativeLoc > 1) )
        {
            continue;
        }

        // Check the interface on which side
        // relativeLoc -> 0时靠近Owner, relativeLoc -> 1时靠近neighbour
        vector interLoc = Own.center + (Nei.center - Own.center)*relativeLoc;

        // 靠近Owner时, side为正, 靠近Neighbour时, side为负
        scalar whichSide = 
            (facesInfo_[faceI].center - interLoc) & facesInfo_[faceI].normal;

        // Mark if current cell contains interface
        // 如果side >=0 把owner作为interface cell, 如果side <= 0 把neighbour作为interface cell
        if ( whichSide >= 0 )
        {
            interCells.append(facesInfo_[faceI].owner);
        }

        if ( whichSide <= 0 )
        {
            interCells.append(facesInfo_[faceI].neighbour);
        }
    }
}


void Foam::interfaceMarker::getDoubleLayerInterfaceCells
(
    labelList& interCells,
    const scalar interVal
)
{
    label nFaces = facesInfo_.size();

    for (label faceI = 0; faceI < nFaces; faceI++)
    {
        // Get two connecting cells' cellInfo
        cellInfo& Own = cellsInfo_[facesInfo_[faceI].owner];
        cellInfo& Nei = cellsInfo_[facesInfo_[faceI].neighbour];

        // Get values and difference
        scalar OwnVal = Own.alpha1Val;
        scalar NeiVal = Nei.alpha1Val;
        scalar dVal = (NeiVal - OwnVal) + SMALL;

        // 非相界面
        if (dVal == 0.0)
        {
            continue;
        }

        // interVal在OwnVal和NeiVal的相对位置
        scalar relativeLoc = (interVal - OwnVal)/dVal;

        // interVal不在OwnVal和NeiVal之间
        if ( (relativeLoc < 0) || (relativeLoc > 1) )
        {
            continue;
        }

        interCells.append(facesInfo_[faceI].owner);
        interCells.append(facesInfo_[faceI].neighbour);
    }
}


void Foam::interfaceMarker::getInterfaceCellPairs
(
    List<cellPairInfo>& cellPairs,
    const scalar interVal
)
{
    label nFaces = facesInfo_.size();

    for (label faceI = 0; faceI < nFaces; faceI++)
    {
        // Get two connecting cells' cellInfo
        cellInfo& Own = cellsInfo_[facesInfo_[faceI].owner];
        cellInfo& Nei = cellsInfo_[facesInfo_[faceI].neighbour];

        // Get values and difference
        scalar OwnVal = Own.alpha1Val;
        scalar NeiVal = Nei.alpha1Val;
        scalar dVal = (NeiVal - OwnVal) + SMALL;

        // 非相界面
        if (dVal == 0.0)
        {
            continue;
        }

        // interVal在OwnVal和NeiVal的相对位置
        scalar relativeLoc = (interVal - OwnVal)/dVal;

        // interVal不在OwnVal和NeiVal之间
        if ((relativeLoc < 0) || (relativeLoc > 1))
        {
            continue;
        }

        cellPairInfo interCellPair;
        interCellPair.faceIndex = faceI;

        if (Own.alpha1Val >= Nei.alpha1Val)
        {
            // 注意这里没有按照owner/neighbour准则标记
            interCellPair.cell1 = facesInfo_[faceI].owner;
            interCellPair.cell2 = facesInfo_[faceI].neighbour;
            interCellPair.val1 = Own.alpha1Val;
            interCellPair.val2 = Nei.alpha1Val;
        }
        else
        {
            interCellPair.cell1 = facesInfo_[faceI].neighbour;
            interCellPair.cell2 = facesInfo_[faceI].owner;
            interCellPair.val1 = Nei.alpha1Val;
            interCellPair.val2 = Own.alpha1Val;
        }

        cellPairs.append(interCellPair);
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
