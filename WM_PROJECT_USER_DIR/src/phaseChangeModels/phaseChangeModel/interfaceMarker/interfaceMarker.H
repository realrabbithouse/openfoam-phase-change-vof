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

Class
    Foam::interfaceMarker

Description
    Helper class for marking interface cells for interface equilibrium phase
    change model.

SourceFiles
    interfaceMarker.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceMarker_H
#define interfaceMarker_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class interfaceMarker Declaration
\*---------------------------------------------------------------------------*/

class interfaceMarker
{
    // Private Data

        //- const reference to mesh
        const fvMesh& mesh_;

        //- const reference to volume fraction field
        const volScalarField& alpha1_;

        struct cellInfo
        {
            scalar alpha1Val; // volume fraction of a cell
            vector center; // cell center
        };

        struct faceInfo
        {
            label owner;
            label neighbour;
            vector normal;
            vector center;
        };

        List<cellInfo> cellsInfo_;
        List<faceInfo> facesInfo_;


public:

    // Public Member Data
        struct cellPairInfo
        {
            label faceIndex;
            label cell1, cell2; // cell index
            scalar val1, val2;
        };

    // Constructors

        //- Disallow default null constructor
        interfaceMarker() = delete;

        //- Construct from components
        interfaceMarker(const fvMesh& mesh, const volScalarField& alpha1);

        //- Disallow default bitwise copy construction
        interfaceMarker(const interfaceMarker&) = delete;


    //- Destructor
    ~interfaceMarker();


    // Member Functions
        //- reset volume fraction field
        void reset();

        //- One cell thick interface-cell layer
        void getInterfaceCells(labelList& interCells, const scalar interVal);

        //- two cells thick interface-cell layer
        void getDoubleLayerInterfaceCells(labelList& interCells, const scalar interVal);

        // Another double-layer method
        void getInterfaceCellPairs(List<cellPairInfo>& cellPairs, const scalar interVal);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "interfaceMarkerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
