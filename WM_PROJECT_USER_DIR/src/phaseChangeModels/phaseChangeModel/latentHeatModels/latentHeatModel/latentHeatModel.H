/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 AUTHOR,AFFILIATION
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
    Foam::latentHeatModel

Description

SourceFiles
    latentHeatModel.C
    latentHeatModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef latentHeatModel_H
#define latentHeatModel_H

#include "dictionary.H"
#include "twoPhaseMixtureThermo.H"
#include "dimensionedScalarFwd.H"
#include "volFieldsFwd.H"
#include "primitiveFieldsFwd.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class latentHeatModel Declaration
\*---------------------------------------------------------------------------*/

class latentHeatModel
:
    public dictionary
{

protected:

    // Protected Data

        //- Reference to twoPhaseMixtureThermo object
        const twoPhaseMixtureThermo& mixture_;

        //- Reference to the mesh objectRegistry
        const fvMesh& mesh_;

        //- Specific latent heat of vaporization, units [J/kg]
        volScalarField Lv_;

public:

    //- Runtime type information
    TypeName("latentHeatModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            latentHeatModel,
            dictionary,
            (
                const twoPhaseMixtureThermo& mixture
            ),
            (mixture)
        );


    // Static Data Members

        //- Dictionary name for accessing latent heat properties
        static const word latentHeatDictName;


    // Constructors

        //- Delete null constructor
        latentHeatModel() = delete;

        //- No copy construct
        latentHeatModel(const latentHeatModel&) = delete;

        //- No copy assignment
        latentHeatModel& operator=(const latentHeatModel&) = delete;

        //- Construct from components
        latentHeatModel(const twoPhaseMixtureThermo& mixture);


    // Selectors

        static autoPtr<latentHeatModel> New
        (
            const twoPhaseMixtureThermo& mixture
        );


    //- Destructor
    virtual ~latentHeatModel();


    // Member Functions

    volScalarField& Lv()
    {
        return Lv_;
    }

    const volScalarField& Lv() const
    {
        return Lv_;
    }

    virtual void correct() = 0;

    virtual bool read(const dictionary& thermoDict);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "latentHeatModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
