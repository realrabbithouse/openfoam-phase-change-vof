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
    Foam::saturationTemperatureModel

Description

SourceFiles
    saturationTemperatureModel.C
    saturationTemperatureModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef saturationTemperatureModel_H
#define saturationTemperatureModel_H

#include "latentHeatModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class saturationTemperatureModel Declaration
\*---------------------------------------------------------------------------*/

class saturationTemperatureModel
:
    public dictionary
{

protected:

    // Protected Data

        //- Reference to the twoPhaseMixtureThermo
        const twoPhaseMixtureThermo& mixture_;

        //- Reference to the mesh objectRegistry
        const fvMesh& mesh_;

        //- Reference to the pressure field p
        const volScalarField& p_;

        //- Pointer to the latentHeatModel type
        autoPtr<latentHeatModel> latentHeatPtr_;

        //- Saturation tempurature, units [T]
        volScalarField TSat_;


public:

    //- Runtime type information
    TypeName("saturationTemperatureModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            saturationTemperatureModel,
            dictionary,
            (
                const twoPhaseMixtureThermo& mixture
            ),
            (mixture)
        );

    // Static Data Members

        //- Dictionary for accessing saturation temperature model properties
        static const word satTempDictName;


    // Constructors

        //- Construct from components
        saturationTemperatureModel(const twoPhaseMixtureThermo& mixture);

        //- Delete null constructor
        saturationTemperatureModel() = delete;

        //- Disallow default bitwise copy construction
        saturationTemperatureModel(const saturationTemperatureModel&) = delete;

        //- Disallow copy assignment operator
        saturationTemperatureModel& operator=(const saturationTemperatureModel&) = delete;


    // Selectors

        //- Selector prototype
        static autoPtr<saturationTemperatureModel> New
        (
            const twoPhaseMixtureThermo& mixture
        );


    //- Destructor
    virtual ~saturationTemperatureModel();


    // Member Functions

        //- Const access to the latentHeatModel
        const latentHeatModel& latentHeat() const;

        //- Return saturation tempurature field
        volScalarField& TSat()
        {
            return TSat_;
        }

        const volScalarField& TSat() const
        {
            return TSat_;
        }

        //- Correct saturation temperature
        virtual void correct() = 0;

        //- Read saturation properties from dictionary
        virtual bool read(const dictionary& thermoDict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "saturationTemperatureModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
