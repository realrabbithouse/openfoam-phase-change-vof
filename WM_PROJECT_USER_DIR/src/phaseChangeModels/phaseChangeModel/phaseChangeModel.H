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
    Foam::phaseChangeModel

Description
    吾魂兮无求乎永生 竭尽兮人事之所能

SourceFiles
    phaseChangeModel.C
    phaseChangeModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeModel_H
#define phaseChangeModel_H

#include "saturationTemperatureModel.H"
#include "fvMatrices.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class phaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeModel
:
    public dictionary
{

protected:

    // Protected Data

        //- Reference to twoPhaseMixtureThermo
        const twoPhaseMixtureThermo& mixture_;

        //- Reference to fvMesh
        const fvMesh& mesh_;

        //- Mass condensation and vaporisation rates of phase change [kg/(s.m^3)]
        //  "> 0" for condensation, "< 0" for evaporation
        volScalarField mDot_;

        //- Construct a saturationTemperatureModel by RTS
        autoPtr<saturationTemperatureModel> saturationTemperaturePtr_;

        // Switch isExplicit_;

        // Switch includeVolChange_;


public:

    //- Runtime type information
    TypeName("phaseChangeModel");

    // Dictionary name for constructing phaseChangeModel
    static const word phaseChangeDictName;

    // Constructors

        //- Construct from phase change model name, twoPhaseMixtureThermo and dictioanry
        phaseChangeModel
        (
            const twoPhaseMixtureThermo& mixture
        );

        //- Disallow default empty constructor
        phaseChangeModel() = delete;

        //- Disallow default bitwise copy construction
        phaseChangeModel(const phaseChangeModel&) = delete;


    // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseChangeModel,
            dictionary,
            (
                const twoPhaseMixtureThermo& mixture
            ),
            (mixture)
        );


    // Selectors
        static autoPtr<phaseChangeModel> New
        (
            const twoPhaseMixtureThermo& mixture
        );


    //- Destructor
    virtual ~phaseChangeModel();


    // Member Functions

        //- Return a const reference to the saturationTemperatureModel
        const saturationTemperatureModel& saturationTemperature() const;

        //- Return a const reference to the latentHeatModel
        const latentHeatModel& latentHeat() const;

        // -Return saturation temperature
        const volScalarField& TSat() const
        {
            return saturationTemperaturePtr_->TSat();
        }

        //- Return latent heat of vaporization
        const volScalarField& Lv() const
        {
            return saturationTemperaturePtr_->latentHeat().Lv();
        }

        virtual bool isExplicit() const
        {
            return false;
        }

        //- Update interface area density
        virtual void updateInterface(const volVectorField& normal);


        //- 隐式地实现相变源项

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate [kg/(s.m^3)]
        //  condensation >= 0, evaporation <= 0
        //  mDot = mDotAlphal()[0]*(1 - alphal) + mDotAlphal()[1]*alphal
        virtual Pair<tmp<volScalarField>> mDotAlphal() const = 0;

        //- Return the mass condensation and vaporisation rates as coefficients [kg/(s.m^3)]
        //  condensation >= 0, evaporation <= 0
        //  mDot = mDot()[0] + mDot()[1]
        virtual Pair<tmp<volScalarField>> mDot() const = 0;

        //- Return the volumetric condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate [1/s]
        virtual Pair<tmp<volScalarField>> vDotAlphal() const;

        //- Return the volumetric condensation and vaporisation rates as
        //  coefficients [1/s]
        virtual Pair<tmp<volScalarField>> vDot() const;

        //- Return the explicit phase change heat [W/m^3] -> [1 -1 -3 0 0]
        //  qDot = qDot()[0] + qDot()[1]
        virtual Pair<tmp<volScalarField>> qDot() const;

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (Tsat - T) for the condensation rate
        //  and a coefficient to multiply  (T - Tsat) for the vaporisation rate
        //  mDot = mDotDeltaT()[0]*(TSat - T) + mDotDeltaT()[1]*(TSat - T)
        //  mDotDeltaT()[0]*(TSat - T) >= 0; mDotDeltaT()[1]*(TSat - T) <= 0
        virtual Pair<tmp<volScalarField>> mDotDeltaT() const = 0;

        //- Source for T equation, implicit energy source
        virtual tmp<fvScalarMatrix> TSource() const = 0;


        //- 显式地实现相变源项

        //- Mass phase change source [kg/(s.m^3)]
        const volScalarField& Sm() const
        {
            return mDot_;
        }

        //- Phase change source term in phase fraction transport equation explicitly
        virtual tmp<volScalarField> Salphal() const;

        //- 更好的方式来处理相体积分数传输方程的源项, 使得矩阵求解时更加的robust
        virtual Pair<tmp<volScalarField>> alphalSuSp() const;

        //- Phase change source term in pressure Poisson's equation 
        virtual tmp<volScalarField> Sv() const;

        //- Phase change source term in energy equation
        virtual tmp<volScalarField> Sq() const;


        //- Correct the phase change properties
        virtual void correct() = 0;

        //- Read the phaseChnageModel dictionary and update
        virtual bool read(const dictionary& thermoDict);

        //- Is a interface area required phase change method?
        virtual bool isCalcInterface() const
        {
            return false;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseChangeModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "phaseChangeModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
