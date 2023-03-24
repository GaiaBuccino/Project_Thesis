/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------
License
    This file is part of ITHACA-FV
    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
Class
    unsteadyNS
Description
    Reduction class for a NON-Stationary NS problem
SourceFiles
    unsteadyNS.C
\*---------------------------------------------------------------------------*/


#ifndef unsteadyNSWPsi_H
#define unsteadyNSWPsi_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "steadyNSWPsi.H"
#include "UnsteadyProblem.H"
//#include "unsteadyNS.H"
#include <iostream>


class unsteadyNSWPsi: public steadyNSWPsi, public UnsteadyProblem
{
    public:
        // Constructors
        /// Construct Null
        unsteadyNSWPsi();

        /// Construct with argc and argv
        unsteadyNSWPsi(int argc, char* argv[]);

        // Dummy variables to transform pimplefoam into a class
        /// pimpleControl
        autoPtr<pimpleControl> _pimple;
        /// Flux
        //autoPtr<surfaceScalarField> _phi;

        /// Turbulence model
        autoPtr<incompressible::turbulenceModel> turbulence;

        /// adjustTimeStep
        bool adjustTimeStep;

        /// maxCourant
        scalar maxCo;

        /// maxDeltaT
        scalar maxDeltaT;

        // Counter used for the time-dependent boundary conditions
        label counter2 = 1;

        /// Time-dependent Boundary Method
        word timedepbcMethod;

        // Matrix containing the time-dependent boundary conditions for each time step
        Eigen::MatrixXd timeBCoff;

        // Time derivative numerical scheme order
        word timeDerivativeSchemeOrder;

        void truthSolve(List<scalar> mu_now, fileName folder = "./ITHACAoutput/Offline/");

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif










