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
    psizeta
Description
    Reduction class for a NON-Stationary NS problem
SourceFiles
    psizeta.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the psizeta class.
/// \dir
/// Directory containing the header and source files for the psizeta class.

#ifndef psizeta_H
#define psizeta_H
#include "fvCFD.H"
//#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "steadyNS.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyNS Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady NS problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady NS problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the steadyNS class */
class psizeta: public steadyNS
{
    public:
        // Constructors
        /// Construct Null
        psizeta();

        /// Construct with argc and argv
        psizeta(int argc, char* argv[]);

        // Specific variable for the unstationary case
        /// Start Time (initial time to start storing the snapshots)
        scalar startTime;

        /// Final time (final time of the simulation and consequently of the acquisition of the snapshots)
        scalar finalTime;

        /// Time step of the simulation
        scalar timeStep;

        /// Time step of the writing procedure
        scalar writeEvery = timeStep;

        /// Auxiliary variable to store the next writing instant
        scalar nextWrite;

        // Dummy variables to transform pimplefoam into a class
        /// pimpleControl
        autoPtr<pimpleControl> _pimple;

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


       //dimensionedScalar nu;
        autoPtr<dimensionedScalar> _nu;

        /// List of pointers used to form the psi snapshots matrix
        PtrList<volScalarField> psifield;

        /// List of pointers used to form the zeta snapshots matrix
        PtrList<volScalarField> zetafield;

       /// List of pointers used to form the psivec snapshots matrix
        PtrList<volVectorField> psivecfield;

        /// List of pointers used to form the psi modes
        volScalarModes psimodes;

        /// List of pointers used to form the forcing modes
        PtrList<volScalarField> forcmodes;

        /// List of pointers used to form the psi vec modes
        volVectorModes psivecmodes;

        /// List of pointers used to form the zeta modes
        volScalarModes zetamodes;

        /// Number of psi modes used for the projection
        label Npsimodes;

        /// Number of psi modes used for the projection
        label Npsivecmodes;

        /// Number of zeta modes used for the projection
        label Nzetamodes;

        autoPtr<volScalarField> _zeta;

        autoPtr<volScalarField> _zeta0;

        autoPtr<volScalarField> _psi;

        autoPtr<volScalarField> _psi0;

        autoPtr<volVectorField> _psi_vec;

        autoPtr<volVectorField> _psi_vec0;


        /// Divergence of momentum zeta
        Eigen::Tensor<double, 3 > zetaTensor;

        /// Diffusion term for zeta
        Eigen::MatrixXd BZeta_matrix;

        /// Diffusion term for psi
        Eigen::MatrixXd BPsi_matrix;

        /// Mass Matrix for zeta
        Eigen::MatrixXd Zeta_matrix;

        /// Hybrid Mass Matrix for zeta
        Eigen::MatrixXd Mixed_matrix;

        /// forcing Matrix 
        Eigen::MatrixXd forcing_Matrix;



        /// Diffusion term for zeta
        Eigen::MatrixXd diffusive_term_zeta(label Nzetamodes);

        /// Diffusion term for psi
        Eigen::MatrixXd diffusive_term_psi(label Npsimodes);

        /// Mass term for zeta
        Eigen::MatrixXd mass_term_zeta(label Nzetamodes);

        /// Hybrid Mass 
        Eigen::MatrixXd hybrid_term(label Nzetamodes, label Npsimodes);

        /// convective term
        Eigen::Tensor<double, 3 > divMomentum_zeta(label Nzetamodes, label Npsivecmodes);

        /// forcing term
        Eigen::MatrixXd forcingMatrix(label Nzetamodes);


        void restart();

        // Functions

        //--------------------------------------------------------------------------
        /// @brief      Perform a truthsolve
        /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
        /// to construct mu_interp matrix which is written out in a specified folder, also for par
        /// file in the Parameters folder.
        ///
        void truthSolve(List<scalar> mu_now,
                        fileName folder = "./ITHACAoutput/Offline/");


       void project(fileName folder, label Npsimodes, label Nzetamodes, label Npsivecmodes);

        //--------------------------------------------------------------------------
        /// Function to check if the solution must be exported.
        ///
        /// @param      timeObject  The time object of OpenFOAM.
        ///
        /// @return     1 if we must write 0 elsewhere.
        ///
        bool checkWrite(Time& timeObject);
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif










