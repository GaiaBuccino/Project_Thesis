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
Description
    Example of an unsteady NS Reduction Problem
SourceFiles
    04psizeta.C
\*---------------------------------------------------------------------------*/

#include "psizeta.H"
#include "ITHACAPOD.H"
#include "ReducedPsiZeta.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>

/*#include <datatable.h>
#include <bspline.h>
#include <bsplinebuilder.h>
#include <spline.h>
#include <Eigen/Dense>
#include "rbfspline.h"
#include "linearsolvers.h"
#include <Eigen/Eigen>

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "steadyNS.H"
#include <iostream>
#include "scalarIOField.H"
#include <datatable.h>
#include <bspline.h>
#include <bsplinebuilder.h>
#include <rbfspline.h>
#include <spline.h>
*/

class tutorial28: public psizeta
{
    public:
        explicit tutorial28(int argc, char* argv[])
            :
            psizeta(argc, argv),
            U(_U()),
            psi(_psi()),
            zeta(_zeta())
        {}

        // Fields To Perform
        volScalarField& zeta;
        volScalarField& psi;
        volVectorField& U;
        //volScalarField& lambda;

        void offlineSolve()
        {
            Vector<double> inl(0, 0, 0);
            scalar pinl = 0.0;
            List<scalar> mu_now(1);


            if (offline)
            {
                ITHACAstream::read_fields(psifield, psi, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(zetafield, zeta, "./ITHACAoutput/Offline/");
                ITHACAstream::read_fields(Ufield, U, "./ITHACAoutput/Offline/");
                //mu_samples =
                  //  ITHACAstream::readMatrix("./ITHACAoutput/Offline/mu_samples_mat.txt");
                //Lambdafield = Lambdafield*1e3;
            }
            else
            {
                for (label i = 0; i < mu.cols(); i++)
                {
                    //inl[0] = mu(0, i);
                    mu_now[0] = mu(0, i);
                     //assignBC(U, 4, inl);
                     //assignBC(Uevolve, 4, inl);
                     /*for (label i = 0; i < U.boundaryField()[4].size(); i++)
                     {
                       U.boundaryFieldRef()[4][i] = inl;
                     }

                     for (label i = 0; i < Uevolve.boundaryField()[4].size(); i++)
                     {
                       Uevolve.boundaryFieldRef()[4][i] = inl;
                     }

                     assignIF(Uevolve, inl);
                     assignIF(U, inl);
                     assignIF(p, pinl);
                     assignIF(lambda, pinl);*/
                    // change_viscosity( mu(0, i));
                    //std::cout << "viscosity:" << mu(0, i) << endl;
                    //change_filtering_radius( mu(0, i));
                    //assignIF(U,inl);
                    /*fvMesh& mesh = _mesh();

                    for (label i = 0; i < zeta.internalField().size(); i++)
                    {      
                           scalar x = mesh.C()[i].x();
                           scalar y = mesh.C()[i].y();
                           zeta.ref()[i] = Foam::exp(-3.14*(sqr(x-2.355) + sqr(y-3.14))) + Foam::exp(-3.14*(sqr(x-3.925) + sqr(y-3.14)));//vector((6.0/sqr(0.41))*y*(0.41-y), 0, 0); 
                    }

                    for (label i = 0; i < psi.internalField().size(); i++)
                    {
                           //scalar x = mesh().C()[faceI].x();
                           //scalar y = mesh().C()[faceI].y();
                           psi.ref()[i] = 0.0; 
                    }
 */                   //assignIF(p,pinl);
                    //assignBC(U,4,inl);
                    truthSolve(mu_now);
                }
            }
        }
};

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // Construct the tutorial20 object
    tutorial28 example(argc, argv);
    // Read parameters from ITHACAdict file
    ITHACAparameters* para = ITHACAparameters::getInstance(example._mesh(),
                             example._runTime());
    int Nmodeszetaout = para->ITHACAdict->lookupOrDefault<int>("Nmodeszetaout", 30);
    int Nmodespsiout = para->ITHACAdict->lookupOrDefault<int>("Nmodespsiout", 30);
    int Nmodespsivecout = para->ITHACAdict->lookupOrDefault<int>("Nmodespsivecout", 30);
    //int NmodesLambdaout = para.ITHACAdict->lookupOrDefault<int>("NmodesLambdaout", 10);
   
    int Nmodeszetaproj = para->ITHACAdict->lookupOrDefault<int>("Nmodeszetaproj", 14);
    int Nmodespsiproj = para->ITHACAdict->lookupOrDefault<int>("Nmodespsiproj", 6);
    int Nmodespsivecproj = para->ITHACAdict->lookupOrDefault<int>("Nmodespsivecproj", 6);
    //int NmodesLambdaout = para.ITHACAdict->lookupOrDefault<int>("NmodesLambdaout", 10);

    /// Set the number of parameters
    example.Pnumber = 1;
    /// Set samples
    example.Tnumber = 1;
    /// Set the parameters infos
    example.setParameters();
    // Set the parameter ranges
    example.mu_range(0, 0) = 0.00125; //0,000666667; //0.00125;
    example.mu_range(0, 1) = 0.00125; //0,000666667; //0.00125; //0.00125;
    // Generate equispaced samples inside the parameter range
    example.genEquiPar();
    // Set the inlet boundaries where we have non homogeneous boundary conditions
    //example.inletIndex.resize(1, 2);
    //example.inletIndex(0, 0) = 4;
    //example.inletIndex(0, 1) = 0;
    //example.inletIndex(1, 0) = 2;
    //example.inletIndex(1, 1) = 0;
    // Time parameters
    example.startTime = 0;
    example.finalTime = 20;
    example.timeStep = 0.01; //4e-4; //4e-4; //4e-4; //0.0004;
    example.writeEvery = 0.08; //0.08; //0.08; //4e-3;
    // Perform The Offline Solve;
    example.offlineSolve();
    // Solve the supremizer problem
    // Search the lift function
    //example.liftSolve();

    //example.solvesupremizer("snapshots");
    //example.solvesupremizer_evolve("snapshots");
    //cout << "norma: " << ITHACAutilities::L2norm(example.liftfield[0]) << "\n";
    // Normalize the lifting function
    //ITHACAutilities::normalizeFields(example.liftfield);
    // Create homogeneous basis functions for velocity
    //example.computeLift(example.Umedfield, example.liftfield, example.Uomfield);
    //example.computeLift_evolve(example.Uevolvefield, example.liftfield, example.Uomevolvefield);
    //example.computeLift(exam####### Performing the POD using EigenDecomposition Uevolve_om #######
    //ple.Ufield, example.liftfield, example.Uomfield);
    //example.computeLift_evolve(example.Uevolvefield, example.liftfield, example.Uomevolvefield);
    //Normalize the lifting function
    //ITHACAutilities::normalizeFields(example.liftfield);
    // Perform a POD decomposition for velocity and pressure
    //    std::clock_t c_start = std::clock();

    ITHACAPOD::getModes(example.zetafield, example.zetamodes, example._zeta().name(), example.podex, 0, 0,Nmodeszetaout);
    ITHACAPOD::getModes(example.psifield, example.psimodes, example._psi().name(), example.podex, 0, 0,Nmodespsiout);
    ITHACAPOD::getModes(example.Ufield, example.psivecmodes, example._U().name(), example.podex, 0, 0,Nmodespsivecout);

    //std::clock_t c_end = std::clock();
    //scalar time_elapsed_ms =  std::milli>c_end-c_start).count(); //(c_end-c_start); /// CLOCKS_PER_SEC;
    //std::cout << "CPU time used: " <<  std::chrono::duration<double, std::milli>(c_end-c_start).count() << " ms\n";
    //

    //example.solvesupremizer("modes");
    //example.solvesupremizer_evolve("modes");

    //example.solvesupremizer("modes");
    //example.solvesupremizer_evolve("modes");
    
    
    /*ITHACAPOD::getModes(example.supevolvefield, example.supevolvemodes, example.podex,
                        example.supex, 1, NmodesSUPevolveout);*/

    //cout<<"prova"<<endl;    /*ITHACAPOD::getModes(example.supfield, example.supmodes, example.podex,
                       //example.supex, 1, NmodesSUPout);
   
    example.project("./Matrices", Nmodespsiproj, Nmodeszetaproj, Nmodespsivecproj);

    //cout<<"prova2"<<endl;
    reducedPsiZeta reduced(example);
    // Set values of the reduced stuff
    reduced.nu = 0.00125; //0.001;
    reduced.tstart = 0;
    reduced.finalTime = 20;
    reduced.dt = 0.01; //4e-4; //4e-4;//2.33657e-5;
    reduced.storeEvery = 0.08; //2e-3;
    reduced.exportEvery = 0.08; //2e-3;
    // Set the online velocity
    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = 0;
    std::clock_t c_start = std::clock();
    reduced.solveOnline(vel_now,0);
    std::clock_t c_end = std::clock();
    //scalar time_elapsed_ms =  std::milli>(c_end-c_start).count(); //(c_end-c_start); /// CLOCKS_PER_SEC;
    std::cout << "CPU time used: " <<  std::chrono::duration<double, std::milli>(c_end-c_start).count() << " ms\n";
    //
    // Reconstruct the solution and export it
    reduced.reconstruct("./ITHACAoutput/ReconstructionSUP/");
    exit(0);
}



/// \dir 04psizeta Folder of the turorial 4
/// \file
/// \brief Implementation of tutorial 4 for an unsteady Navier-Stokes problem

/// \example 04psizeta.C
/// \section intro_unsreadyNS Introduction to tutorial 4
/// In this tutorial we implement a parametrized unsteady Navier-Stokes 2D problem where the parameter is the kinematic viscosity.
/// The physical problem represents an incompressible flow passing around a very long cylinder. The simulation domain is rectangular
/// with spatial bounds of [-4, 30], and [-5, 5] in the X and Y directions, respectively. The cylinder has a radius of
/// 0.5 unit length and is located at the origin. The system has a prescribed uniform inlet velocity of 1 m/s which is constant through the whole simulation.
///
/// The following image illustrates the simulated system at time = 50 s and Re = 100.
/// \image html cylinder.png
///
/// \section code04 A detailed look into the code
///
/// In this section we explain the main steps necessary to construct the tutorial N°4
///
/// \subsection header ITHACA-FV header files
///
/// First of all let's have a look at the header files that need to be included and what they are responsible for.
///
/// The header files of ITHACA-FV necessary for this tutorial are: <psizeta.H> for the full order unsteady NS problem,
/// <ITHACAPOD.H> for the POD decomposition, <reducedpsizeta.H> for the construction of the reduced order problem,
/// and finally <ITHACAstream.H> for some ITHACA input-output operations.
///
/// \dontinclude 04psizeta.C
/// \skip psizeta
/// \until ITHACAstream
///4
/// \subsection classtutorial20 Definition of the tutorial20 class
///
/// We define the tutorial20 class as a child of the psizeta class.
/// The constructor is defined with members that are the fields need to be manipulated
/// during the resolution of the full order problem using pimpleFoam. Such fields are
/// also initialized with the same initial conditions in the solver.
/// \skipline tutorial20
/// \until {}
///
/// Inside the tutorial20 class we define the offlineSolve method according to the
/// specific parametrized problem that needs to be solved. If the offline solve has
/// been previously performed then the method just reads the existing snapshots from the Offline directory.
/// Otherwise it loops over all the parameters, changes the system viscosity with the iterable parameter
/// then performs the offline solve.
///
/// \skipline offlineSolve
/// \until }
/// \skipline else
/// \until }
/// \skipline }
/// \skipline }
///
/// We note that in the commented line we show that it is possible to parametrize the boundary conditions.
/// For further details we refer to the classes: reductionProblem, and psizeta.
///
/// \subsection main Definition of the main function
///
/// In this section we show the definition of the main function.
/// First we construct the object "example" of type tutorial20:
///4
/// \skipline example
///
/// Then we parse the ITHACAdict file to determine the number of modes
/// to be written out and also the ones to be used for projection of
/// the velocity, pressure, and the supremizer:
/// \skipline ITHACAparameters
/// \until NmodesSUPproj
///
/// we note that a default value can be assigned in case the parser did
/// not find the corresponding string in the ITHACAdict file.
///
/// Now we would like to perform 10 parametrized simulations where the kinematic viscosity
/// is the sole parameter to change, and it lies in the range of {0.1, 0.01} m^2/s equispaced.
/// Alternatively, we can also think of those simulations as that they are performed for fluid
/// flow that has Re changes from Re=10 to Re=100 with step size = 10. In fact, both definitions
/// are the same since the inlet velocity and the domain geometry are both kept fixed through all
/// simulations.
///
/// In our implementation, the parameter (viscosity) can be defined by specifying that
/// Nparameters=1, Nsamples=10, and the parameter ranges from 0.1 to 0.01 equispaced, i.e.
///
/// \skipline example.Pnumber
/// \until example.genEquiPar()
///
/// After that we set the inlet boundaries where we have the non homogeneous BC:
///
/// \skipline example.inlet
/// \until example.inletIndex(0, 1) = 0;
///
/// And we set the parameters for the time integration, so as to simulate 20 seconds for each
/// simulation, with a step size = 0.01 seconds, and the data are dumped every 1.0 seconds, i.e.
///
/// \skipline example.startTime
/// \until example.writeEvery
///
/// Now we are ready to perform the offline stage:
///
/// \skipline Solve()
///
/// and to solve the supremizer problem:
///
/// \skipline supremizer()
///
/// In order to search and compute the lifting function (which should be a step function of value
/// equals to the unitary inlet velocity), we perform the following:
///
/// \skipline liftSolve()
///
/// Then we create homogenuous basis functions for the velocity:
///
/// \skipline computeLift
///
/// After that, the modes for velocity, pressure and supremizers are obtained:
///
/// \skipline getModes
/// \until supfield
///
/// then the projection onto the POD modes is performed with:
///
/// \skipline projectSUP
///
/// Now that we obtained all the necessary information from the POD decomposition and the reduced matrices,
/// we are now ready to construct the dynamical system for the reduced order model (ROM). We proceed
/// by constructing the object "reduced" of type reducedpsizeta:
///
/// \skipline reducedpsizeta
///
/// And then we can use the new constructed ROM to perform the online procedure, from which we can simulate the
/// problem at new set of parameters. For instance, we solve the problem with a viscosity=0.055 for a 15
/// seconds of physical time:
///
/// \skipline reduced.nu
/// \until reduced.dt
///
/// and then the online solve is performed. In this tutorial, the value of the online velocity
/// is in fact a multiplication factor of the step lifting function for the unitary inlet velocity.
/// Therefore the online velocity sets the new BC at the inlet, hence we solve the ROM at new BC.
///
/// \skipline Eigen::
/// \until solveOnline_sup
///
/// Finally the ROM solution is reconstructed and exported:
///
/// \skipline reconstruct_sup
///
/// We note that all the previous evaluations of the pressure were based on the supremizers approach.
/// We can also use the Pressure Poisson Equation (PPE) instead of SUP so as to be implemented for the
/// projections, the online solve, and the fields reconstructions.
///
///
/// \section plaincode The plain program
/// Here there's the plain code
///
