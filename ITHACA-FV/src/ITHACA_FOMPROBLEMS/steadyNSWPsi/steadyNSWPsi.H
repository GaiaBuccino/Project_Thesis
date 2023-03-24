#ifndef steadyNSWPsi_H
#define steadyNSWPsi_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "reductionProblem.H"
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#if OPENFOAM >= 1812
#include "ITHACAforces18.H"
#else
#include "ITHACAforces.H"
#endif
#include "volFields.H"
#include <iostream>
#include "IPstream.H"
#include "OPstream.H"
#include "Modes.H"
#include "steadyNS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class \ Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> steady NS problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a steady NS problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the reduction problem class */
class steadyNSWPsi: public steadyNS
{


    public:
        // Constructors
        /// Null constructor
        steadyNSWPsi();

        /// Construct with argc and argv
        steadyNSWPsi(int argc, char* argv[]);
        ~steadyNSWPsi() {};


        /// List of pointers used to form the vorticity snapshots matrix
        PtrList<volScalarField> Wfield;

        /// List of pointers used to form the Streamlines scalar snapshots matrix
        PtrList<volScalarField> Psi_zfield;

        /// List of pointers used to form the Streamlines vector snapshots matrix
        PtrList<volVectorField> Psifield;

        /// List of pointers used to form the temporal vector snapshots matrix
        PtrList<volVectorField> temp;

        /// List of pointers used to form the vorticity modes
        volScalarModes Wmodes;

        /// List of pointers used to form the streamlines scalar modes
        volScalarModes Psi_zmodes;

        /// List of pointers used to form the streamlines vector modes
        volVectorModes Psimodes;

        /* /// List of pointers used to the projection modes of W and Psi_z    
        volScalarModes W_proj_modes;
        volScalarModes Psi_z_proj_modes;
        volVectorModes Psi_proj_modes; */


        /// Number of vorticity modes used for the projection
        label NWmodes;

        /// Number of streamlines scalar modes used for the projection
        label NPsi_zmodes;

        /// Number of streamlines vector modes used for the projection
        label NPsimodes;

        /// Diffusion term
        Eigen::MatrixXd AWPsi_matrix;
        Eigen::MatrixXd BWPsi_matrix;

        /// Mass Matrix
        Eigen::MatrixXd MW_matrix;
        Eigen::MatrixXd MPsi_matrix;
        Eigen::Tensor<double, 3 > GWPsi_tensor;

        /// Vorticity field
        autoPtr<volScalarField> _W;

        autoPtr<volScalarField> _T;

        /// Streamlines scalar field
        autoPtr<volScalarField> _Psi_z;

        /// Streamlines vector field
        autoPtr<volVectorField> _Psi;

        /// temporal field
        autoPtr<volVectorField> _temp;

        /// Initial Vorticity field (for restart purposes)
        autoPtr<volScalarField> _W0;

        /// Initial Streamlines field (for restart purposes)
        autoPtr<volScalarField> _Psi_z0;

        /// Initial Streamlines field (for restart purposes)
        autoPtr<volVectorField> _Psi0;
        
        void truthSolve(List<scalar> mu_now);

        //Method overloading
        void projectSUP(fileName folder, label NWmodes, label NPsi_zmodes);
        Eigen::MatrixXd  massW_term(label NWmodes, label NPsi_zmodes);
        Eigen::MatrixXd  massPsi_z_term(label NWmodes, label NPsi_zmodes);
        Eigen::MatrixXd diffusiveW_term(label NWmodes, label NPsi_zmodes);
        Eigen::MatrixXd diffusivePsi_z_term(label NWmodes, label NPsi_zmodes);
        Eigen::Tensor<double, 3 > convective_term_tens(label NWmodes, label NPsi_zmodes);
        Eigen::Tensor<double, 3 > convective_term_tens(label NUmodes,label NPmodes,label NSUPmodes);
};
#endif


