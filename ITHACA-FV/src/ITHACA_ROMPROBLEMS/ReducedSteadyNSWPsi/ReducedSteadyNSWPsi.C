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

\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the reducedSteadyNS class

#include "ReducedSteadyNSWPsi.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedSteadyNSWPsi::reducedSteadyNSWPsi()
{
    para = ITHACAparameters::getInstance();
}

reducedSteadyNSWPsi::reducedSteadyNSWPsi(steadyNSWPsi& FOMproblem)
    :
    problem(&FOMproblem)
{
    //N_BC = problem->inletIndex.rows();
    Nphi_w = problem->AWPsi_matrix.rows();
    Nphi_psi_z = problem->BWPsi_matrix.cols();

    /* for (int k = 0; k < problem->liftfield.size(); k++)
    {
        Umodes.append((problem->liftfield[k]).clone());
    }

    for (int k = 0; k < problem->NUmodes; k++)
    {
        Umodes.append((problem->Umodes[k]).clone());
    }

    for (int k = 0; k < problem->NSUPmodes; k++)
    {
        Umodes.append((problem->supmodes[k]).clone());
    }
 */
    newton_object = newton_steadyNSWPsi(Nphi_w + Nphi_psi_z, Nphi_w + Nphi_psi_z, FOMproblem);
}

int newton_steadyNSWPsi::operator()(const Eigen::VectorXd& x,
                                Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_w);
    Eigen::VectorXd a_tmp(Nphi_w);
    Eigen::VectorXd b_tmp(Nphi_psi_z);
    a_tmp = x.head(Nphi_w);
    b_tmp = x.tail(Nphi_psi_z);
    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // MomW Term
    Eigen::VectorXd MA = problem->AWPsi_matrix * a_tmp * nu;
    // MomW Term
    Eigen::VectorXd MB = problem->BWPsi_matrix * b_tmp;
    // MassW Term
    Eigen::VectorXd MMW = problem->MW_matrix * a_dot;
    // MassPsi Term
    Eigen::VectorXd MMPsi = problem->MPsi_matrix * a_dot;
    // Pressure Term
    Eigen::VectorXd M3 = problem->P_matrix * a_tmp;
    // Penalty term
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_w, N_BC);

    // Term for penalty method
    /* if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = BC(l) * problem->bcVelVec[l] - problem->bcVelMat[l] *
                              a_tmp;
        }
    }*/

    for (int i = 0; i < Nphi_w; i++)
    {
        cc = a_tmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0,
                i) * a_tmp;
        fvec(i) = - MMW(i) + MA(i) - cc(0, 0);
    } 

    for (int j = 0; j < Nphi_psi_z; j++)
    {
        int k = j + Nphi_w;
        fvec(k) = -MMPsi(j) - MB(j);
    }

   /*  if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - BC(j);
        }
    }*/

    return 0; 
}


int newton_steadyNSWPsi::df(const Eigen::VectorXd& x,
                        Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_steadyNSWPsi> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}


// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //

/* void reducedSteadyNSWPsi::solveOnline_PPE(Eigen::MatrixXd vel_now)
{
    Info << "This function is still not implemented for the stationary case" <<
         endl;
    exit(0);
}

void reducedSteadyNSWPsi::solveOnline_sup(Eigen::MatrixXd vel)
{
    if (problem->bcMethod == "lift")
    {
        vel_now = setOnlineVelocity(vel);
    }
    else if (problem->bcMethod == "penalty")
    {
        vel_now = vel;
    }
    else
    {
        M_Assert(false,
                 "The BC method must be set to lift or penalty in ITHACAdict");
    }

    y.resize(Nphi_w + Nphi_psi_z, 1);
    y.setZero();

    // Change initial condition for the lifting function
    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
    }

    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    Eigen::HybridNonLinearSolver<newton_steadyNSWPsi> hnls(newton_object);
    newton_object.BC.resize(N_BC);
    newton_object.tauU = tauU;

    for (int j = 0; j < N_BC; j++)
    {
        newton_object.BC(j) = vel_now(j, 0);
    }

    newton_object.nu = nu;
    hnls.solve(y);
    Eigen::VectorXd res(y);
    newton_object.operator()(y, res);
    Info << "################## Online solve N° " << count_online_solve <<
         " ##################" << endl;

    if (Pstream::master())
    {
        std::cout << "Solving for the parameter: " << vel_now << std::endl;
    }

    if (res.norm() < 1e-5 && Pstream::master())
    {
        std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                  hnls.iter << " iterations " << def << std::endl << std::endl;
    }
    else if (Pstream::master())
    {
        std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                  hnls.iter << " iterations " << def << std::endl << std::endl;
    }

    count_online_solve += 1;
}
 */

// * * * * * * * * * * * * * * * Jacobian Evaluation  * * * * * * * * * * * * * //

/* void reducedSteadyNSWPsi::reconstruct_PPE(fileName folder, int printevery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextwrite = 0;

    for (int i = 0; i < online_solution.size(); i++)
    {
///////////modificare da qui
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", Umodes[0] * 0);

            for (int j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * online_solution[i](j + 1, 0);
            }

            ITHACAstream::exportSolution(U_rec, name(online_solution[i](0, 0)), folder);
            volScalarField P_rec("P_rec", problem->Pmodes[0] * 0);

            for (int j = 0; j < Nphi_p; j++)
            {
                P_rec += problem->Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }

            ITHACAstream::exportSolution(P_rec, name(online_solution[i](0, 0)), folder);
            nextwrite += printevery;
        }

        counter++;
    }
} */

void reducedSteadyNSWPsi::reconstruct(bool exportFields, fileName folder,
                                  int printevery)
{
    if (exportFields)
    {
        mkDir(folder);
        ITHACAutilities::createSymLink(folder);
    }

    int counter = 0;
    int nextwrite = 0;
    List <Eigen::MatrixXd> CoeffW;
    List <Eigen::MatrixXd> CoeffPsi;
    CoeffW.resize(0);
    CoeffPsi.resize(0);

    for (int i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            Eigen::MatrixXd currentWCoeff;
            Eigen::MatrixXd currentPsiCoeff;
            currentWCoeff = online_solution[i].block(1, 0, Nphi_w, 1);
            currentPsiCoeff = online_solution[i].bottomRows(Nphi_psi_z);
            CoeffW.append(currentWCoeff);
            CoeffPsi.append(currentPsiCoeff);
            nextwrite += printevery;
        }

        counter++;
    }

    volScalarField wRec("wRec", Wmodes[0]);
    volScalarField psiRec("psiRec", problem->Psi_zmodes[0]);
    wRecFields = problem->Wmodes.reconstruct(wRec, CoeffW, "wRec");
    psi_zRecFields = problem->Psi_zmodes.reconstruct(psiRec, CoeffPsi, "psiRec");

    if (exportFields)
    {
        ITHACAstream::exportFields(wRecFields, folder,
                                   "wRec");
        ITHACAstream::exportFields(psi_zRecFields, folder,
                                   "psi_zRec");
    }
}

/* double reducedSteadyNS::inf_sup_constant()
{
    double a;
    Eigen::VectorXd sup(Nphi_u);
    Eigen::VectorXd inf(Nphi_p);

    for (int i = 0; i < Nphi_p; i++)
    {
        for (int j = 0; j < Nphi_u; j++)
        {
            sup(j) = fvc::domainIntegrate(fvc::div(Umodes[j]) * Pmodes[i]).value() /
                     ITHACAutilities::H1Seminorm(Umodes[j]) / ITHACAutilities::L2Norm(Pmodes[i]);
        }

        inf(i) = sup.maxCoeff();
    }

    a = inf.minCoeff();
    return a;
} */


/* void reducedSteadyNS::reconstructLiftAndDrag(steadyNS& problem,
        fileName folder)
{
    mkDir(folder);
    system("ln -s ../../constant " + folder + "/constant");
    system("ln -s ../../0 " + folder + "/0");
    system("ln -s ../../system " + folder + "/system");
    int NUmodes = problem.NUmodes;
    int NSUPmodes = problem.NSUPmodes;
    int NPmodes = problem.NPmodes;
    int liftfieldSize = problem.liftfield.size();
    int totalSize = NUmodes + NSUPmodes + liftfieldSize;
    Eigen::VectorXd cl(online_solution.size());
    Eigen::VectorXd cd(online_solution.size());
    //Read FORCESdict
    IOdictionary FORCESdict
    (
        IOobject
        (
            "FORCESdict",
            "./system",
            Umodes[0].mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    fTau.setZero(online_solution.size(), 3);
    fN.setZero(online_solution.size(), 3);

    for (int i = 0; i < online_solution.size(); i++)
    {
        for (int j = 0; j < totalSize; j++)
        {
            fTau.row(i) += problem.tauMatrix.row(j) * online_solution[i](j + 1, 0);
        }

        for (int j = 0; j < NPmodes; j++)
        {
            fN.row(i) += problem.nMatrix.row(j) * online_solution[i](j + Nphi_u + 1, 0);
        }
    }

    // Export the matrices
    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(fTau, "fTau", "python", folder);
        ITHACAstream::exportMatrix(fN, "fN", "python", folder);
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(fTau, "fTau", "matlab", folder);
        ITHACAstream::exportMatrix(fN, "fN", "matlab", folder);
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(fTau, "fTau", "eigen", folder);
        ITHACAstream::exportMatrix(fN, "fN", "eigen", folder);
    }
}
 */
/* Eigen::MatrixXd reducedSteadyNS::setOnlineVelocity(Eigen::MatrixXd vel)
{
    assert(problem->inletIndex.rows() == vel.rows()
           && "Imposed boundary conditions dimensions do not match given values matrix dimensions");
    Eigen::MatrixXd vel_scal;
    vel_scal.resize(vel.rows(), vel.cols());

    for (int k = 0; k < problem->inletIndex.rows(); k++)
    {
        int p = problem->inletIndex(k, 0);
        int l = problem->inletIndex(k, 1);
        scalar area = gSum(problem->liftfield[0].mesh().magSf().boundaryField()[p]);
        scalar u_lf = gSum(problem->liftfield[k].mesh().magSf().boundaryField()[p] *
                           problem->liftfield[k].boundaryField()[p]).component(l) / area;
        vel_scal(k, 0) = vel(k, 0) / u_lf;
    }

    return vel_scal;
}
 */