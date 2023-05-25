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
/// Source file of the reducedPsiZeta class


#include "ReducedPsiZeta.H"
#include "psizeta.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor initialization
reducedPsiZeta::reducedPsiZeta()
{
}

reducedPsiZeta::reducedPsiZeta(psizeta& FOMproblem)
    :
    problem(&FOMproblem)
{
    //N_BC = problem->inletIndex.rows();
    Nphi_psi = problem->BPsi_matrix.rows();
    Nphi_zeta = problem->BZeta_matrix.rows();
    //Nphi_psivec = Nphi_psi;
    //interChoice = problem->interChoice;



    // Create locally the velocity modes


    for (label k = 0; k < problem->Npsimodes; k++)
    {
        psimodes.append(problem->psimodes[k].clone());
    }

    for (label k = 0; k < problem->Nzetamodes; k++)
    {
        zetamodes.append(problem->zetamodes[k].clone());
    }

    /*for (label k = 0; k < problem->Npsimodes; k++)
    {
        psivecmodes.append(problem->psivecmodes[k]);
    }*/

    /*for (label k = 0; k < 1; k++)
    {
        forcmodes.append(problem->forcmodes[k]);
    }*/


    newton_object = newton_psizeta(Nphi_psi + Nphi_zeta, Nphi_psi + Nphi_zeta,
                        FOMproblem); //TODO 

}

// * * * * * * * * * * * * * * * Operators PPE * * * * * * * * * * * * * * * //

// Operator to evaluate the residual for the Pressure Poisson Equation (PPE) approach
int newton_psizeta::operator()(const Eigen::VectorXd& x,
                                      Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd beta_psi(Nphi_psi);
    Eigen::VectorXd beta_zeta(Nphi_zeta);

    beta_zeta = x.head(Nphi_zeta);
    beta_psi = x.tail(Nphi_psi);
   
    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // EVOLVE: laplacian term
    Eigen::VectorXd Lap_Zeta = problem->BZeta_matrix *nu*beta_zeta;
    Eigen::VectorXd Lap_Psi = problem->BPsi_matrix *beta_psi;
    Eigen::VectorXd Mass_Zeta = problem->Zeta_matrix * beta_zeta;
    Eigen::VectorXd Mass_Zeta_old = problem->Zeta_matrix * y_old.head(Nphi_zeta);
    Eigen::VectorXd Mass_Zeta_old_old = problem->Zeta_matrix * yOldOld.head(Nphi_zeta);
    Eigen::VectorXd Mass_Hyb = problem->Mixed_matrix *y_old.head(Nphi_zeta); //beta_zeta;
    Eigen::VectorXd forcing = 0.1*problem->forcing_Matrix*Foam::exp(-nu*time); //beta_zeta;
    
    for (label i = 0; i < Nphi_zeta; i++)
    {         
         cc = y_old.tail(Nphi_psi).transpose()*Eigen::SliceFromTensor(problem->zetaTensor, 0, i)*beta_zeta;
         //fvec(i) = Mass_Zeta(i)*(1.5/dt)- Mass_Zeta_old(i)*(2.0/dt) + Mass_Zeta_old_old(i)*(0.5/dt)  - Lap_Zeta(i) + cc(0, 0);
         fvec(i) = Mass_Zeta(i)*(1.0/dt)- Mass_Zeta_old(i)*(1.0/dt) - Lap_Zeta(i) + cc(0, 0); //+ forcing(i);

    }


    for (label j = 0; j < Nphi_psi; j++)
    {
        label k = j + Nphi_zeta;
        fvec(k) = Lap_Psi(j) + Mass_Hyb(j);

    }
    

    /*for (label j = 0; j < N_BC; j++)
    {
        fvec(j) = x(j) - BC(j);
    }*/


    /*for (label j = Nphi_u_evolve + Nphi_p; j < Nphi_u_evolve + Nphi_p + N_BC; j++)
    {
        fvec(j) = x(j) - BC(j - Nphi_u_evolve - Nphi_p);
    }*/

    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newton_psizeta::df(const Eigen::VectorXd& x,
                              Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_psizeta> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}



// * * * * * * * * * * * * * * * Solve Functions PPE * * * * * * * * * * * * * //

void reducedPsiZeta::solveOnline(Eigen::MatrixXd vel,
                                        label startSnap)
{

    M_Assert(exportEvery >= dt,
             "The time step dt must be smaller than exportEvery.");
    M_Assert(storeEvery >= dt,
             "The time step dt must be smaller than storeEvery.");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt) == true,
             "The variable storeEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt) == true,
             "The variable exportEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / storeEvery) == true,
             "The variable exportEvery must be an integer multiple of the variable storeEvery.");
    int numberOfStores = round(storeEvery / dt);

    vel_now = vel;
    //vel_now = setOnlineVelocity(vel);
    // Create and resize the solution vector
    y.resize(Nphi_zeta + Nphi_psi, 1);
    //y.resize(Nphi_u + Nphi_u, 1);
    y.setZero();
     //startSnap = 2;

    y.head(Nphi_zeta) = ITHACAutilities::getCoeffs(problem->zetafield[startSnap],
                     zetamodes); //TODO
    y.tail(Nphi_psi) = ITHACAutilities::getCoeffs(problem->psifield[startSnap],
                     psimodes);
     //y.tail(Nphi_u) = y.head(Nphi_u); //ITHACAutilities::get_coeffs(problem->Ufield[startSnap], Umodes); //TODO
     /*y.tail(Nphi_lambda) = ITHACAutilities::get_coeffs(problem->Lambdafield[startSnap],
                     Lambdamodes);*/
    int nextStore = 0;
    int counter2 = 0;

    time = tstart;

    // Set some properties of the newton object
    newton_object.nu = nu;
    newton_object.y_old = y;
    newton_object.dt = dt;
    //newton_object.BC.resize(N_BC);
    newton_object.tauU = tauU;
    newton_object.yOldOld = newton_object.y_old;
    //newton_object_sup.aNut = anut0;

    // Set number of online solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores);
    online_solution.resize(onlineSize);
    //rbfCoeffMat.resize(nphi_a + 1, onlineSize + 3);
    // Set the initial time
    time = tstart;
    // Counting variable
    int counter = 0;
    // Create vector to store temporal solution and save initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_zeta + Nphi_psi  + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    online_solution[counter] = tmp_sol;
    counter ++;
    counter2++;
    nextStore += numberOfStores;
   
    // Create nonlinear solver object
    Eigen::HybridNonLinearSolver<newton_psizeta> hnls(newton_object);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    // Start the time loop
    while (time < finalTime)
    { 
        time = time + dt;
        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);


        newton_object.operator()(y, res);
        newton_object.yOldOld = newton_object.y_old;
        newton_object.y_old = y;

        std::cout << "################## Online solve N° " << counter <<
                  " ##################" << std::endl;
        Info << "Time = " << time << endl;

        if (res.norm() < 1e-6)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                      hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
                      hnls.iter << " iterations " << def << std::endl << std::endl;
        }
        //count_online_solve +=1; 
        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;

        if (counter == nextStore)
        {
            if (counter2 >= online_solution.size())
            {
                online_solution.append(tmp_sol);
            }
            else
            {
                online_solution[counter2] = tmp_sol;
            }

            //rbfCoeffMat(0, counter2) = time;
            //rbfCoeffMat.block(1, counter2, nphi_a, 1) = newton_object_PPE.aNut;
            nextStore += numberOfStores;
            counter2 ++;
        }

        counter ++;
    }

    // Export the solution
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "python",
                               "./ITHACAoutput/red_coeff");
    ITHACAstream::exportMatrix(online_solution, "red_coeff", "matlab",
                               "./ITHACAoutput/red_coeff");
    //count_online_solve +=1; 
}
void reducedPsiZeta::reconstruct(fileName folder, label startSnap)
{

    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextwrite = 0;
    int counter2 = 1;
    int exportEveryIndex = round(exportEvery / storeEvery);

    int ind = 1;

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            volScalarField zeta_rec("zeta", zetamodes[0] * 0);

            for (label j = 0; j < Nphi_zeta; j++)
            {
                zeta_rec += zetamodes[j] * online_solution[i](j + 1, 0);
                 // U_rec += Umodes[0];
            }

            ITHACAstream::exportSolution(zeta_rec,  name(counter2), folder);


           volScalarField psi_rec("psi", psimodes[0] * 0);

            for (label j = 0; j < Nphi_psi; j++)
            {
                
                label k = j + Nphi_zeta;
                psi_rec += psimodes[j] * online_solution[i](k + 1, 0);
            }

            ITHACAstream::exportSolution(psi_rec,  name(counter2), folder);

            nextwrite += exportEveryIndex;
            double timenow = online_solution[i](0, 0);
            std::ofstream of(folder + name(counter2) + "/" + name(timenow));
            counter2 ++;
            //ind++;
            ZETAREC.append(zeta_rec.clone());
            PSIREC.append(psi_rec.clone());


            //std::ofstream file;
            //file.open ("FOM_ROM_L2norm.txt", std::ofstream::out | std::ofstream::app);

            //scalar L2_U = 0.0;
            //scalar L2_U_FOM = 0.0;
            //scalar L2_Uevolve = 0.0;
            //scalar L2_Uevolve_FOM = 0.0;
            //scalar L2_p = 0.0;
            //scalar L2_p_FOM = 0.0;
            //scalar L2_lambda = 0.0;
            //scalar L2_lambda_FOM = 0.0;

            //scalar L2_K1 = 0.0;
            //scalar L2_K2 = 0.0;



           
            /*forAll(U_rec, cellI)
            {of new numerical methods for ocean flows. Thus, tounderstand the potential of using ROMs for the effi

                L2_U += sqr(mag(U_rec[cellI] - problem->Ufield[ind][cellI]))*problem->_mesh().V()[cellI];
                L2_U_FOM += sqr(mag(problem->Ufield[ind][cellI]))*problem->_mesh().V()[cellI];
                L2_Uevolve += sqr(mag(Uevolve_rec[cellI] - problem->Uevolvefield[ind][cellI]))*problem->_mesh().V()[cellI];
                L2_Uevolve_FOM += sqr(mag(problem->Uevolvefield[ind][cellI]))*problem->_mesh().V()[cellI];
                L2_p += sqr(mag(P_rec[cellI] - problem->Pfield[ind][cellI]))*problem->_mesh().V()[cellI];
                L2_p_FOM += sqr(mag(problem->Pfield[ind][cellI]))*problem->_mesh().V()[cellI];
                L2_lambda += sqr(mag(Lambda_rec[cellI] - problem->Lambdafield[ind][cellI]))*problem->_mesh().V()[cellI];
                L2_lambda_FOM += sqr(mag(problem->Lambdafield[ind][cellI]))*problem->_mesh().V()[cellI];
            }

                L2_U = L2_U/(L2_U_FOM + 1e-16);
                L2_Uevolve = L2_Uevolve/(L2_Uevolve_FOM + 1e-16);
                L2_p = L2_p/(L2_p_FOM + 1e-16);
                L2_lambda = L2_lambda/(L2_lambda_FOM + 1e-16);
           */



                /*L2_U = ITHACAutilities::L2norm(U_rec - problem->Ufield[ind])/(ITHACAutilities::L2norm(problem->Ufield[ind]) + 1e-16);
                L2_Uevolve = ITHACAutilities::L2norm(Uevolve_rec - problem->Uevolvefield[ind])/(ITHACAutilities::L2norm(problem->Uevolvefield[ind]) + 1e-16);
                L2_p = ITHACAutilities::L2norm(P_rec - problem->Pfield[ind])/(ITHACAutilities::L2norm(problem->Pfield[ind]) + 1e-16);*/
                //L2_lambda = ITHACAutilities::L2norm(Lambda_rec - problem->Lambdafield[ind])/(ITHACAutilities::L2norm(problem->Lambdafield[ind]) + 1e-16);

               // L2_K1 = (sqr(ITHACAutilities::L2norm(Uevolve_rec)) - sqr(ITHACAutilities::L2norm(problem->Uevolvefield[ind])))/(sqr(ITHACAutilities::L2norm(problem->Uevolvefield[ind])) + 1e-16);
                //L2_K2 = (sqr(ITHACAutilities::L2norm(U_rec)) - sqr(ITHACAutilities::L2norm(problem->Ufield[ind])))/(sqr(ITHACAutilities::L2norm(problem->Ufield[ind])) + 1e-16);

                /*L2_U = ITHACAutilities::L2norm(U_rec - problem->Umedfield[ind])/(ITHACAutilities::L2norm(problem->Umedfield[ind]) + 1e-16);
                L2_Uevolve = ITHACAutilities::L2norm(Uevolve_rec - problem->Uevolvefield[ind])/(ITHACAutilities::L2norm(problem->Uevolvefield[ind]) + 1e-16);
                L2_p = ITHACAutilities::L2norm(P_rec - problem->Pfield[ind])/(ITHACAutilities::L2norm(problem->Pfield[ind]) + 1e-16);
                L2_lambda = ITHACAutilities::L2norm(Lambda_rec - problem->Lambdafield[ind])/(ITHACAutilities::L2norm(problem->Lambdafield[ind]) + 1e-16);

                L2_K1 = (sqr(ITHACAutilities::L2norm(Uevolve_rec)) - sqr(ITHACAutilities::L2norm(problem->Uevolvefield[ind])))/(sqr(ITHACAutilities::L2norm(problem->Uevolvefield[ind])) + 1e-16);
                L2_K2 = (sqr(ITHACAutilities::L2norm(U_rec)) - sqr(ITHACAutilities::L2norm(problem->Umedfield[ind])))/(sqr(ITHACAutilities::L2norm(problem->Umedfield[ind])) + 1e-16);*/
       
           /*if (Pstream::master())
            {
                file << L2_U << "\t" << L2_Uevolve << "\t" << L2_p << std::endl;
            }*/
           
        }

        counter++;
        //ind++;
    }
}


Eigen::MatrixXd reducedPsiZeta::setOnlineVelocity(Eigen::MatrixXd vel)
{
    assert(problem->inletIndex.rows() == vel.rows()
           && "Imposed boundary conditions dimensions do not match given values matrix dimensions");
    Eigen::MatrixXd vel_scal;
    vel_scal.resize(vel.rows(), vel.cols());

    for (int k = 0; k < problem->inletIndex.rows(); k++)
    {
        label p = problem->inletIndex(k, 0);
        label l = problem->inletIndex(k, 1);
        scalar area = gSum(problem->liftfield[0].mesh().magSf().boundaryField()[p]);
        scalar u_lf = gSum(problem->liftfield[k].mesh().magSf().boundaryField()[p] *
                           problem->liftfield[k].boundaryField()[p]).component(l) / area;

        for (int i = 0; i < vel.cols(); i++)
        {
            vel_scal(k, i) = vel(k, i) / u_lf;
        }
    }

    return vel_scal;
}


//************************************************************************* //
