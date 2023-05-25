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

#include "psizeta.H"

/// \file
/// Source file of the psizeta class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
psizeta::psizeta() {}

// Construct from zero
psizeta::psizeta(int argc, char* argv[])
{
    _args = autoPtr<argList>
            (
                new argList(argc, argv)
            );

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
    _pimple = autoPtr<pimpleControl>
              (
                  new pimpleControl
                  (
                      mesh
                  )
              );
    ITHACAdict = new IOdictionary
    (
        IOobject
        (
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
#include "createFields.H"
#include "createFvOptions.H"
    para = ITHACAparameters::getInstance(mesh, runTime);
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "lift");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty",
             "The BC method must be set to lift or penalty in ITHACAdict");
    timedepbcMethod = ITHACAdict->lookupOrDefault<word>("timedepbcMethod", "no");
    M_Assert(timedepbcMethod == "yes" || timedepbcMethod == "no",
             "The BC method can be set to yes or no");
    timeDerivativeSchemeOrder =
        ITHACAdict->lookupOrDefault<word>("timeDerivativeSchemeOrder", "second");
    M_Assert(timeDerivativeSchemeOrder == "first"
             || timeDerivativeSchemeOrder == "second",
             "The time derivative approximation must be set to either first or second order scheme in ITHACAdict");
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void psizeta::truthSolve(List<scalar> mu_now, fileName folder)
{
    Time& runTime = _runTime();
    dimensionedScalar& nu = _nu();
    surfaceScalarField& phi = _phi();
    fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    volScalarField& zeta = _zeta();
    volScalarField& psi = _psi();
    volVectorField& U = _U();
    volVectorField& psi_vec = _psi_vec();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;

    // Set time-dependent velocity BCs for initial condition
    /*if (timedepbcMethod == "yes")
    {
        for (label i = 0; i < inletPatch.rows(); i++)
        {
            Vector<double> inl(0, 0, 0);

            for (label j = 0; j < inl.size(); j++)
            {
                inl[j] = timeBCoff(i * inl.size() + j, 0);
            }

            assignBC(U, inletPatch(i, 0), inl);
        }
    }*/

    // Export and store the initial conditions for velocity and pressure
    // int counter = 3766;
    ITHACAstream::exportSolution(psi, name(counter), folder);
    ITHACAstream::exportSolution(zeta, name(counter), folder);
    ITHACAstream::exportSolution(U, name(counter), folder);
    std::ofstream of(folder + name(counter) + "/" +
                     runTime.timeName());
    zetafield.append(zeta.clone());
    psifield.append(psi.clone());
    psivecfield.append(U.clone());
    counter++;
    nextWrite += writeEvery;

    /*const label patch1 = mesh.boundaryMesh().findPatchID("top");
    const label patch2 = mesh.boundaryMesh().findPatchID("bottom");
    const label patch3 = mesh.boundaryMesh().findPatchID("left");
    const label patch4 = mesh.boundaryMesh().findPatchID("right");*/
    const label patch = mesh.boundaryMesh().findPatchID("sides");


    volScalarField ycord
(
    IOobject
    (
        "ycord",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField xcord
(
    IOobject
    (
        "xcord",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

    dimensionedScalar tempo2("tempo", dimTime, scalar(1.0));

    forAll(xcord,faceI)
    {           
       xcord[faceI] = mesh.C()[faceI].component(0);
    }  

    forAll(ycord,faceI)
    {           
       ycord[faceI] = mesh.C()[faceI].component(1);
    }  

    //cout << "provaxy" <<endl;
    //forcmodes[0] = 1; //ycord; //Foam::cos(3*ycord)*Foam::cos(3*xcord);
    //cout << "provaaaaaaa" <<endl;

    //ITHACAstream::exportSolution(forcmodes[0], "forcmodes", "forc");
    //cout << "XXXXXXXXXX" <<endl;

    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime);
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;

        // Set time-dependent velocity BCs
        /*if (timedepbcMethod == "yes")
        {
            for (label i = 0; i < inletPatch.rows(); i++)
            {
                Vector<double> inl(0, 0, 0);

                for (label j = 0; j < inl.size(); j++)
                {
                    inl[j] = timeBCoff(i * inl.size() + j, counter2);
                }

                assignBC(U, inletPatch(i, 0), inl);
            }

            counter2 ++;
        }*/


        // Vorticity equation

        fvScalarMatrix zetaEqn
        (
            fvm::ddt(zeta)
          + fvm::div(phi, zeta)
          - fvm::laplacian(nu, zeta) //+ 0.1*Foam::exp(-runTime.value()*nu.value())*Foam::sin(3.14*ycord)/sqr(tempo2) //0.1*Foam::exp(-runTime.value()*nu.value())*Foam::sin(3.14*ycord)/sqr(tempo2) //0.1*Foam::exp(-runTime.value()*nu.value())*
        );

        zetaEqn.solve();

        // Stream function equation

        fvScalarMatrix psiEqn
        (
            fvm::laplacian(dimensionedScalar("1", dimless, 1), psi) + zeta //+ ycord
        );

       psiEqn.solve();

       // update phi
        //volVectorField psi_vec("psi_vec", U*0);
        //psi_vec.component(0) = 0; //(o component(2))
        //psi_vec.component(1) = 0; //(o component(2))

        forAll(psi_vec,faceI)
           {           psi_vec[faceI].component(0) = 0; 
                       psi_vec[faceI].component(1) = 0; 
                       psi_vec[faceI].component(2) = psi[faceI]; 
           }  

forAll(psi_vec.boundaryField()[patch], faceI)
           {           psi_vec.boundaryFieldRef()[patch][faceI].component(0) = 0; 
                       psi_vec.boundaryFieldRef()[patch][faceI].component(1) = 0; 
                       psi_vec.boundaryFieldRef()[patch][faceI].component(2) = psi.boundaryField()[patch][faceI]; 
           }  


/*forAll(psi_vec.boundaryField()[patch2], faceI)
           {           psi_vec.boundaryFieldRef()[patch2][faceI].component(0) = 0; 
                       psi_vec.boundaryFieldRef()[patch2][faceI].component(1) = 0; 
                       psi_vec.boundaryFieldRef()[patch2][faceI].component(2) = psi.boundaryField()[patch2][faceI];
           }  

forAll(psi_vec.boundaryField()[patch3], faceI)
           {           psi_vec.boundaryFieldRef()[patch3][faceI].component(0) = 0; 
                       psi_vec.boundaryFieldRef()[patch3][faceI].component(1) = 0; 
                       psi_vec.boundaryFieldRef()[patch3][faceI].component(2) = psi.boundaryField()[patch3][faceI];
           }  

forAll(psi_vec.boundaryField()[patch4], faceI)
           {           psi_vec.boundaryFieldRef()[patch4][faceI].component(0) = 0; 
                       psi_vec.boundaryFieldRef()[patch4][faceI].component(1) = 0; 
                       psi_vec.boundaryFieldRef()[patch4][faceI].component(2) = psi.boundaryField()[patch4][faceI];
           } */
      

        U = fvc::curl(psi_vec);

                 forAll(U.boundaryField()[patch], faceI) //top
           {           
                     U.boundaryFieldRef()[patch][faceI].y() = 0; 
                     U.boundaryFieldRef()[patch][faceI].x() = 0; 
           } 

            /*     forAll(U.boundaryField()[patch2], faceI) //bottom
           {           
                     U.boundaryFieldRef()[patch2][faceI].y() = 0; 
                     //U.boundaryFieldRef()[patch2][faceI].x() = 0; 
           } 

                 forAll(U.boundaryField()[patch3], faceI) //left
           {           
                     //U.boundaryFieldRef()[patch3][faceI].y() = 0; 
                     U.boundaryFieldRef()[patch3][faceI].x() = 0; 
           } 

                 forAll(U.boundaryField()[patch4], faceI) //right
           {           
                     //U.boundaryFieldRef()[patch4][faceI].y() = 0; 
                     U.boundaryFieldRef()[patch4][faceI].x() = 0; 
           } */
        // --- Pressure-velocity PIMPLE corrector loop
        phi = fvc::interpolate(U) & mesh.Sf();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        if (checkWrite(runTime))
        {
            ITHACAstream::exportSolution(psi, name(counter), folder);
            ITHACAstream::exportSolution(zeta, name(counter), folder);
            ITHACAstream::exportSolution(U, name(counter), folder);
            std::ofstream of(folder + name(counter) + "/" +
                             runTime.timeName());
            psifield.append(psi.clone());
            zetafield.append(zeta.clone());
            psivecfield.append(U.clone());
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
            // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
            mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
            mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

            for (int i = 0; i < mu_now.size(); i++)
            {
                mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];
            }
        }
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == counter * mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   folder);
    }
}


bool psizeta::checkWrite(Time& timeObject)
{
    scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));
    scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) -
                          timeObject.deltaTValue());

    if ( diffnow < diffnext)
    {
        return true;
    }
    else
    {
        return false;
    }
}


void psizeta::project(fileName folder, label Npsi, label Nzeta, label Npsivec)
{
    Nzetamodes = Nzeta;
    Npsimodes = Npsi;
    Npsivecmodes = Npsivec;
    //NSUPmodes = 0;
    //L_U_SUPmodes.resize(0);

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        word Bpsi_str = "Bpsi_" + name(Npsimodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + Bpsi_str))
        {
            ITHACAstream::ReadDenseMatrix(BPsi_matrix, "./ITHACAoutput/Matrices/", Bpsi_str);
        }
        else
        {
            BPsi_matrix = diffusive_term_psi(Npsimodes);
        }

        word Bzeta_str = "Bzeta_" + name(Nzetamodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + Bzeta_str))
        {
            ITHACAstream::ReadDenseMatrix(BZeta_matrix, "./ITHACAoutput/Matrices/", Bzeta_str);
        }
        else
        {
            BZeta_matrix = diffusive_term_zeta(Nzetamodes);
        }

        word Zeta_str = "Zeta_" + name(Nzetamodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + Zeta_str))
        {
            ITHACAstream::ReadDenseMatrix(Zeta_matrix, "./ITHACAoutput/Matrices/", Zeta_str);
        }
        else
        {
           Zeta_matrix = mass_term_zeta(Nzetamodes);
        }

        word Mixed_str = "Mixed_" + name(Nzetamodes) + name(Npsimodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + Mixed_str))
        {
            ITHACAstream::ReadDenseMatrix(Mixed_matrix, "./ITHACAoutput/Matrices/", Mixed_str);
        }
        else
        {
            Mixed_matrix = hybrid_term(Nzetamodes, Npsimodes);
        }

        word G_str = "G_" + name(Nzetamodes) + name(Npsimodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + G_str))
        {
            ITHACAstream::ReadDenseTensor(zetaTensor, "./ITHACAoutput/Matrices/", G_str);
        }
        else
        {
            zetaTensor = divMomentum_zeta(Nzetamodes, Npsivecmodes);
        }

    }
    else
    {    
        BZeta_matrix = diffusive_term_zeta(Nzetamodes);
        BPsi_matrix = diffusive_term_psi(Npsimodes);
        Zeta_matrix = mass_term_zeta(Nzetamodes);
        Mixed_matrix = hybrid_term(Nzetamodes, Npsimodes);
        zetaTensor = divMomentum_zeta(Nzetamodes, Npsivecmodes);
        zetaTensor = divMomentum_zeta(Nzetamodes, Npsivecmodes);
        forcing_Matrix = forcingMatrix(Nzetamodes);
        /*if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }*/
    }

    // Export the matrices
    /*if (para->exportPython)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC4_matrix, "BC4", "python",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(gTensor, "G", "python", "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC4_matrix, "BC4", "matlab",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(gTensor, "G", "matlab", "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(D_matrix, "D", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC3_matrix, "BC3", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(BC4_matrix, "BC4", "eigen",
                                   "./ITHACAoutput/Matrices/");
        ITHACAstream::exportTensor(C_tensor, "C", "eigen",
                                   "./ITHACAoutput/Matrices/C");
        ITHACAstream::exportTensor(gTensor, "G", "eigen",
                                   "./ITHACAoutput/Matrices/G");
    }*/ //TODO finalizzare la function...
}


void psizeta::restart()
{
    volScalarField& psi = _psi();
    volVectorField& psi_vec = _psi_vec();
    volScalarField& psi0 = _psi0();
    volVectorField& psi_vec0 = _psi_vec0();
    volScalarField& zeta = _zeta();
    volScalarField& zeta0 = _zeta0();
    volVectorField& U = _U();
    volVectorField& U0 = _U0();
    surfaceScalarField& phi = _phi();
    surfaceScalarField& phi0 = _phi0();
    psi_vec = psi_vec0;
    U = U0;
    psi = psi0;
    zeta = zeta0;
    phi = phi0;
    //turbulence.reset(
    //    (incompressible::turbulenceModel::New(U, phi, _laminarTransport())).ptr()
    //);
}

Eigen::MatrixXd psizeta::diffusive_term_zeta(label Nzetamodes)
{
    label Bsize = Nzetamodes; // //
    //Info << "Bsize" << Bsize <<endl;
    Eigen::MatrixXd BZeta_matrix;
    BZeta_matrix.resize(Bsize, Bsize);

    // Project everything
    for (label i = 0; i < Bsize; i++)
    {   Info << "here" <<endl;
        for (label j = 0; j < Bsize; j++)
        {   Info << "here1" <<endl;
            BZeta_matrix(i, j) = fvc::domainIntegrate(zetamodes[i] * fvc::laplacian(
                    dimensionedScalar("1", dimless, 1),zetamodes[j])).value();
                    Info << "here2" <<endl;
        }
    }

    if (Pstream::parRun())
    {
        reduce(BZeta_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(BZeta_matrix, "./ITHACAoutput/Matrices/",
                                  "Bzeta_" + name(Npsimodes));
    return BZeta_matrix;
}

Eigen::MatrixXd psizeta::diffusive_term_psi(label Npsimodes)
{
    label Bsize = Npsimodes; // //
    Eigen::MatrixXd BPsi_matrix;
    BPsi_matrix.resize(Bsize, Bsize);

    // Project everything
    for (label i = 0; i < Bsize; i++)
    {
        for (label j = 0; j < Bsize; j++)
        {
            BPsi_matrix(i, j) = fvc::domainIntegrate(psimodes[i] * fvc::laplacian(
                    dimensionedScalar("1", dimless, 1),psimodes[j])).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(BPsi_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(BPsi_matrix, "./ITHACAoutput/Matrices/",
                                  "Psizeta_" + name(Npsimodes));
    return BPsi_matrix;
}

Eigen::MatrixXd psizeta::mass_term_zeta(label Nzetamodes)
{
    label Bsize = Nzetamodes; // //
    Eigen::MatrixXd Zeta_matrix;
    Zeta_matrix.resize(Bsize, Bsize);

    // Project everything
    for (label i = 0; i < Bsize; i++)
    {
        for (label j = 0; j < Bsize; j++)
        {
            Zeta_matrix(i, j) = fvc::domainIntegrate(zetamodes[i] * zetamodes[j]).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(Zeta_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(Zeta_matrix, "./ITHACAoutput/Matrices/",
                                  "Zeta_" + name(Nzetamodes));
    return Zeta_matrix;
}

Eigen::MatrixXd psizeta::hybrid_term(label Nzetamodes, label Npsimodes)
{
    label B1size = Npsimodes; // //
    label B2size = Nzetamodes; // //
    Eigen::MatrixXd Mixed_matrix;
    Mixed_matrix.resize(B1size, B2size);

    // Project everything
    for (label i = 0; i < B1size; i++)
    {
        for (label j = 0; j < B2size; j++)
        {
            Mixed_matrix(i, j) = fvc::domainIntegrate(psimodes[i] * zetamodes[j]).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(Mixed_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(Mixed_matrix, "./ITHACAoutput/Matrices/",
                                  "Mixed_" + name(Nzetamodes));
    return Mixed_matrix;
}

Eigen::Tensor<double, 3> psizeta::divMomentum_zeta(label Nzetamodes, label Npsivecmodes) //TODO da modificare ad hoc!!
{

    fvMesh& mesh = _mesh();
    const label patch = mesh.boundaryMesh().findPatchID("sides");

    //volVectorField psivecModes = psivecmodes[0] * 0;
    //volVectorField psivecModes2 = psivecModes;
    label g1Size = Nzetamodes;
    label g2Size = Npsivecmodes;
    label g3Size = Nzetamodes;
    Eigen::Tensor<double, 3> zetaTensor;
    zetaTensor.resize(g1Size, g2Size, g3Size);

    for (label i = 0; i < g1Size; i++)
    {
        for (label j = 0; j < g2Size; j++)
        {
            dimensionedScalar lung("lung", dimLength, scalar(1.0));
           volVectorField psivecModes ("psivecmodes",  psivecmodes[0]*0);
            //volVectorField psivecModes2 = psivecModes; //psivecmodes[0] * 0);
           volVectorField psivecModes2 ("psivecmodes",  psivecmodes[0]*0);
                   //psivecModes.component(2) = 0*psimodes[j];
        forAll (psivecModes, xx)
        {
             psivecModes[xx].component(2) = psimodes[j][xx];
                   //psivecModes.replace(2, psimodes[j]);
        }
                    psivecModes2 = fvc::curl(psivecModes)*lung;

                 forAll(psivecModes2.boundaryField()[patch], faceI) //top
           {           
                     psivecModes2.boundaryFieldRef()[patch][faceI].y() = 0; 
                     psivecModes2.boundaryFieldRef()[patch][faceI].x() = 0; 
           } 
            ITHACAstream::exportSolution(psivecModes2, name(j), "./ITHACAoutput/Matrices/");

            for (label k = 0; k < g3Size; k++)
            {
               // psivecmodes[j] = vector (0 0, psimodes[j]);
                zetaTensor(i, j, k) = fvc::domainIntegrate(zetamodes[i] * fvc::div(
                        linearInterpolate(psivecModes2) & psivecModes2.mesh().Sf(),
                        zetamodes[k])).value();

/*fvc::domainIntegrate(L_Uevolve_SUPmodes[i] & fvc::div(
                                        linearInterpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(),
                                        L_Uevolve_SUPmodes[k])).value();*/
                /*zetaTensor(i, j, k) = fvc::domainIntegrate(zetamodes[i] * fvc::div(
                                        linearInterpolate(psivecmodes[j]) & psivecmodes[j].mesh().Sf(),
                                        zetamodes[k])).value();*/
            }
        }
     }

    if (Pstream::parRun())
    {
        reduce(zetaTensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    // Export the tensor
    ITHACAstream::SaveDenseTensor(zetaTensor, "./ITHACAoutput/Matrices/",
                                  "zeta_Tensor_" +  name(Nzetamodes) + "_" + name(
                                      Npsivecmodes) + "_t");
    return zetaTensor;
}


        //zetaTensor = divMomentum_zeta(Nzetamodes, Npsimodes);


Eigen::MatrixXd psizeta::forcingMatrix(label Nzetamodes)
{
    label B1size = Nzetamodes; // //
    label B2size = 1; // //
    Eigen::MatrixXd forc_matrix;
    forc_matrix.resize(B1size, B2size);
    fvMesh& mesh = _mesh();
    Time& runTime = _runTime();

volScalarField ycord
(
    IOobject
    (
        "ycord",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField xcord
(
    IOobject
    (
        "xcord",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

    dimensionedScalar tempo2("tempo", dimTime, scalar(1.0));

    forAll(xcord,faceI)
    {           
       xcord[faceI] = mesh.C()[faceI].component(0);
    }  

    forAll(ycord,faceI)
    {           
       ycord[faceI] = mesh.C()[faceI].component(1);
    }  

    // Project everything
    for (label i = 0; i < B1size; i++)
    {
        for (label j = 0; j < B2size; j++)
        {
            forc_matrix(i, j) = fvc::domainIntegrate(zetamodes[i] * (Foam::sin(3.14*ycord))).value();//(Foam::cos(3*ycord)*Foam::cos(3*xcord))).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(forc_matrix, sumOp<Eigen::MatrixXd>());
    }

    ITHACAstream::SaveDenseMatrix(forc_matrix, "./ITHACAoutput/Matrices/",
                                  "forc_" + name(Nzetamodes));
    return forc_matrix;
}
