
#include "unsteadyNSWPsi.H"

/// Source file of the unsteadyNSPsi class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
unsteadyNSWPsi::unsteadyNSWPsi() {}

// Construct from zero
unsteadyNSWPsi::unsteadyNSWPsi(int argc, char* argv[])
    :
    UnsteadyProblem()
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
    //supex = ITHACAutilities::check_sup();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void unsteadyNSWPsi::truthSolve(List<scalar> mu_now, fileName folder)
{
    cout<< "post createFields"<<endl; //da SPOSTARE IN BASSO FINCHE NON SMETTE DI9 STAMPARE
    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    //cout<< "siamo prima dei campi"<<endl;
    fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    //volScalarField& p = _p();
    
    volVectorField& U = _U(); 
    volScalarField& W = _W();
    volScalarField& Psi_z = _Psi_z();
    volVectorField& Psi = _Psi();
    volVectorField& temp = _temp();
    IOMRFZoneList& MRF = _MRF();
    //cout<< "siamo dopo i campi"<<endl;
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    // Perform a TruthSolve
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;

    // Set time-dependent velocity BCs for initial condition
    /* if (timedepbcMethod == "yes")
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
    } */

    // Export and store the initial conditions for velocity and pressure
    ITHACAstream::exportSolution(U, name(counter), folder);
    //ITHACAstream::exportSolution(p, name(counter), folder); 
    ITHACAstream::exportSolution(W, name(counter), folder);
    ITHACAstream::exportSolution(Psi_z, name(counter), folder);
    ITHACAstream::exportSolution(Psi, name(counter), folder);
    std::ofstream of(folder + name(counter) + "/" +
                     runTime.timeName());
    Ufield.append(U.clone());
    //Pfield.append(p.clone());
    Wfield.append(W.clone());
    Psi_zfield.append(Psi_z.clone());
    Psifield.append(Psi.clone());
    counter++;
    nextWrite += writeEvery;

    
    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime);
        runTime++;
        Info << "Time = " << runTime.timeName() << nl << endl;

        /* // Set time-dependent velocity BCs
        if (timedepbcMethod == "yes")
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
        } */

#include "W-PsiEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
/*         while (pimple.loop())
        {
#include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
#include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }  
 */
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        //Info<<"il valore di check runtime e':"<<checkWrite(runTime)<<"\n"<<endl;
        if (checkWrite(runTime))
        {
            //Info<< "PRE-EXPORT"<<endl;
            ITHACAstream::exportSolution(U, name(counter), folder);
            //ITHACAstream::exportSolution(p, name(counter), folder);
            ITHACAstream::exportSolution(W, name(counter), folder);
            ITHACAstream::exportSolution(Psi_z, name(counter), folder);
            ITHACAstream::exportSolution(Psi, name(counter), folder);
            Ufield.append(U.clone());
            //Info<<"U APPENDED"<<endl;
            //Pfield.append(p.clone());
            Wfield.append(W.clone());
            Psi_zfield.append(Psi_z.clone());
            Psifield.append(Psi.clone());
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
            // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
            mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
            mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

            for (label i = 0; i < mu_now.size(); i++)
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

