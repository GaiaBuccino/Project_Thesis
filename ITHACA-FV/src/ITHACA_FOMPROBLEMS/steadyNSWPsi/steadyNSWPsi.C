#include "steadyNSWPsi.H"
#include "steadyNS.H"
#include "viscosityModel.H"


steadyNSWPsi::steadyNSWPsi() : steadyNS() {}
steadyNSWPsi::steadyNSWPsi(int argc, char* argv[])  :
    steadyNS()
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
    _simple = autoPtr<simpleControl>
              (
                  new simpleControl
                  (
                      mesh
                  )
              );
    simpleControl& simple = _simple();
#include "createFields.H"
#include "createFvOptions.H"
    turbulence->validate();
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
    tolerance = ITHACAdict->lookupOrDefault<scalar>("tolerance", 1e-5);
    maxIter = ITHACAdict->lookupOrDefault<scalar>("maxIter", 1000);
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "lift");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty" || bcMethod == "none",
             "The BC method must be set to lift or penalty or none in ITHACAdict");
    fluxMethod = ITHACAdict->lookupOrDefault<word>("fluxMethod", "inconsistent");
    M_Assert(fluxMethod == "inconsistent" || bcMethod == "consistent",
             "The flux method must be set to inconsistent or consistent in ITHACAdict");
    para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    //supex = ITHACAutilities::check_sup();
};

// Method to perform a truthSolve
void steadyNSWPsi::truthSolve(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U(); 
    volScalarField& W = _W();
    volScalarField& Psi_z = _Psi_z();
    volVectorField& Psi = _Psi();
    volVectorField& temp = _temp();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    simpleControl& simple = _simple();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    //#include "NLsolvesteadyNS.H"
    ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/"); 
    ITHACAstream::exportSolution(W, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(Psi_z, name(counter), "./ITHACAoutput/Offline/");
    //Psi = Psi_z * temp;
    ITHACAstream::exportSolution(Psi, name(counter), "./ITHACAoutput/Offline/");
    Wfield.append(W.clone());
    Psi_zfield.append(Psi_z.clone());
    Psifield.append(Psi.clone());
    counter++;
    writeMu(mu_now);
    // --- Fill in the mu_samples with parameters (mu) to be used for the PODI sample points
    mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size());

    for (label i = 0; i < mu_now.size(); i++)
    {
        mu_samples(mu_samples.rows() - 1, i) = mu_now[i];
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   "./ITHACAoutput/Offline");
    }
}


///////////////////////////////
void steadyNSWPsi::projectSUP(fileName folder, label NW, label NPsi_z)
{
    NWmodes = NW;
    NPsi_zmodes = NPsi_z;
    NPsimodes = NPsi_z;

    //Info << "\nIN PROJECTSUP NWmodes = "<< NWmodes << endl;


    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {

        word AWPsi_str = "AWPsi_" + name(NWmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + AWPsi_str))
        {
            ITHACAstream::ReadDenseMatrix(AWPsi_matrix, "./ITHACAoutput/Matrices/", AWPsi_str);
        }
        else
        {
            AWPsi_matrix = diffusiveW_term(NWmodes, NPsi_zmodes);
        }

        word BWPsi_str = "BWPsi_" + name(NPsi_zmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + BWPsi_str))
        {
            ITHACAstream::ReadDenseMatrix(BWPsi_matrix, "./ITHACAoutput/Matrices/", BWPsi_str);
        }
        else
        {
            BWPsi_matrix = diffusivePsi_z_term(NWmodes, NPsi_zmodes);
        }

        word MW_str = "MW_" + name(NWmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + MW_str))
        {
            ITHACAstream::ReadDenseMatrix(MW_matrix, "./ITHACAoutput/Matrices/", MW_str);
        }
        else
        {
            MW_matrix = massW_term(NWmodes, NPsi_zmodes);
        }

        word MPsi_str = "MPsi_" + name(NPsi_zmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + MPsi_str))
        {
            ITHACAstream::ReadDenseMatrix(MPsi_matrix, "./ITHACAoutput/Matrices/", MPsi_str);
        }
        else
        {
            MPsi_matrix = massPsi_z_term(NWmodes, NPsi_zmodes);
        }

        word GWPsi_str = "GWPsi_" + name(NWmodes)+ "_" +name(NPsimodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + GWPsi_str))
        {
            ITHACAstream::ReadDenseTensor(GWPsi_tensor, "./ITHACAoutput/Matrices/", GWPsi_str);
        }
        else
        {
            GWPsi_tensor = convective_term_tens(NWmodes, NPsi_zmodes);
        }


       /*  if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }  */
    }
    else
    {
 
        AWPsi_matrix = diffusiveW_term(NWmodes, NPsi_zmodes);
        BWPsi_matrix = diffusivePsi_z_term(NWmodes, NPsi_zmodes);
        GWPsi_tensor = convective_term_tens(NWmodes, NPsi_zmodes);
        //KWPsi_matrix = pressure_gradient_term(NWmodes, NPsi_zmodes);
        //PWPsi_matrix = divergence_term(NWmodes, NPsi_zmodes);
        MW_matrix = massW_term(NWmodes, NPsi_zmodes);
        MPsi_matrix = massPsi_z_term(NWmodes, NPsi_zmodes);

        /* if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        } */
    }

}

Eigen::MatrixXd steadyNSWPsi::diffusiveW_term(label NWmodes, label NPsi_zmodes)
{
    Info << "\nNWmodes = " << NWmodes << endl;
    Info << "\nNPsi_zmodes = " << NPsi_zmodes << endl;
    label AWPsi_size = NWmodes;
    Info << "\nsize Wmodes = "<< Wmodes.size() << endl;
    Eigen::MatrixXd AWPsi_matrix;
    AWPsi_matrix.resize(AWPsi_size, AWPsi_size);


    for (label i = 0; i < AWPsi_size; i++)     //L_U_SUPmodes = modes W
    {   //Info << "\nsono in diffusive W ciclo esterno" << endl;
        for (label j = 0; j < AWPsi_size; j++)
        {
            //Info << "here1" <<endl;
            //Info << "\nMatrix A constructed " << NWmodes << endl;
            AWPsi_matrix(i, j) = fvc::domainIntegrate(Wmodes[i] * fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), Wmodes[j])).value();
            //Info << "\nsono in diffusive W ho fatto AWPSI" << endl;
        }
    }


    return AWPsi_matrix;
}


/// @brief 
/// @param NWmodes 
/// @param NPsi_zmodes 
/// @return 
Eigen::MatrixXd steadyNSWPsi::diffusivePsi_z_term(label NWmodes, label NPsi_zmodes)
{
    //Info << "\nsono in diffusive PSI " << endl;
    label BWPsi_size = NPsi_zmodes;
    Eigen::MatrixXd BWPsi_matrix;
    BWPsi_matrix.resize(BWPsi_size, BWPsi_size);

    for (size_t i = 0; i < Psi_zmodes.size(); i++)
    {   
        //Info<<"\n Psi_zmodes size =" << Psi_zmodes.size()<<endl;
        //Info<<"\n Psi_zmodes[i] size =" << Psi_zmodes[i].size()<<endl;
        //Info<<"\n value of i =" << i <<endl;

        for (size_t j = 0; j < Psi_zmodes[i].size(); j++)
        { 
            //Info << "\ncalcolo Psi field vettoriale"<< endl;
            //Info<<"\n value of j =" << j <<endl;

            Psimodes[i][j][2] = Psi_zmodes[i][j];
            //Psimodes[i][2] = Psi_zmodes[i];
            //Info<<"\n valore test di Psi =" << Psimodes[0][0][2]<<endl;

        }
        
    } 
    

    for (label i = 0; i < BWPsi_size; i++)     //L_U_SUPmodes = modes W
    {   //Info << "\nsono in diffusive PSI PRIMO CICLO " << endl;
        for (label j = 0; j < BWPsi_size; j++)
        {
            //Info << "\nsono in diffusive PSI SECONDO CICLO " << endl;
            BWPsi_matrix(i, j) = fvc::domainIntegrate(Psi_zmodes[i] * fvc::laplacian(
                    dimensionedScalar("1", dimless, 1), Psi_zmodes[j])).value();
            //Info << "\nB constructed " << endl;
        }
    }

    return BWPsi_matrix;
}


Eigen::Tensor<double, 3> steadyNSWPsi::convective_term_tens(label NWmodes,
        label NPsi_zmodes)

{ //creare psi_modes
    label CW_size = NWmodes;
    label CPsi_size = NPsi_zmodes;
    Eigen::Tensor<double, 3> C_tensor;
    C_tensor.resize(CW_size, CPsi_size, CW_size);
    //const fvMesh& mesh = L_U_SUPmodes[0].mesh();   

    /* for (size_t i = 0; i < Psi_zmodes.size(); i++)
    {
        Info<<"\n Psi_zmodes[i] size =" << Psi_zmodes[i].size()<<endl;

        for (size_t j = 0; j < Psi_zmodes[i].size(); j++)
        { 
            Info << "\ncalcolo Psi field vettoriale"<< endl;
            Psimodes[i][j][2] = Psi_zmodes[i][j];
            //Psimodes[i][2] = Psi_zmodes[i];
            //Info<<"\n valore test di Psi =" << Psimodes[0][0][2]<<endl;

        }
        
    } */
   
    for (label i = 0; i < CW_size; i++)
    {
        //Info << "\n sono nel primo ciclo per GWPsi "<<endl;
        for (label j = 0; j < CPsi_size; j++)
        {
            
            //Psimodes[j] = Psi_zmodes[j]*temp[j];   
            for (label k = 0; k < CW_size; k++)
            {
                //Info << "\n sono nel terzo ciclo per GWPsi "<<endl;
                
                //if (fluxMethod == "consistent")
                //{
                    
                    volVectorField curl_Psi = fvc::curl(Psimodes[j]);
                    C_tensor(i, j, k) = fvc::domainIntegrate(Wmodes[i] * fvc::div(
                                            fvc::flux(curl_Psi),
                                            Wmodes[k])).value();
                //}
                //Info << "\nC constructed" << endl;
                /* else
                {
                    C_tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::div(
                                            linearInterpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(),
                                            L_U_SUPmodes[k])).value();
                } */
            }
        }
    }

 
    return C_tensor;
}


Eigen::MatrixXd steadyNSWPsi::massW_term(label NWmodes, label NPsi_zmodes)
//modified
{
    label MW_size = NWmodes;
    Eigen::MatrixXd MW_matrix(MW_size, MW_size);

    // Project everything
    for (label i = 0; i < MW_size; i++)
    {
        for (label j = 0; j < MW_size; j++)
        {
            MW_matrix(i, j) = fvc::domainIntegrate(Wmodes[i] *
                                                  Wmodes[j]).value();
            //Info << "\nMW constructed " << endl;
        }
    }

     
    return MW_matrix;
}


Eigen::MatrixXd steadyNSWPsi::massPsi_z_term(label NWmodes, label NPsi_zmodes)
{
    label MW_size = NWmodes;
    label MPsi_size = NPsi_zmodes;
    Eigen::MatrixXd MPsi_matrix(MPsi_size, MW_size);

    // Project everything
    for (label i = 0; i < MPsi_size; i++)
    {
        //Psimodes[i] = Psi_zmodes[i]*temp[i];
        for (label j = 0; j < MW_size; j++)
        {
            MPsi_matrix(i, j) = fvc::domainIntegrate(Psi_zmodes[i] * 
                                                   Wmodes[j]).value();
            //Info << "\nMPsi constructed " << endl;
        }
    }

  

    return MPsi_matrix;
}
