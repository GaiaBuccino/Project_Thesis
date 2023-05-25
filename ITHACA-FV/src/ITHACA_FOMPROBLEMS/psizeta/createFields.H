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

Info << "Reading thermophysical properties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Laminar viscocity [m2/s]
_nu = autoPtr<dimensionedScalar>
      (
          new dimensionedScalar
          (
              "nu",
              dimViscosity,
              transportProperties.lookup("nu")
          )
      );
dimensionedScalar& nu = _nu();

Info << "Reading field zeta\n" << endl;
_zeta = autoPtr<volScalarField>
     (
         new volScalarField
         (
             IOobject
             (
                 "zeta",
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh
         )
     );
volScalarField& zeta = _zeta();

volScalarField zeta0(zeta);

_zeta0 = autoPtr<volScalarField>
      (
          new volScalarField(zeta0)
      );

Info << "Reading field psi\n" << endl;
_psi= autoPtr<volScalarField>
     (
         new volScalarField
         (
             IOobject
             (
                 "psi",
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh
         )
     );
volScalarField& psi = _psi();

volScalarField psi0(psi);

_psi0 = autoPtr<volScalarField>
      (
          new volScalarField(psi0)
      );

Info << "Reading field psivec\n" << endl;
_psi_vec = autoPtr<volVectorField>
     (
         new volVectorField
         (
             IOobject
             (
                 "psi_vec",
                 runTime.timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh
         )
     );
volVectorField& psi_vec = _psi_vec();

volVectorField psi_vec0(psi_vec);

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

_psi_vec0 = autoPtr<volVectorField>
      (
          new volVectorField(psi_vec0)
      );

Info << "Reading field U\n" << endl;
_U = autoPtr<volVectorField>
     (
         new volVectorField
         (
             IOobject
             (
                 "U",
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh
         )
     );
volVectorField& U = _U();

volVectorField U0(U);

_U0 = autoPtr<volVectorField>
      (
          new volVectorField(U0)
      );

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
/*setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);*/

_laminarTransport = autoPtr<singlePhaseTransportModel>
                    (
                        new singlePhaseTransportModel( U, phi )
                    );
singlePhaseTransportModel& laminarTransport = _laminarTransport();

turbulence = autoPtr<incompressible::turbulenceModel>
             (
                 incompressible::turbulenceModel::New(U, phi, laminarTransport)
             );

_MRF = autoPtr<IOMRFZoneList>
       (
           new IOMRFZoneList(mesh)
       );
