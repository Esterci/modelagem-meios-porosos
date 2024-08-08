/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "simpleControl.H"

#include "RelPerm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"


    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info<< "\nCalculating p and Sw field\n" << endl;

        fvModels.correct();

        while (simple.correctNonOrthogonal())
        {    
            Info<< "\n->> Sw field\n" << endl;
            // printField(Sw);

            krw = correct_krw(Sw);

            // Info<< "\n->> krw field\n" << endl;
            // printField(krw);

            kro = correct_kro(Sw);
            
            // Info<< "\n->> kro field\n" << endl;
            // printField(kro);

            alpha = myScalar*(-K * (kro/mu_o + krw/mu_w));
            
            // Info<< "\n->> alpha field\n" << endl;
            // printField(alpha);

            // Solve pressure equation
            fvScalarMatrix pEqn
            (
                -fvm::laplacian(alpha, p)
                ==
                fvModels.source(p)
            );
            // fvConstraints.constrain(pEqn);
            pEqn.solve();
            // fvConstraints.constrain(p);

            ut = (alpha)*fvc::grad(p);

            phi = fvc::flux(ut);

            

            // // // Solve water saturation equation
            // // dimensionedScalar myScalar2("myScalar2", dimensionSet(0,0,1,0,0,0,0), 1.0);
            // // dimensionedScalar myScalar3("myScalar3", dimensionSet(0,2,0,0,0,0,0), 1.0);
            
            fw = correct_fw(krw,kro);

            fvScalarMatrix SwEqn
            (
                prst*fvm::ddt(Sw) + fvc::div(phi,fw) // 
                // 0.15*fvm::ddt(Sw) + fvc::div(ut) // 
                ==
                fvModels.source(Sw)
            );
            // fvConstraints.constrain(SwEqn);
            SwEqn.solve();
            // fvConstraints.constrain(Sw);
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
