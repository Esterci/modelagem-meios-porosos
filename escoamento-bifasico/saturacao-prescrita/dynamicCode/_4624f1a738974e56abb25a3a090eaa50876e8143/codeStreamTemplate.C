/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 35 "/home/thiago/Documentos/gitHub/modelagem-meios-porosos/escoamento-bifasico/saturacao-prescrita/0/p/boundaryField/left/#codeStream"
#include "fvCFD.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    void codeStream_4624f1a738974e56abb25a3a090eaa50876e8143
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 46 "/home/thiago/Documentos/gitHub/modelagem-meios-porosos/escoamento-bifasico/saturacao-prescrita/0/p/boundaryField/left/#codeStream"
const IOdictionary& d = static_cast<const IOdictionary&> 
            ( 
                dict.parent().parent() 
            );

            const fvMesh& mesh = refCast<const fvMesh>(d.db()); 
            const label id = mesh.boundary().findPatchID("velocity-inlet-5"); 
            const fvPatch& patch = mesh.boundary()[id]; 
            vectorField p(patch.size(), 0);

            const scalar pi = constant::mathematical::pi; 
            const scalar U_0   = 2.; //maximum velocity
            const scalar p_ctr = 8.; //patch center 
            const scalar p_r   = 8.; //patch radius
            
            forAll(U, i) //equivalent to for (int i=0; patch.size()<i; i++)
            {
                const scalar y = patch.Cf()[i][1];
                U[i] = vector(U_0*(1-(pow(y -p_ctr,2))/(p_r*p_r)), 0., 0.);
            } 
            
            writeEntry(os, "", U);
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

