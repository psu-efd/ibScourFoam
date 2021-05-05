/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"


#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef IncompressibleTurbulenceModel<transportModel>
                     transportModelIncompressibleTurbulenceModel;
    typedef RASModel<transportModelIncompressibleTurbulenceModel>
                  RAStransportModelIncompressibleTurbulenceModel;
    typedef LESModel<transportModelIncompressibleTurbulenceModel>
                  LEStransportModelIncompressibleTurbulenceModel;
}

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (transportModelIncompressibleTurbulenceModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (transportModelIncompressibleTurbulenceModel, LES, Type)

#include "kEpsilonIB/kEpsilonIB.H"
makeRASModel(kEpsilonIB);

#include "kOmegaIB/kOmegaIB.H"
makeRASModel(kOmegaIB);

//#include "Base/kOmegaSSTBaseIB/kOmegaSSTBaseIB.H"
//makeRASModel(kOmegaSSTBaseIB);

//#include "kOmegaSSTIB/kOmegaSSTIB.H"
//makeRASModel(kOmegaSSTIB);

#include "kOmegaSSTLMIB/kOmegaSSTLMIB.H"
makeRASModel(kOmegaSSTLMIB);

#include "kOmegaSSTSASIB/kOmegaSSTSASIB.H"
makeRASModel(kOmegaSSTSASIB);

#include "SmagorinskyIB/SmagorinskyIB.H"
makeLESModel(SmagorinskyIB);

#include "dynamicKEqnIB/dynamicKEqnIB.H"
makeLESModel(dynamicKEqnIB);

#include "WALEIB/WALEIB.H"
makeLESModel(WALEIB);
// ************************************************************************* //
