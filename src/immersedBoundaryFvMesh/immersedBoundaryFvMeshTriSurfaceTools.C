/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

SourceFiles
    immersedBoundaryFvMeshTriSurfaceTools.C
\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvMesh.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::immersedBoundaryFvMesh::writeVTK
(
    const fileName &surfaceName,
    const word &fieldName,
    const Field<Type>& psi,
    label objectID
) const
{
    const triSurface& ts = objectsList()[objectID];
    if (psi.size() != ts.size())
    {
        word names = objectNames(objectID);
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::Field<Type>\n"
            "immersedBoundaryFvMesh::writeVTK\n"
            "(\n"
            "    const Field<Type>& ibValues,\n"
            "    label objectID\n"
            ") const"
        )   << "Field size does not correspond to size of "
            << "triangulated surface for object " << objectNames(objectID) << nl
            << "Field size = " << psi.size()
            << " triSurface size = " << ts.size()
            << abort(FatalError);
    }

    fileName path;
    if (Pstream::parRun())
    {
        path = this->time().path()/".."/"postProcessing"/"VTK";
    }
    else
    {
        path = this->time().path()/"postProcessing"/"VTK";
    }

    autoPtr <surfaceWriter> writerPtr =
        surfaceWriter::New("vtk");

    faceList f(ts.size());
    forAll (ts, faceI)
    {
        f[faceI] = ts[faceI].triFaceFace();
    }
    fileName outputDir = path/fieldName;

    mkDir(outputDir);

    writerPtr->write
    (
        outputDir,
        surfaceName+"_"+ this->time().timeName(),
        ts.points(),
        f,
        fieldName,
        psi,
        false
    );
}
// ************************************************************************* //

