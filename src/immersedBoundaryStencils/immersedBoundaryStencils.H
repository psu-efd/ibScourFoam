 /*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
     \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
      \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

 Class
     Foam::immersedBoundaryStencils

 Description
     Provide stencil and stencil weights used for immersed boundary method interpolation


 SourceFiles
     immersedBoundaryStencils.C

 \*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryStencils_H
#define immersedBoundaryStencils_H


#include "DynamicList.H"
#include "meshSearch.H"
#include "SubField.H"
#include "volFields.H"
#include "fvCFD.H"
#include "IOmanip.H"
#include "vectorList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class immersedBoundaryStencils Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryStencils
{
typedef List<vectorList> vectorListList;

protected:
// Protected data

    //- Reference to the mesh
    const fvMesh& mesh_;

    //- Reference to the mesh
    const meshSearch& queryMesh_;

    //- Reference to the gamma
    const volScalarField& gamma_;

    //- Reference to the sampling points
    const pointField& points_;

    //- Reference to the the cell for each sampling points
    labelList cellIDList_;

private:
// Private member data

    //- name
    mutable word* namePtr_;

    mutable labelListList* cellCellsPtr_;

    mutable vectorListList* procCentresPtr_;

    mutable scalarListList* procGammaPtr_;

    mutable List<List<labelPair>>* cellProcCellsPtr_;

    mutable labelListList* procCellsPtr_;

    mutable scalarListList* cellWeightsPtr_;

    mutable scalarListList* cellProcWeightsPtr_;

    mutable dictionary* dictPtr_;

    mutable scalar maxCellCellRows_;

    mutable scalar radiusFactor_;

    mutable scalar gammaInclude_c_;

    mutable scalar gammaWeight_c_;

    //- Disallow default bitwise copy construct
    immersedBoundaryStencils(const immersedBoundaryStencils&);

    //- Disallow default bitwise assignment
    void operator=(const immersedBoundaryStencils&);

    void makeStencils() const;
    void makeStencilsWeights() const;

    void findCellCells
    (
        const vector& pt,
        const label cellID,
        labelList& cellCells
    )const;

public:

    //- Runtime type information
    TypeName("immersedBoundaryStencils");

 // Constructors

    //- Construct from triSurface
    immersedBoundaryStencils
    (
        const word& name,
        const pointField& points,
        const labelList& cellIDList,
        const volScalarField& gamma,
        const meshSearch& queryMesh,
        const dictionary& dict
    );

    scalar cellSize(label cellID) const;

//- Destructor
virtual ~immersedBoundaryStencils();

    //- Clear storage
    void clearOut();


// Member Functions
    //- name
    const word& name() const
    {
        return *namePtr_;
    }

    const fvMesh& mesh() const
    {
        return mesh_;
    }

    const meshSearch& queryMesh() const
    {
        return queryMesh_;
    }

    const pointField& points() const
    {
        return points_;
    }

    const labelList& cellIDList() const
    {
        return cellIDList_;
    }

    const volScalarField& gamma() const
    {
        return gamma_;
    }

    scalar& maxCellCellRows() const
    {
        return maxCellCellRows_;
    }

    scalar& radiusFactor() const
    {
        return radiusFactor_;
    }

    scalar& gammaInclude_c() const
    {
        return gammaInclude_c_;
    }

    scalar& gammaWeight_c() const
    {
        return gammaWeight_c_;
    }

    const labelListList& cellCells() const
    {
        if(!cellCellsPtr_)
        {
            makeStencils();
        }
        return *cellCellsPtr_;
    };

    const vectorListList& procCentres() const
    {
        if(!procCentresPtr_)
        {
            makeStencils();
        }
        return *procCentresPtr_;
    };

    const scalarListList& procGamma() const
    {
        if(!procGammaPtr_)
        {
            makeStencils();
        }
        return *procGammaPtr_;
    };

    const List<List<labelPair>>& cellProcCells() const
    {
        if(!cellProcCellsPtr_)
        {
            makeStencils();
        }
        return *cellProcCellsPtr_;
    };

    const labelListList& procCells() const
    {
        if(!procCellsPtr_)
        {
            makeStencils();
        }
        return *procCellsPtr_;
    };

    const scalarListList& weights() const
    {
        if(!cellWeightsPtr_)
        {
            makeStencilsWeights();
        }
        return *cellWeightsPtr_;
    };

    const scalarListList& procWeights() const
    {
        if(!cellProcWeightsPtr_)
        {
            makeStencilsWeights();
        }
        return *cellProcWeightsPtr_;
    };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
