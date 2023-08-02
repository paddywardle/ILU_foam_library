/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    ILUC0After

Description
    ILU preconditioning without fill in based on Crout algorithm. L and U are
    calculated and stored.

    Reference: Saad, Y.: Iterative Methods for Sparse Linear Systems (2nd
    Edition), SIAM, 2003.

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved

SourceFiles
    ILUC0After.C

\*---------------------------------------------------------------------------*/

#ifndef ILUC0After_H
#define ILUC0After_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ILUC0After Declaration
\*---------------------------------------------------------------------------*/

class ILUC0After
:
    public lduPreconditioner
{
    // Private Data

        //- Preconditioned diagonal
        scalarField preconDiag_;

        //- Strictly lower part
        scalarField preconLower_;

        //- Strictly upper part
        scalarField preconUpper_;

        // Coupled Boundaries 
        FieldField<Field, scalar> preconCoupleBouCoeffs_;

        //- Temporary working diagonal
        scalar zDiag_;

        //- Temporary working row field
        scalarField z_;

        //- Temporary Working column field
        scalarField w_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ILUC0After(const ILUC0After&);

        //- Disallow default bitwise assignment
        void operator=(const ILUC0After&);

        //- Calculate LU factorization
        void calcFactorization();


public:

    //- Runtime type information
    TypeName("ILUC0After");


    // Constructors

        //- Construct from matrix and dictionary
        ILUC0After
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );

        //- Construct from matrix as a smoother
        ILUC0After
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );


    //- Destructor
    virtual ~ILUC0After();


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt = 0
        ) const;

        //- Execute preconditioning with matrix transpose
        virtual void preconditionT
        (
            scalarField& x,
            const scalarField& b,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
