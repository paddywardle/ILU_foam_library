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
    ILU0After

Description
    ILU preconditioning with no fill-in

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    ILU0After.C

\*---------------------------------------------------------------------------*/

#ifndef ILU0After_H
#define ILU0After_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ILU0After Declaration
\*---------------------------------------------------------------------------*/

class ILU0After
:
    public lduPreconditioner
{
    // Private Data

        //- Preconditioned diagonal
        scalarField preconDiag_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ILU0After(const ILU0After&);

        //- Disallow default bitwise assignment
        void operator=(const ILU0After&);

        //- Precondition the diagonal
        void calcPreconDiag();


public:

    //- Runtime type information
    TypeName("ILU0After");


    // Constructors

        //- Construct from matrix and dictionary
        ILU0After
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& dict
        );

        //- Construct from matrix as a smoother
        ILU0After
        (
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& coupleBouCoeffs,
            const FieldField<Field, scalar>& coupleIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );


    //- Destructor
    virtual ~ILU0After();


    // Member Functions

        //- Execute preconditioning
        virtual void precondition
        (
            scalarField& wA,
            const scalarField& rA,
            const direction cmpt = 0
        ) const;

        //- Execute preconditioning with matrix transpose
        virtual void preconditionT
        (
            scalarField& wT,
            const scalarField& rT,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
