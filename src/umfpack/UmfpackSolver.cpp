
/*
 *  Copyright (C) 2016 DRG. All rights reserved.
 *
 *  Erik Jan Lingen
 *
 *  Copyright (C) 2017 TU Delft. All rights reserved.
 *  
 *  Frank Everdij, Januari 2017
 *  
 *  - Fixed transpose bug, since JemJive handles Sparse Matrices in
 *    row format (CRS) and not column format (CCS), like Umfpack does
 *  - Changed array pointers and Umfpack calling functions type from integer
 *    to long
 *  - Added event functions
 *  - Added debugging output
 *  - Removed storing of symbolic_ data structure: this clears up memory at
 *    the cost of having to recompute factorization.
 *
 */

#include <cmath>
#include <fstream>
#include <jem/base/assert.h>
#include <jem/base/limits.h>
#include <jem/base/Array.h>
#include <jem/base/Float.h>
#include <jem/base/System.h>
#include <jem/base/ArithmeticException.h>
#include <jem/base/IllegalArgumentException.h>
#include <jem/util/Event.h>
#include <jive/util/utilities.h>
#include <jive/util/Constraints.h>
#include <jive/algebra/SparseMatrixObject.h>
#include <jive/solver/SolverParams.h>
#include <jive/solver/SolverFactory.h>
#include <jive/solver/SolverException.h>
#include <jive/solver/StdConstrainer.h>
#include <jive/solver/DummyConstrainer.h>
#include "UmfpackSolver.h"


using jem::newInstance;
using jive::algebra::SparseMatrixExt;

//=======================================================================
//   class UmfpackSolver
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  UmfpackSolver::TYPE_NAME = "Umfpack";
const int    UmfpackSolver::MAX_ITER  = 10;

const int    UmfpackSolver::NEW_VALUES_        = 1 << 0;
const int    UmfpackSolver::NEW_STRUCT_        = 1 << 1;



//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


UmfpackSolver::UmfpackSolver

  ( const String&        name,
    Ref<AbstractMatrix>  matrix,
    Ref<Constraints>     cons ) :

    Super ( name )

{
  using jem::util::connect;
  using jive::util::joinNames;
  using jive::solver::StdConstrainer;
  using jive::solver::DummyConstrainer;

  jem::System::debug ( myName_ )  << myName_ << "> Umfpack constructor\n";

  JEM_PRECHECK ( matrix != NIL );

  String  conName = joinNames ( myName_, "constrainer" );


  // Check whether the matrix is compatible with this solver.

  if ( ! matrix->hasExtension<SparseMatrixExt>() )
  {
    throw jem::IllegalArgumentException (
      JEM_FUNC,
      matrix->getContext() +
      " does not implement the sparse matrix extension"
    );
  }

  if ( matrix->isDistributed() )
  {
    throw jem::IllegalInputException (
      getContext (),
      getContext () + " does not support distributed matrices"
    );
  }

  // Create a Constrainer that applies the constraints to the
  // system matrix.

  if ( cons == NIL )
  {
    conman_ =

      newInstance<DummyConstrainer> ( conName, matrix );
  }
  else
  {
    conman_ =

      newInstance<StdConstrainer>   ( conName, cons, matrix );
  }

  matrix_    = conman_->getOutputMatrix ();
  numeric_   = NULL;
  mode_      = 0;
  precision_ = PRECISION;
  events_    = ~0x0;

  // Set default umfpack parameters

  umfpack_dl_defaults( control_ );

  // Delete the symbolic/numeric factorization if the matrix is
  // modified.

  connect ( matrix_->newValuesEvent, this, & Self::valuesChanged_ );
  connect ( matrix_->newStructEvent, this, & Self::structChanged_ );
}


UmfpackSolver::~UmfpackSolver ()
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack destructor\n";

  freeNumeric_ ();
}


//-----------------------------------------------------------------------
//   improve  
//-----------------------------------------------------------------------


void UmfpackSolver::improve

  ( const Vector&  lhs,
    const Vector&  rhs )

{
  using jem::dot;
  using jem::max;
  using jem::Float;
  using jem::ArithmeticException;
  using jive::solver::SolverException;

  jem::System::debug ( myName_ )  << myName_ << "> Umfpack improve\n";

  JEM_PRECHECK ( lhs.size() == rhs.size() );

  Matrix  buf;
  Vector  du;
  Vector  r;
  Vector  u;
  Vector  f;

  double  rscale;
  int     umferr;

  // Special event call declaration during solve

  jem::util::Event<double,Self&>       solveEvent;

  // Apply the constraints and update the symbolic/numeric
  // factorization if necessary.

  conman_->update ();

  if (events_)
  {
    update_ ();
  }

  const idx_t  dofCount =  matrix_->size(0);

  if ( lhs.size() != dofCount )
  {
    jive::util::sizeError ( getContext(),
                      "lhs vector", lhs.size(), dofCount );
  }

  // Allocate local vectors.

  if ( mode_ & PRECON_MODE )
  {
    Matrix  vbuf ( dofCount, 2 );

    r .ref ( vbuf[0] );
    du.ref ( vbuf[1] );
    u .ref ( lhs );
    f .ref ( rhs );
  }
  else
  {
    Matrix  vbuf ( dofCount, 4 );

    r .ref ( vbuf[0] );
    du.ref ( vbuf[1] );
    u .ref ( vbuf[2] );
    f .ref ( vbuf[3] );

    conman_->getRhs  ( f, rhs );
    conman_->initLhs ( u, lhs );
  }

  iiter_ = 0;
  error_  = 0.0;
  rscale = std::sqrt ( dot( f, f ) );

  if ( Float::isNaN( rscale ) )
  {
    throw ArithmeticException (
      getContext (),
      "invalid norm of right-hand vector: NaN"
    );
  }

  // Nothing to do if the right-hand side vector is zero.

  if ( Float::isTiny( rscale ) )
  {
    u = 0.0;

    if ( ! (mode_ & PRECON_MODE) )
    {
      conman_->getLhs ( lhs, u );
    }

    return;
  }

  // Compute the initial residual vector.

  rscale = 1.0 / rscale;

  matrix_->matmul ( r, u );

  r = f - r;

  for ( iiter_ = 0; iiter_ < MAX_ITER; iiter_++ )
  {
    // Compute the new soluton and residual vector.

    umferr = umfpack_dl_solve ( UMFPACK_At,
                                offsets_.addr (),
                                indices_.addr (),
                                values_ .addr (),
                                du      .addr (),
                                r       .addr (),
                                numeric_,
                                control_,
                                info_ );
    if ( umferr < 0 )
    {
      umfpack_dl_report_info(control_, info_);
      umfpack_dl_report_status(control_, umferr);

//      throw ArithmeticException (
//        getContext (),
//        String::format ("error in umfpack_di_solve: %d", umferr )
//      );
    }

    u += du;

    matrix_->matmul ( r, u );
    
    r     = f - r;

    error_ = rscale * std::sqrt ( dot( r, r ) );

    solveEvent.emit ( error_, *this );

    if ( Float::isNaN( error_) )
    {
      throw ArithmeticException (
        getContext (),
        "invalid norm of residual vector: NaN"
      );
    }

    if ( error_ <= precision_ || error_ > 1.0e5 )
    {
      break;
    }
  }

  // Check if the solution has been obtained.

  if ( (error_ > max( 1.0, precision_ )) ||
       (error_ > precision_ && ! (mode_ & LENIENT_MODE)) )
  {
    throw SolverException (
      getContext     (),
      String::format ( "residual norm too large: %e", error_ )
    );
  }

  if ( ! (mode_ & PRECON_MODE) )
  {
    conman_->getLhs ( lhs, u );
  }
}


//-----------------------------------------------------------------------
//   setMode
//-----------------------------------------------------------------------


void UmfpackSolver::setMode ( int mode )
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack setMode\n";

  mode_ = mode;
}


//-----------------------------------------------------------------------
//   getMode
//-----------------------------------------------------------------------


int UmfpackSolver::getMode () const
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack getMode\n";

  return mode_;
}


//-----------------------------------------------------------------------
//   setPrecision
//-----------------------------------------------------------------------


void UmfpackSolver::setPrecision ( double eps )
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack setPrecision\n";

  JEM_PRECHECK ( eps >= 0.0 );

  precision_ = eps;
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


double UmfpackSolver::getPrecision () const
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack getPrecision\n";

  return precision_;
}


//-----------------------------------------------------------------------
//   getMatrix
//-----------------------------------------------------------------------


AbstractMatrix* UmfpackSolver::getMatrix () const
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack getMatrix\n";

  return conman_->getInputMatrix ();
}


//-----------------------------------------------------------------------
//   getConstraints
//-----------------------------------------------------------------------


Constraints* UmfpackSolver::getConstraints () const
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack getConstraints\n";

  return conman_->getConstraints ();
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Solver> UmfpackSolver::makeNew

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::solver::SolverParams;

  jem::System::debug () << "> Umfpack makeNew\n";

  Ref<AbstractMatrix>  mat;
  Ref<Constraints>     cons;


  params.find ( mat,  SolverParams::MATRIX      );
  params.find ( cons, SolverParams::CONSTRAINTS );

  if ( mat != NIL && mat->hasExtension<SparseMatrixExt>() )
  {
    return jem::newInstance<Self> ( name, mat, cons );
  }
  else
  {
    return NIL;
  }
}


//-----------------------------------------------------------------------
//   declare
//-----------------------------------------------------------------------


void UmfpackSolver::declare ()
{
  using jive::solver::SolverFactory;

  jem::System::debug () << "> Umfpack declare\n";

  SolverFactory::declare ( TYPE_NAME,  & makeNew );
  SolverFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   update_
//-----------------------------------------------------------------------


void UmfpackSolver::update_ ()
{
  using jem::maxOf;
  using jem::castTo;
  using jem::System;
  using jem::ArithmeticException;
  using jive::SparseMatrix;
  using jem::slice;

  jem::System::debug ( myName_ )  << myName_ << "> Umfpack update_\n";

  SparseMatrixExt*  sx  =

    matrix_->getExtension<SparseMatrixExt>  ();

  SparseMatrix      sm = sx->toSparseMatrix ();

  const idx_t       n   = sm.size (0);
  const idx_t       nnz = sm.nonZeroCount (); // Number of non-zeroes.
  int               umferr;
  Vector            accu ( n );


  JEM_PRECHECK ( sm.size(0) == sm.size(1) );
  JEM_PRECHECK ( n <= maxOf<int>() );

  offsets_.resize ( n + 1 );  // Ap
  indices_.resize ( nnz );    // Ai
  values_ .resize ( nnz );    // Ax

  offsets_ = castTo<int> ( sm.getRowOffsets() );
  indices_ = castTo<int> ( sm.getColumnIndices() );
  values_  =               sm.getValues ();



  if ( ! numeric_ )
  {
  // Sort the matrix entries according to their column index.

  for ( idx_t irow = 0; irow < n; irow++ )
  {
    int  i = offsets_[irow];
    int  k = offsets_[irow + 1];

    for ( int j = i; j < k; j++ )
    {
      accu[indices_[j]] = values_[j];
    }

    sort ( indices_[slice(i,k)] );

    for ( int j = i; j < k; j++ )
    {
      values_[j] = accu[indices_[j]];
    }
  }

    umferr = umfpack_dl_symbolic ( (int) n, (int) n,
                                   offsets_.addr (),
                                   indices_.addr (),
                                   values_ .addr (),
                                   &symbolic_,
                                   control_,
                                   info_ );

    if ( umferr < 0 )
    {
      umfpack_dl_report_info(control_, info_);
      umfpack_dl_report_status(control_, umferr);

//      throw ArithmeticException (
//        getContext (),
//        String::format ("error in umfpack_di_symbolic: %d", umferr )
//      );
    }

    umferr = umfpack_dl_numeric ( offsets_.addr (),
                                  indices_.addr (),
                                  values_ .addr (),
                                  symbolic_,
                                  &numeric_,
                                  control_,
                                  info_ );

    if ( umferr < 0 )
    {
      umfpack_dl_report_info(control_, info_);
      umfpack_dl_report_status(control_, umferr);

//      throw ArithmeticException (
//        getContext (),
//        String::format ("error in umfpack_di_numeric: %d", umferr )
//      );
    }
    freeSymbolic_();

  }
//  factorEvent.emit ( 100, *this );

  matrix_->resetEvents ();

  events_ = 0;

}


//-----------------------------------------------------------------------
//   valuesChanged_
//-----------------------------------------------------------------------


void UmfpackSolver::valuesChanged_ ()
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack valuesChanged_\n";

  events_ |= NEW_VALUES_;

  UmfpackSolver::freeNumeric_();
}


//-----------------------------------------------------------------------
//   structChanged_
//-----------------------------------------------------------------------


void UmfpackSolver::structChanged_ ()
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack structChanged_\n";

  events_ |= NEW_STRUCT_;

  UmfpackSolver::freeNumeric_();
}


//-----------------------------------------------------------------------
//   freeSymbolic_
//-----------------------------------------------------------------------


void UmfpackSolver::freeSymbolic_ ()
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack freeSymbolic\n";

  if ( symbolic_ )
  {
    umfpack_dl_free_symbolic ( &symbolic_ );

    symbolic_ = NULL;
  }
}


//-----------------------------------------------------------------------
//   freeNumeric_
//-----------------------------------------------------------------------


void UmfpackSolver::freeNumeric_ ()
{
  jem::System::debug ( myName_ )  << myName_ << "> Umfpack freeNumeric\n";

  if ( numeric_ )
  {
    umfpack_dl_free_numeric ( &numeric_ );

    numeric_ = NULL;
  }
}

//-----------------------------------------------------------------------
//   dump_
//-----------------------------------------------------------------------

void UmfpackSolver::dump_ ()
{
  using jem::maxOf;
  using jem::castTo;
  using jem::System;
  using jem::ArithmeticException;
  using jive::SparseMatrix;
  using jem::slice;

  jem::System::debug ( myName_ )  << myName_ << "> Umfpack dump_\n";

  SparseMatrixExt*  sx  =

    matrix_->getExtension<SparseMatrixExt>  ();

  SparseMatrix      sm = sx->toSparseMatrix ();

  const idx_t       n   = sm.size (0);
  const idx_t       nnz = sm.nonZeroCount (); // Number of non-zeroes.
  Vector            accu ( n );
  std::ofstream f;


  JEM_PRECHECK ( sm.size(0) == sm.size(1) );
  JEM_PRECHECK ( n <= maxOf<int>() );

  offsets_.resize ( n + 1 );  // Ap
  indices_.resize ( nnz );    // Ai
  values_ .resize ( nnz );    // Ax

  offsets_ = castTo<int> ( sm.getRowOffsets() );
  indices_ = castTo<int> ( sm.getColumnIndices() );
  values_  =               sm.getValues ();



  // Sort the matrix entries according to their column index.

  for ( idx_t irow = 0; irow < n; irow++ )
  {
    int  i = offsets_[irow];
    int  k = offsets_[irow + 1];

    for ( int j = i; j < k; j++ )
    {
      accu[indices_[j]] = values_[j];
    }

    sort ( indices_[slice(i,k)] );

    for ( int j = i; j < k; j++ )
    {
      values_[j] = accu[indices_[j]];
    }
  }

  f.open("M");

  f << n << " " << n << std::endl;

  for ( idx_t irow = 0; irow <= n; irow++) f << offsets_[irow] << std::endl;
  for ( idx_t i = 0; i < nnz; i++) f << indices_[i] << " " << values_[i] << std::endl;

  f.close();
}

