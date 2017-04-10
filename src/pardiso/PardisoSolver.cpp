
/*
 *  Adapted from UmfpackSolver.cpp:
 *  
 *  Copyright (C) 2016 DRG. All rights reserved.
 *
 *  Erik Jan Lingen
 *
 *  - Fixed transpose bug, since JemJive handles Sparse Matrices in
 *    row format (CRS) and not column format (CCS).
 *  - Added debugging output
 *
 *  Copyright (C) 2017 TU Delft. All rights reserved.
 *  
 *  Frank Everdij, 3rd March 2017
 *
 */

#include <cmath>
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

#include "PardisoSolver.h"

#include <cstdio>

using jem::newInstance;
using jive::algebra::SparseMatrixExt;


//=======================================================================
//   class PardisoSolver
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


const char*  PardisoSolver::TYPE_NAME = "Pardiso";
const int    PardisoSolver::MAX_ITER  = 1;


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


PardisoSolver::PardisoSolver

  ( const String&        name,
    Ref<AbstractMatrix>  matrix,
    Ref<Constraints>     cons ) :

    Super ( name )

{
  using jem::util::connect;
  using jive::util::joinNames;
  using jive::solver::StdConstrainer;
  using jive::solver::DummyConstrainer;

  jem::System::debug ( myName_ )  << myName_ << "> Pardiso constructor\n";

  JEM_PRECHECK ( matrix != NIL );

  String    conName = joinNames ( myName_, "constrainer" );

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
  mode_      = 0;
  precision_ = PRECISION;
  updated_   = false;
  error_     = 0;
  solver_    = 0;
  mtype_     = 11;  // Real Unsymmetrical Matrix
  nrhs_      = 1;   // Number of right hand sides

  // Initialize Pardiso and check license

  pardisoinit (pt_,  &mtype_, &solver_, iparm_, dparm_, &error_);

  if ( error_ != 0 )
  {
    switch(error_)
    {
      case -10 :
        throw jem::IllegalInputException (
          getContext (),
          getContext () + " no license file found"
        );
        break;
      case -11 :
        throw jem::IllegalInputException (
          getContext (),
          getContext () + " license is expired"
        );
        break;
      case -12 : 
        throw jem::IllegalInputException (
          getContext (),
          getContext () + " username or hostname not found"
        );
        break;
      default :
        throw jem::IllegalInputException (
          getContext (),
          getContext () + " init error"
        );
        break;
    }
  }
  const char * var = getenv("OMP_NUM_THREADS");
  if(var != NULL)
  {
    sscanf( var, "%d", &iparm_[2] );
  }
  else
  {
    throw jem::IllegalInputException (
      getContext (),
      getContext () + " OMP_NUM_THREADS variable not set"
    );
  }
    
  maxfct_ = 1;
  mnum_   = 1;
    
  msglvl_ = 0;
  error_  = 0;   /* Initialize error flag */

  // Delete the work array if the matrix is modified.
  connect ( matrix_->newValuesEvent, this, & Self::freePt_ );
  connect ( matrix_->newStructEvent, this, & Self::freePt_ );
}


PardisoSolver::~PardisoSolver ()
{
  jem::System::debug ( myName_ )  << myName_ << "> Pardiso destructor\n";

  freePt_ ();
}


//-----------------------------------------------------------------------
//   improve  
//-----------------------------------------------------------------------


void PardisoSolver::improve

  ( const Vector&  lhs,
    const Vector&  rhs )

{
  using jem::dot;
  using jem::max;
  using jem::Float;
  using jem::System;
  using jem::ArithmeticException;
  using jive::solver::SolverException;

  jem::System::debug ( myName_ )  << myName_ << "> Pardiso improve\n";

  JEM_PRECHECK ( lhs.size() == rhs.size() );

  Matrix  buf;
  Vector  du;
  Vector  r;
  Vector  u;
  Vector  f;

  double  rscale;
  double  res_error;


  // Apply the constraints and update the symbolic/numeric
  // factorization if necessary.

  conman_->update ();

  if ( updated_ == false )
  {
    update_ ();
  }

  // Allocate local vectors.

  buf.resize ( matrix_->size(0), 4 );

  r  .ref    ( buf[0] );
  du .ref    ( buf[1] );
  u  .ref    ( buf[2] );
  f  .ref    ( buf[3] );

  if ( mode_ & PRECON_MODE )
  {
    f = rhs;
    u = lhs;
  }
  else
  {
    conman_->getRhs  ( f, rhs );
    conman_->initLhs ( u, lhs );
  }

  res_error  = 0.0;
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

    if ( mode_ & PRECON_MODE )
    {
      lhs = u;
    }
    else
    {
      conman_->getLhs ( lhs, u );
    }

    return;
  }

  // Compute the initial residual vector.

  rscale = 1.0 / rscale;

  matrix_->matmul ( r, u );

  r = f - r;

#ifndef NDEBUG
    System::debug ( myName_ )
      << myName_ << ": Checking matrix (2) ...\n";

    pardiso_chkmatrix  (&mtype_, &n_, values_.addr (), offsets_.addr (),
                        indices_.addr (), &error_);
    if (error_ != 0) {
      throw ArithmeticException (
        getContext (),
        " error in improve, matrix failed check"
      );
    }

    System::debug ( myName_ )
      << myName_ << ": Checking rhs ...\n";

    pardiso_chkvec (&n_, &nrhs_, r.addr (), &error_);
    if (error_ != 0) {
      throw ArithmeticException (
        getContext (),
        " error in improve, rhs vector failed check"
      );
    }
#endif

  for ( int iiter = 0; iiter < MAX_ITER; iiter++ )
  {
    // Compute the new soluton and residual vector.

    phase_ = 33;

    iparm_[7] = MAX_ITER;  // Max numbers of iterative refinement steps
    iparm_[11] = 0;        // Solving with untransposed matrix
   
    pardiso (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
             &n_, values_.addr (), offsets_.addr (), indices_.addr (), &idum_, &nrhs_,
             iparm_, &msglvl_, r.addr (), du.addr (), &error_,  dparm_);
   
    if ( error_ != 0 )
    {
      throw ArithmeticException (
        getContext (),
        String::format ("error during solution: %d", error_ )
      );
    }

    u += du;

    matrix_->matmul ( r, u );

    r     = f - r;
    res_error = rscale * std::sqrt ( dot( r, r ) );

    if ( Float::isNaN( res_error ) )
    {
      throw ArithmeticException (
        getContext (),
        "invalid norm of residual vector: NaN"
      );
    }

    if ( res_error <= precision_ || res_error > 1.0e5 )
    {
      break;
    }
  }

  // Check if the solution has been obtained.

  if ( (res_error > max( 1.0, precision_ )) ||
       (res_error > precision_ && ! (mode_ & LENIENT_MODE)) )
  {
    throw SolverException (
      getContext     (),
      String::format ( "residual norm too large: %e", res_error )
    );
  }

  if ( mode_ & PRECON_MODE )
  {
    lhs = u;
  }
  else
  {
    conman_->getLhs ( lhs, u );
  }
}


//-----------------------------------------------------------------------
//   setMode
//-----------------------------------------------------------------------


void PardisoSolver::setMode ( int mode )
{
  jem::System::debug ( myName_ )  << myName_ << "> Pardiso setMode\n";

  mode_ = mode;
}


//-----------------------------------------------------------------------
//   getMode
//-----------------------------------------------------------------------


int PardisoSolver::getMode () const
{
  jem::System::debug ( myName_ )  << myName_ << "> Pardiso getMode\n";

  return mode_;
}


//-----------------------------------------------------------------------
//   setPrecision
//-----------------------------------------------------------------------


void PardisoSolver::setPrecision ( double eps )
{
  jem::System::debug ( myName_ )  << myName_ << "> Pardiso setPrecision\n";

  JEM_PRECHECK ( eps >= 0.0 );

  precision_ = eps;
}


//-----------------------------------------------------------------------
//   getPrecision
//-----------------------------------------------------------------------


double PardisoSolver::getPrecision () const
{
  jem::System::debug ( myName_ )  << myName_ << "> Pardiso getPrecision\n";

  return precision_;
}


//-----------------------------------------------------------------------
//   getMatrix
//-----------------------------------------------------------------------


AbstractMatrix* PardisoSolver::getMatrix () const
{
  jem::System::debug ( myName_ )  << myName_ << "> Pardiso getMatrix\n";

  return conman_->getInputMatrix ();
}


//-----------------------------------------------------------------------
//   getConstraints
//-----------------------------------------------------------------------


Constraints* PardisoSolver::getConstraints () const
{
  jem::System::debug ( myName_ )  << myName_ << "> Pardiso getConstraints\n";

  return conman_->getConstraints ();
}


//-----------------------------------------------------------------------
//   makeNew
//-----------------------------------------------------------------------


Ref<Solver> PardisoSolver::makeNew

  ( const String&      name,
    const Properties&  conf,
    const Properties&  props,
    const Properties&  params,
    const Properties&  globdat )

{
  using jive::solver::SolverParams;

  Ref<AbstractMatrix>  mat;
  Ref<Constraints>     cons;

  jem::System::debug ()  << "> Pardiso makeNew\n";


  params.find ( mat,  SolverParams::MATRIX      );
  params.find ( cons, SolverParams::CONSTRAINTS );

  if ( mat != NIL && mat->hasExtension<SparseMatrixExt>() )
  {
    return newInstance<Self> ( name, mat, cons );
  }
  else
  {
    return NIL;
  }
}


//-----------------------------------------------------------------------
//   declare
//-----------------------------------------------------------------------


void PardisoSolver::declare ()
{
  using jive::solver::SolverFactory;

  jem::System::debug ()  << "> Pardiso declare\n";

  SolverFactory::declare ( TYPE_NAME,  & makeNew );
  SolverFactory::declare ( CLASS_NAME, & makeNew );
}


//-----------------------------------------------------------------------
//   update_
//-----------------------------------------------------------------------


void PardisoSolver::update_ ()
{
  using jem::maxOf;
  using jem::castTo;
  using jem::System;
  using jem::ArithmeticException;
  using jive::SparseMatrix;
  using jem::slice;

  jem::System::debug ( myName_ )  << myName_ << "> Pardiso update_\n";

  SparseMatrixExt*  sx  =

    matrix_->getExtension<SparseMatrixExt>  ();

  SparseMatrix      sm = sx->toSparseMatrix ();

  n_ = sm.size (0);
  const idx_t       nnz = sm.nonZeroCount (); // Number of non-zeroes.
  Vector            accu ( n_ );


  JEM_PRECHECK ( sm.size(0) == sm.size(1) );
  JEM_PRECHECK ( n_ <= maxOf<int>() );

  offsets_.resize ( n_ + 1 );  // Ap
  indices_.resize ( nnz );    // Ai
  values_ .resize ( nnz );    // Ax

  offsets_ = castTo<int> ( sm.getRowOffsets() );
  indices_ = castTo<int> ( sm.getColumnIndices() );
  values_  =               sm.getValues ();

  // since indexing needs to be changed, all update commands are
  // scoped in the following 'if'

  if ( ! updated_ )
  {
    // Sort the matrix entries according to their column index.

    for ( idx_t irow = 0; irow < n_; irow++ )
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

    // switch to Fortran indexing, since that is what Pardiso wants

    for ( idx_t i = 0; i < n_ + 1; i++ ) offsets_[i]++;
    for ( idx_t i = 0; i < nnz; i++ ) indices_[i]++;

#ifndef NDEBUG
    System::debug ( myName_ )
      << myName_ << ": Checking matrix (1) ...\n";

    pardiso_chkmatrix  (&mtype_, &n_, values_.addr (), offsets_.addr (),
                        indices_.addr (), &error_);
    if (error_ != 0) {
      throw ArithmeticException (
        getContext (),
        " error, matrix failed first check"
      );
    }
#endif

    System::debug ( myName_ )
      << myName_ << ": performing symbolic factorization ...\n";

    phase_ = 11; 

    pardiso (pt_, &maxfct_, &mnum_, &mtype_, &phase_, &n_, values_.addr (),
             offsets_.addr (), indices_.addr (), &idum_, &nrhs_,
             iparm_, &msglvl_, &ddum_, &ddum_, &error_,  dparm_);

    if ( error_ != 0 )
    {
      throw ArithmeticException (
        getContext (),
        String::format ("error during symbolic factorization: %d", error_ )
      );
    }

    System::debug ( myName_ )
      << myName_ << ": performing numeric factorization ...\n";

    phase_ = 22;

    pardiso (pt_, &maxfct_, &mnum_, &mtype_, &phase_, &n_, values_.addr (),
             offsets_.addr (), indices_.addr (), &idum_, &nrhs_,
             iparm_, &msglvl_, &ddum_, &ddum_, &error_, dparm_);

    if ( error_ != 0 )
    {
      throw ArithmeticException (
        getContext (),
        String::format ("error during numerical factorization: %d", error_ )
      );
    }
  }

  updated_ = true; // signal that pardiso factorization is up to date

  matrix_->resetEvents ();
}


//-----------------------------------------------------------------------
//   freePt_
//-----------------------------------------------------------------------


void PardisoSolver::freePt_ ()
{
  jem::System::debug ( myName_ )  << myName_ << "> Pardiso freePt_\n";

  phase_ = -1; // Release internal memory

  pardiso (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
           &n_, &ddum_, offsets_.addr (), indices_.addr (), &idum_, &nrhs_,
           iparm_, &msglvl_, &ddum_, &ddum_, &error_,  dparm_);

  updated_ = false;
}
