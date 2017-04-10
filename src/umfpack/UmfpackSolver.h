#ifndef UMFPACKSOLVER_H
#define UMFPACKSOLVER_H

#include <jive/util/error.h>
#include <jive/solver/Solver.h>
#include <jive/solver/Constrainer.h>
#include <umfpack.h>
#include <SuiteSparse_config.h>

using jem::NIL;
using jem::idx_t;
using jem::Ref;
using jem::String;
using jem::util::Properties;
using jive::util::Constraints;
using jive::algebra::AbstractMatrix;
using jive::solver::Solver;
using jive::solver::Constrainer;
using jive::Vector;
using jive::Matrix;

//-----------------------------------------------------------------------
//   class Solver
//-----------------------------------------------------------------------


class UmfpackSolver : public Solver
{
 public:

  typedef UmfpackSolver     Self;
  typedef Solver      Super;

  static const char*        TYPE_NAME;
  static const int          MAX_ITER;


  explicit                  UmfpackSolver

    ( const String&           name,
      Ref<AbstractMatrix>     matrix,
      Ref<Constraints>        cons = NIL );

  virtual void              improve

   ( const Vector&            lhs,
     const Vector&            rhs );

  virtual void              setMode

    ( int                     mode );

  virtual int               getMode         () const;

  virtual void              setPrecision

    ( double                  eps );

  virtual double            getPrecision    () const;
  virtual AbstractMatrix*   getMatrix       () const;
  virtual Constraints*      getConstraints  () const;

  static Ref<Solver>        makeNew

    ( const String&           name,
      const Properties&       conf,
      const Properties&       props,
      const Properties&       params,
      const Properties&       globdat );

  static void               declare         ();


 protected:

  virtual                  ~UmfpackSolver   ();


 private:

  void                      update_         ();
  void                      dump_           ();
  void                      freeSymbolic_   ();
  void                      freeNumeric_    ();
  void                      valuesChanged_  ();
  void                      structChanged_  ();


 private:

  static const int          NEW_VALUES_;
  static const int          NEW_STRUCT_;

  Ref<AbstractMatrix>       matrix_;
  Ref<Constrainer>          conman_;
  void*                     symbolic_;
  void*                     numeric_;
  double                    control_[UMFPACK_CONTROL];
  double                    info_[UMFPACK_INFO];
  jem::Array<long>          offsets_;
  jem::Array<long>          indices_;
  Vector                    values_;

  int                       mode_;
  double                    precision_;
  bool                      updated_;

  idx_t                     iiter_;
  double                    error_;
  int                       events_;

};


#endif

