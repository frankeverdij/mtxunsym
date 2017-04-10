#ifndef PARDISOSOLVER_H
#define PARDISOSOLVER_H

#include <jive/solver/Solver.h>
#include <jive/solver/Constrainer.h>

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
//   Pardiso C prototype declaration
//-----------------------------------------------------------------------

extern "C" {
  void pardisoinit (  void * pt   ,     int * mtype,   int * solver,
                       int * iparm,  double * dparm,   int * error  );
  void pardiso     (  void * pt   ,     int * maxfct,  int * mnum,
                       int * mtype,     int * phase,   int * n,
                    double * values,    int * icol,    int * rows,
                       int * idummy,    int * nrhs,    int * iparm,
                       int * msglvl, double * rhs,  double * x,
                       int * error,  double * dparm                 );

  void pardiso_chkmatrix  (   int * mtype, int * n,    double * values,
                              int * icol,  int * rows,    int * error  );
  void pardiso_chkvec     (   int * n,     int * nrhs, double * rghs,
                              int * error                              );
  void pardiso_printstats (   int * mtype, int * n,    double * values,
                              int * icol,  int * rows,    int * nrhs,
                           double * rhs,   int * error                 );
}


//-----------------------------------------------------------------------
//   class Solver
//-----------------------------------------------------------------------


class PardisoSolver : public Solver
{
 public:

  typedef PardisoSolver     Self;
  typedef Solver            Super;

  static const char*        TYPE_NAME;
  static const int          MAX_ITER;


  explicit                  PardisoSolver

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

  virtual                  ~PardisoSolver   ();


 private:

  void                      update_         ();
  void                      freePt_    ();


 private:

  Ref<AbstractMatrix>       matrix_;
  Ref<Constrainer>          conman_;
  jem::Array<int>           offsets_;
  jem::Array<int>           indices_;
  Vector                    values_;

  int                       mode_;
  double                    precision_;
  bool                      updated_;

  void*    pt_[64];   // Internal solver memory pointer
  double   dparm_[64];
  int      iparm_[64];

  int      n_;        // Matrix dimension
  int      solver_;   // Solver selection flag
  int      mtype_;    // Matrix type and symmetry
  int      nrhs_;     // Number of right hand sides
  int      maxfct_;   // Maximum number of numerical factorizations
  int      mnum_;     // Factorization choice
  int      phase_;    // Selection number for solver task
  int      error_;    // Error flag
  int      msglvl_;   // Print statistical information

  double   ddum_;     // Double dummy
  int      idum_;     // Integer dummy

};


#endif

