/*
 * Copyright (C) 2017 TU Delft. All rights reserved.
 *  
 * Frank Everdij, March 10nd, 2017
 *
 * mtxUnsymSparseSol: Matrix Market Unsymmetric Sparse Solver (mtxuss)
 *
 *  Adapted from poisson.cpp and
 *  jem-2.2/packages/numeric/examples/algebra/sparse.c
 *
 * Reads .mtx files in Matrix Market format and solves the system
 *  for a decreasing integer sequence rhs vector.
 *
 */

#include <jem/base/Array.h>
#include <jem/base/System.h>
#include <jem/util/Timer.h>
#include <jive/algebra/utilities.h>
#include <jive/algebra/SparseMatrixObject.h>
#include <jive/algebra/AbstractMatrix.h>
#include <jive/solver/SkylineLU.h>
#include <jive/solver/SparseLU.h>
#include <jive/util/Constraints.h>

#if defined(WITH_UMFPACK)
#include "UmfpackSolver.h"
#endif

#if defined(WITH_PARDISO)
#include "PardisoSolver.h"
#endif

#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace jem;
using namespace jem::io;
using namespace jem::util;

using jive::SparseMatrix;
using jive::algebra::SparseMatrixObject;
using jive::IdxVector;

typedef Array<double>    Vector;


// fills a Vector with decreasing integer sequence

void fillRhs ( Vector & v )
{
    idx_t n = v.size(0);

    for ( idx_t i = 0; i < n; i++ )
    {
        v[i] = (double) (n - i);
    }
}


// writes a Vector object to a file

void writeVector( const Vector & v , std::string solver , std::string mtxfile)
{
    std::string file(1, 'x');
    std::ofstream f;
    idx_t n = v.size(0);

    file += solver;
    file += '.';
    file += mtxfile;
    file.erase( file.find_last_of( '.' ) , std::string::npos );

    f.open(file.c_str());

    f << std::setprecision(16) << n << " 1" << std::endl;

    for (idx_t i = 0; i< n; i++ )
    {
        f << (double) v[i] << std::endl;
    }

    f.close();
}


// creates row index offset array from row and column index arrays

IdxVector createOffsets ( idx_t n, idx_t nnz,
                          IdxVector &irow, IdxVector &jcol, Vector &val )
{
    IdxVector offset( n + 1 );
    
    // sort the indices row first, then column
    // Shellsort, Knuth sequence

    int cols[] = {1391376, 463792, 198768, 86961, 33936, 13776, 4592,
                    1968, 861, 336, 112, 48, 21, 7, 3, 1};
    idx_t j,iv,jv,h;
    double dv;

    for ( idx_t k = 0; k < 16; k++ )
    {
        h = cols[k];
        for ( idx_t i = h; i < nnz; i++ )
        {
            iv = irow[i];
            jv = jcol[i];
            dv = val[i];
            j = i;
            while (j >= h && ( irow[j-h] > iv ||
                               (irow[j-h] == iv && jcol[j-h] > jv)))
            {
                irow[j] = irow[j-h];
                jcol[j] = jcol[j-h];
                val[j] = val[j-h];
                j = j - h;
            }
            irow[j] = iv;
            jcol[j] = jv;
            val[j] = dv;
        }
    }

    // clear offset vector
    offset = 0;

    // first count how many entries there are for each index
    idx_t check = 0;
    for ( idx_t i = 0; i < nnz; i++ )
    {
        if ( irow[i] >= check )
        {
            offset[irow[i]-1]++;
            check = irow[i];
        }
        else
        {
            // Something went wrong, row index vector wasn't sorted
            System::err() << "createOffsets: row index vector not sorted!" <<
                             endl;
            exit(-1);
        }

        // in the meantime change both column and row index to 0-based
        jcol[i]--;
        irow[i]--;

    }

    // finally, construct the CRS array `offset` from the index count array
    check = 0;
    for ( idx_t i = 0; i <= n; i++ )
    {
        j = offset[i];
        offset[i] = check;
        check += j;
    }

    return offset;
}


// outputs a small sparse matrix, 8 x 8 full rank, identical to unsym.mtx

SparseMatrix unsymSparse( void )
{
    int   io[] = {0, 4, 7, 9, 11, 12, 15, 17, 20};
    int    j[] = {0, 2, 5, 6, 1, 2, 4, 2, 7, 3, 6, 1, 2, 5, 7, 1, 6, 2, 6, 7};
    double v[] = {7.0, 1.0, 2.0, 7.0, -4.0, 8.0, 2.0, 1.0, 5.0, 7.0, 9.0,
                  -4.0, 7.0, 3.0, 8.0, 1.0, 11.0, -3.0, 2.0, 5.0};
    idx_t n = 8;
    idx_t nnz = 20;
    IdxVector Offsets( n + 1 );
    IdxVector jIndices( nnz );
    Vector Values( nnz );

    for (idx_t i=0; i < nnz; i++)
    {
        jIndices[i] = j[i];
        Values[i] = v[i];
    }

    for (idx_t i=0; i <= n; i++)
    {
        Offsets[i] = io[i];
    }

    return SparseMatrix(jem::shape( n, n ),
                        Offsets,
                        jIndices,
                        Values             );
}


// compares solver result with the reference Lhs vector result from
//  x = A\b where A equals the unsymSparse() matrix and
//  b equals [8, 7, 6, 5, 4, 3, 2, 1] in column format.

void writeResidualNormRef( Vector & x )
{
    double reference[] = {2.0136363636363636, -1.0,
                          1.3863636363636363,  0.3636363636363636,
                         -4.0454545454545454, -4.6954545454545454,
                          0.2727272727272727,  0.9227272727272727};
    double norm = 0.0;
    double tol = 1.0e-14;

    for (idx_t i = 0; i < min(8, x.size(0)); i++ )
    {
        norm += ( x[i] - reference[i] ) * ( x[i] - reference[i] ) ;
    }

    norm = sqrt(norm);

    std::cout << std::scientific << "Norm of (x - reference) = " <<
              norm;

    if ( std::abs( norm ) < tol )
    {
        std::cout << " , passed, norm is within " << tol << " tol" <<
                  std::endl;
    }
    else
    {
        std::cout << " , failed, norm is bigger than " << tol << " tol" <<
                  std::endl;
    }
}


// reads Matrix Market file into a SparseMatrix structure

SparseMatrix readMTXtoSparse( std::string s )
{
    std::string     dummy;
    std::ifstream   mtxfile;

    idx_t           i, ir, j, n, nnz;
    double          a;

    // open MTX file for reading, c-style
    mtxfile.open(s.c_str());

    if ( mtxfile.is_open() )
    {
        // skip lines with % as first character
        do
        {
            getline( mtxfile, dummy );
        }
        while (dummy[0] == '%');

        // next line contains m,n, and nonzeros
        std::sscanf( dummy.c_str(), "%d %d %d", &n, &j, &nnz );
        System::out() << "MTX file size " << n << " x " << j << endl;
        if ( n != j )
        {
            System::err() << "MTX file: Matrix is not square!" << endl;
            exit(-1);
        }

        // allocate sparse vectors
        IdxVector iIndices ( nnz );
        IdxVector jIndices ( nnz );
        Vector Values (nnz);

        // read all triplets: row index, column index and values
        ir = 0;
        while (ir < nnz)
        {
            getline( mtxfile, dummy );
            std::sscanf( dummy.c_str(), "%d %d %lg", &i, &j, &a );
            iIndices[ir] = i;
            jIndices[ir] = j;
            Values[ir] = a;
            System::debug() << "triplets " << i << " " << j << " " << a <<
                             endl;
            ir++;
        }

        // make an offset array for the rows
        IdxVector Offsets = createOffsets(n, nnz, iIndices, jIndices, Values);

        for (ir=0;ir<nnz;ir++)
        {
            System::debug() << "sorted triplets " << iIndices[ir] << " " <<
                             jIndices[ir] << " " << Values[ir] << endl;
        }
        for (ir=0;ir<=n;ir++)
        {
            System::debug() << Offsets[ir] << " ";
        }

        // finally create the Sparse Array
        return SparseMatrix( jem::shape( n, n ),
                             Offsets,
                             jIndices,
                             Values             );
    } else {

        System::err() <<
         "file open failed, using internal unittest sparse matrix." << endl;
        return unsymSparse();

    }
}


// template function : solves for a particular solver class the system of
// linear equations and does start-stop timing

template<class T> void doSolve( String sname, Ref<SparseMatrixObject>SP,
                                idx_t n, int argc, std::string filename )
{
    Ref<jive::util::Constraints> cons = jem::NIL;
    Vector b,x;
    Timer tmr;
    double seconds;

    Ref<T> solver = newInstance<T> ( sname, SP, cons );

    System::out() << endl << solver->TYPE_NAME;

    x.resize(n);
    b.resize(n);
    fillRhs(b);

    tmr.reset();
    tmr.start();
    solver->solve ( x, b );
    tmr.stop();
    seconds = tmr.toDouble();
    System::out() << " total run time is : " << seconds  <<
                     "  seconds" << endl;

    if (argc == 1)
    {
        writeResidualNormRef( x );
    }
    else
    {
        writeVector ( x, solver->TYPE_NAME, filename);
    }
}

// Main program

int main (int argc, char ** argv)
{
    std::string filename;
    SparseMatrix A;
    int traits = 0;
    idx_t n;
    Ref<SparseMatrixObject> SP;

    if (argc == 1)
    {
        System::out() <<
         "no MTX file given. Using internal unittest sparse matrix." << endl;
    
        A = unsymSparse();
    }
    else
    {
        filename = argv[1];

        A = readMTXtoSparse( filename );
    }

    n = A.size(0);

    System::out() << endl << "size is " << n << " x "<< A.size(1);
    System::out() << " with " << A.nonZeroCount() << " nonzeros." << endl;

    SP = newInstance<SparseMatrixObject> ("sparse", A, traits);

#if defined(WITH_UMFPACK) 
    doSolve<UmfpackSolver>("Umfpack", SP, n, argc, filename );
#endif

#if defined(WITH_PARDISO) 
    doSolve<PardisoSolver>("Pardiso", SP, n, argc, filename );
#endif

    doSolve<jive::solver::SparseLU>("SparseLU", SP, n, argc, filename);
    doSolve<jive::solver::SkylineLU>("SkylineLU", SP, n, argc, filename);

    System::out() << endl;

    return 0;

}
