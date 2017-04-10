.. Matrix Market Unsymmetric Sparse Solver documentation master file, created by
   sphinx-quickstart on Fri Mar  3 17:31:25 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Matrix Market Unsymmetric Sparse Solver's documentation!
===================================================================

.. toctree::
   :maxdepth: 2

Introduction
============

This document explains the JemJive program **mtxuss** which reads a matrix in Matrix Market format and factorizes and solves the matrix against a predetermined RHS vector.

..  math:: A x = b \rightarrow x = A^{-1} b

It was written to act as a unit-test for the Umfpack and Pardiso solver, because we recently discovered bugs in the implementation of both interface drivers. We concluded we needed a more rigorous test-program and also a comparison against the implemented JemJive solvers Skyline and SparseLU.

The design was to program a small matrix reader for a generally available and simple matrix file format and output the resulting LHS vector so that numerical difference programs like **kdiff3**, **fldiff**, **meld** and **numdiff** can compare the LHS results among each direct solver used.

Solvers
-------

**Umfpack** is an unsymmetric multifrontal sparse solver and is part of the SuiteSparse software suite written by Prof Tim Davis, previously working at University of Florida, now Texas. For references, see http://faculty.cse.tamu.edu/davis/suitesparse.html

**Pardiso** is a parallel unsymmetric sparse solver developed at the University of Lugano in Switzerland by Prof Olaf Schenk and Klaus Gartner. It is closed source and requires a license, but for academic purposes the license per username per year can be obtained free of charge.

See http://pardiso-project.org for code, licensing and references.

I have chosen these two solvers to be included since they are modern fast direct solvers suitable for very large unsymmetrical matrices with large gaps in their singular values. Moreover, Pardiso can run on multiple cores providing a boost in speed compared to a single threaded solver. Both solvers need LAPACK and a fast BLAS library implementation installed, for which i've chosen openBLAS for ease of installation.

The original drivers both suffer from a bug where the sparse matrix input was assumed to be in compressed column format, but JemJive stores the actual sparse structure in compressed row format, effectively transposing the matrix when entering the solver.

For symmetric problems this transpose bug did not affect the outcome and because of this, was unnoticed for some time. The bug is fixed within the drivers for Umfpack and Pardiso in this release.

Matrix Market
-------------

The Matrix Market format stores sparse matrices, either binary or ASCII/UTF-8 files, in triplet format: rows, columns and values. In ASCII/UTF-8 each line is one triplet.

The first line starts with ``%%MatrixMarket`` and a few keywords describing the matrix, like ``matrix coordinate real general`` . All next lines starting with ``%`` are treated as comment lines.

The first line following the comment lines is a triplet storing the number of rows, the number of columns and the number of non-zero entries. The actual triplets are then listed, each on a separate line.

For a reference, see http://math.nist.gov/MatrixMarket/formats.html


Program
=======

OS and software versions
------------------------

* Ubuntu 16.04.2 is used as operating system with gcc/g++ compiler version 5.4.0 20160609.
* JemJive version is 2.2, downloaded from http://jive.dynaflow.com/downloads/
* Umfpack and AMD versions are 5.7.1 and 2.4.1 respectively.
* Pardiso library used is version 5.0.0 compiled with gcc/g++ 4.8.1 toolchain with openMP support.


Directory structure
-------------------

The repository consists of three directories **mtx**, **sphinx** and **src**.

.. made with 'tree --charset ascii'

::

    .
    |-- mtx
    |   |-- aft02.mtx
    |   |-- autounit.py
    |   |-- Baumann.mtx
    |   |-- bfwa62.mtx
    |   |-- mmread.m
    |   |-- mtxuss.m
    |   |-- pardiso.lic
    |   |-- sherman5.mtx
    |   |-- sme3Da.mtx
    |   |-- symsmall.mtx
    |   |-- ted_A.mtx
    |   |-- unsym.mtx
    |   |-- xPardiso.aft02.unittest
    |   |-- xSkylineLU.aft02.unittest
    |   |-- xSparseLU.aft02.unittest
    |   `-- xUmfpack.aft02.unittest
    |-- sphinx
    |   |-- _build
    |   |   |-- doctrees
    |   |   |   |-- environment.pickle
    |   |   |   `-- index.doctree
    |   |   `-- html
    |   |       |-- genindex.html
    |   |       |-- index.html
    |   |       |-- objects.inv
    |   |       |-- search.html
    |   |       |-- searchindex.js
    |   |       |-- _sources
    |   |       |   `-- index.txt
    |   |       `-- _static
    |   |-- conf.py
    |   |-- index.rst
    |   |-- Makefile
    |   `-- _static
    `-- src
        |-- liblapack.a -> /usr/lib/liblapack.a
        |-- libpardiso500-GNU481-X86-64.so
        |-- Makefile
        |-- Makefile.atlas
        |-- Makefile.Pardiso
        |-- Makefile.Umfpack
        |-- mtxUnsymSparseSol.cpp
        |-- pardiso
        |   |-- PardisoSolver.cpp
        |   `-- PardisoSolver.h
        `-- umfpack
            |-- UmfpackSolver.cpp
            `-- UmfpackSolver.h

mtx
---

**mtx** holds a few examples of matrices in the matrix market format from the Sparse Matrix Collection. Only ``symsmall.mtx`` is symmetric, the rest are unsymmetric real matrices of full rank. ``unsym.mtx`` is a small matrix which is also represented internally in the ``mtxUnsymSparseSol.cpp`` code for unit-testing purposes.

The files ``mmread.m`` and ``unsym.m`` are respectively a MATLAB file for reading a MTX file and a MATLAB version of the core of program **mtxuss** by reading ``unsym.mtx`` and producing and displaying a solution LHS vector.

MATLAB internally uses Umfpack for solving square unsymmetric matrices, so the result should be the same as the result from Umfpack.

For checking result vectors of the matrix ``aft02.mtx`` , its unittest files for each solver type have been added so that ``autounit.py -u`` can compare the outcome of the program when running with ``aft02.mtx`` as argument.

``pardiso.lic`` is an empty placeholder file. Substitute it with a real license obtained from the Pardiso-project website.

sphinx
------

This documentation source is in ``sphinx/index.rst`` and the html output in ``sphinx/_build/html/index.html`` can be generated with::

    make html

It uses mathjax to display the formula in the introduction but it is not essential for the rest of the document.


src
---

In **src** the source of the ``mtxuss`` program consists of a single C++ file::

    mtxUnsymSparseSol.cpp

The program has a few subroutines::

 fillRhs ( Vector & v )
     fills a Vector with decreasing integer sequence
 writeVector( const Vector & v , std::string solver , std::string mtxfile)
     writes a Vector object to a file
 createOffsets ( idx_t n, idx_t nnz, IdxVector &irow, IdxVector &jcol, Vector &val ) 
     creates row index offset array from row and column index arrays
 unsymSparse( void )
     outputs a small sparse matrix, 8 x 8 full rank, identical to unsym.mtx
 writeResidualNormRef( Vector & x )
     compares solver result with the reference LHS vector result
 readMTXtoSparse( std::string s )
     reads Matrix Market file into a SparseMatrix structure
 void doSolve( String sname, Ref<SparseMatrixObject>SP, idx_t n, int argc, std::string filename )
     template function for doing the actual solve

The main subroutine is a standard C-style argument count and pointer subroutine to parse input::

 main (int argc, char ** argv)


The solvers are each in their own subdirectory::

    umfpack/UmfpackSolver.h
    umfpack/UmfpackSolver.cpp

and::

    pardiso/PardisoSolver.h
    pardiso/PardisoSolver.cpp

This makes compiling the solvers into the program optional: This choice is done at compile time with the preprocessor macros ``-DWITH_UMFPACK`` and ``-DWITH_PARDISO``

There are four Makefiles::

    Makefile
    Makefile.Umfpack
    Makefile.Pardiso
    Makefile.atlas

The vanilla ``Makefile`` produces ``mtxuss`` and ``mtxuss-opt`` with both external solvers included. ``Makefile.umfpack`` and ``Makefile.pardiso`` are for building the program with only Umfpack and Pardiso respectively.

``Makefile.atlas`` is a makefile primarily for linking with ATLAS on the cluster. See the `ATLAS support`_ subsection in `Notes`_ section for details.


Compile and Run
===============


Compiling with Umfpack Driver
-----------------------------

For Ubuntu systems the easiest way to install Umfpack is via apt::

    sudo apt install liblapack-dev libopenblas-dev libsuitesparse-dev

This will install openBLAS as the fast BLAS library, a (fortran) reference LAPACK library and suitesparse development library which includes Umfpack and its dependencies.

Next, compile the program within the JemJive framework::

    cd src
    make -f makefile.Umfpack opt

This produces an executable ``mtxuss-opt`` with Umfpack solver support.


Compiling with Pardiso Driver
-----------------------------

Similar to Umfpack, install the BLAS and LAPACK libraries via apt:

``sudo apt install liblapack-dev libopenblas-dev``

Then compile the program within the JemJive framework::

    cd src
    make -f makefile.Pardiso opt

This produces an executable ``mtxuss-opt`` with Pardiso solver support.


Running the program
-------------------

After building the program, run it in the **mtx** directory, but if the program has Pardiso support, check that you have a valid license obtained from the Pardiso website and that you have set the number of threads in openMP via an environment variable::

    $ cd ../mtx
    $ ls -l pardiso.lic
    -rw-r--r-- 1 feverdij users 237 mrt  6 16:51 pardiso.lic
    $ export OMP_NUM_THREADS=2

Failure to set the ``OMP_NUM_THREADS`` environment variable will result in an error message from jem::

    terminate called after throwing an instance of 'jem::IllegalInputException'
    Aborted (core dumped)


If all is set, run the program without arguments::

    $ ../src/mtxuss-opt

This produces output with an internal matrix identical to ``unsym.mtx``. For each solver used, the norm of the resulting LHS vector is computed against the reference vector::

    $ ../src/mtxuss-opt
    no MTX file given. Using internal unittest sparse matrix.

    size is 8 x 8 with 20 nonzeros.
    .
    .
    Pardiso total run time is : 0.0009840  seconds
    Norm of (x - reference) = 1.110223e-16 , passed, norm is within 1.000000e-14 tol

    Umfpack total run time is : 0.0001440  seconds
    Norm of (x - reference) = 1.110223e-16 , passed, norm is within 1.000000e-14 tol

    SparseLU total run time is : 7.200e-05  seconds
    Norm of (x - reference) = 4.475452e-16 , passed, norm is within 1.000000e-14 tol

    solver `SkylineLU' : residual norm too large: 3.725144e+00
    SkylineLU : switching to the `SparseLU' solver
    SkylineLU total run time is : 0.003654  seconds
    Norm of (x - reference) = 4.475452e-16 , passed, norm is within 1.000000e-14 tol

When giving a filename to the program as extra option, the file wil be treated as an ASCII/UTF-8 Matrix Market file and the matrix will be factorized and solved::

    $ ../src/mtxuss-opt sme3Da.mtx 
    MTX file size 12504 x 12504

    size is 12504 x 12504 with 874887 nonzeros.

    Pardiso total run time is : 0.4505  seconds

    Umfpack total run time is : 0.3590  seconds

    SparseLU total run time is : 3.840  seconds

    SkylineLU total run time is : 0.9384  seconds


Comparing output
================

Comparing with unittest files
-----------------------------

Additionally, the computed LHS for each solver type will be written to file. The file format will be ``x + <solver name> + . + mtx-base-filename`` , so in the case of MTX file ``aft02.mtx``::

    xPardiso.aft02
    xUmfpack.aft02
    xSparseLU.aft02
    xSkylineLU.aft02

This enables the user to compare the output vectors between solvers and check against the unittest versions in the case of ``aft02.mtx`` with the Python script ``autounit.py``::

    $ ../src/mtxuss-opt aft02.mtx

    $ ./autounit.py -u
    unittest: In .
 +++ Files xPardiso.aft02 and xPardiso.aft02.unittest match. 
 +++ Files xSparseLU.aft02 and xSparseLU.aft02.unittest match. 
 +++ Files xUmfpack.aft02 and xUmfpack.aft02.unittest match. 
 +++ Files xSkylineLU.aft02 and xSkylineLU.aft02.unittest match.


Comparing among solvers
-----------------------

A compare of the output vectors among solver is useful to see if the solvers produce similar output. Re-using the output files of ``aft02.mtx`` we can f.i. apply ``numdiff`` manually on Umfpack and SparseLU output::

    $ numdiff -r 1.0e-7 xUmfpack.aft02 xSparseLU.aft02

    +++  Files "xUmfpack.aft02" and "xSparseLU.aft02" are equal

Only at a much lower tolerance of :math:`10^{-9}` we will see differences between these two vectors::

    $ numdiff -r 1.0e-9 xUmfpack.aft02 xSparseLU.aft02
    ----------------
    ##994     #:1   <== -18.17861607954626
    ##994     #:1   ==> -18.17861611238935
    @ Absolute error = 3.2843090000e-8, Relative error = 1.8066881360e-9


Comparing with MATLAB
---------------------

In the **mtx** directory, a comparison of the ``mtxuss-opt unsym.mtx`` and MATLAB can be made. First, run the program with unsym.mtx as argument::

    $ ../src/mtxuss-opt unsym.mtx

Then display the output file ``xUmfpack.unsym``::

    $ cat xUmfpack.unsym 
    8 1
    2.013636363636364
    -1
    1.386363636363636
    0.3636363636363637
    -4.045454545454545
    -4.695454545454545
    0.2727272727272727
    0.9227272727272727

It should be the same (within double precision tolerance) with the result from ``mtxuss.m`` in MATLAB::

    $ matlab -nodisplay -nosplash < mtxuss.m 

                            < M A T L A B (R) >
                  Copyright 1984-2015 The MathWorks, Inc.
                   R2015a (8.5.0.197613) 64-bit (glnxa64)
                             February 12, 2015

 
    To get started, type one of these: helpwin, helpdesk, or demo.
    For product information, visit www.mathworks.com.
 

	Academic License

    >> >> >> >> 
    x =

       2.013636363636364
      -1.000000000000000
       1.386363636363636
       0.363636363636364
      -4.045454545454546
      -4.695454545454546
       0.272727272727273
       0.922727272727273



Notes
=====

More MTX files
--------------

Users can test other MTX files. These can be found on the website for The SuiteSparse Matrix Collection at https://www.cise.ufl.edu/research/sparse/matrices/

As a test with a really large unsymmetrical matrix, try to run the program with ``stomach.mtx`` as input. The file can be downloaded `here <http://www.cise.ufl.edu/research/sparse/MM/Norris/stomach.tar.gz>`_

For Pardiso and Umfpack this matrix can be solved withint seconds, but for the other solvers the output can take a long time::

    $ ../src/mtxuss-opt stomach.mtx
    MTX file size 213360 x 213360

    size is 213360 x 213360 with 3021648 nonzeros.
    .
    .
    Pardiso total run time is : 6.741  seconds

    Umfpack total run time is : 7.638  seconds
    .
    .


ATLAS support
-------------

ATLAS comes with its own LAPACK library, but it is only a subset of the full LAPACK API and thus not usable with Pardiso. In order to link with a full LAPACK library in the JemJive framework you have to provide a static archive ``liblapack.a`` as the first library when linking, otherwise it will try to link the dynamical library ``/usr/lib/liblapack.so`` and will include ``/usr/lib/libblas.so`` as dependency.
This is bad because that ``libblas.so`` usually is the slow reference library, so we need to find another way of providing LAPACK to the program.

As a workaround for this, add to the library options line an extra dot at the beginning::

    MY_LIBDIRS	= . /opt/lib /usr/local/atlas/lib

and a symlink to /usr/lib/liblapack.a ::

   ln -s /usr/lib/liblapack.a liblapack.a 

This statically links a full ``liblapack.a`` archive before anything else.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

