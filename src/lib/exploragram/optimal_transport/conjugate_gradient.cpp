
/*
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000-2009 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy - levy@loria.fr
 *
 *     Project ALICE
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 *
 * As an exception to the GPL, Graphite can be linked with 
 *     the following (non-GPL) libraries:
 *     Qt, SuperLU, WildMagic and CGAL
 */

#include <exploragram/optimal_transport/conjugate_gradient.h>
#include <geogram/NL/nl.h>
#include <geogram/basic/stopwatch.h>

namespace {
    using namespace GEO;
    Numeric::uint64 flops_ = 0;

    /**
     * \brief Computes the dot product between two vectors.
     * \param[in] n dimension of the vectors
     * \param[in] x,y const pointers to the coefficients of the vectors
     * \return the dot product between the two vectors
     */
    double dot(index_t n, const double* x, const double* y) {
        double result = 0.0;
        for(index_t i=0; i<n; ++i) {
            result += x[i]*y[i];
        }
        flops_ += Numeric::uint64(2*n);
        return result;
    }


    /**
     * \brief Computes the squared norm of a vector.
     * \param[in] n dimension of the vectors
     * \param[in] x a const pointer to the coefficients of the vector
     * \return the squared norm of the vector
     */
    double nrm2(index_t n, const double* x) {
        double result = 0.0;
        for(index_t i=0; i<n; ++i) {
            result += x[i]*x[i];
        }
        flops_ += Numeric::uint64(2*n);        
        return result;
    }

    /**
     * \brief Computes the linear combination between two vectors.
     * \details In formula: \f$ y \leftarrow a x + y \f$
     * \param[in] n the dimension of the vectors
     * \param[in] a a scalar
     * \param[in] x a const pointer to a vector
     * \param[in,out] y a pointer to a vector
     */
    void axpy(
        index_t n, double a, const double* x, double* y
    ) {
        for(index_t i=0; i<n; ++i) {
            y[i] = a*x[i]+y[i];
        }
        flops_ += Numeric::uint64(2*n);                
    }

    /**
     * \brief Copies a vector.
     * \param[in] n the dimension of the vectors
     * \param[in] x the source vector
     * \param[out] y the destination vector
     */
    void copy(index_t n, const double* x, double* y) {
        Memory::copy(y,x,sizeof(double)*n);
    }
    
    /**
     * \brief Scales a vector.
     * \param[in] n the dimension of the vector
     * \param[in] s the scaling factor
     * \param[in,out] x the vector to be scaled
     */
    void scal(index_t n, double s, double* x) {
        for(index_t i=0; i<n; ++i) {
            x[i] *= s;
        }
    }

    /**
     * \brief Computes a sparse matrix-vector product.
     * \param[in] n the dimension of the square matrix
     * \param[in] x the input vector
     * \param[out] y the result
     * \param[in] rowptr, colind, val the sparse matrix in compressed row
     *  storage format
     * \param[in] symmetric_storage true if the matrix is stored with symmetric
     *    storage (i.e. only lower triangular part is stored), false otherwise
     *    (i.e. then all the coefficients of the matrix are stored).
     */
    void spMV(
        index_t n, const double* x, double* y,
        index_t* rowptr, index_t* colind, double* val,
        bool symmetric
    ) {
        // TODO: balanced parallel version.

        index_t nnz = rowptr[n];
        
        if(symmetric) {
            Memory::clear(y,sizeof(double)*n);
            for(index_t i=0; i<n; ++i) {
                for(index_t jj=rowptr[i]; jj<rowptr[i+1]; ++jj) {
                    double a = val[jj];
                    index_t j = colind[jj];
                    y[i] += a * x[j];
                    if(j != i) {
                        y[j] += a * x[i];
                    }
                }
            }
            // There are two additions + two products for each non-zero entry
            flops_ += Numeric::uint64(4*nnz);
            // ... except for diagonal entries.
            flops_ -= Numeric::uint64(2*n);
        } else {
#ifdef GEO_OPENMP	    
#pragma omp parallel for
#endif	    
            // Note: the iteration variable for an OpenMP loop needs 
            // to be a *signed* integers (else it generates a FPE !!!)
            for(int i=0; i<int(n); ++i) {
                y[i] = 0.0;
                for(index_t jj=rowptr[i]; jj<rowptr[i+1]; ++jj) {
                    y[i] += val[jj] * x[colind[jj]];
                }
            }
            // There is an addition + a product for each non-zero entry.
            flops_ += Numeric::uint64(2*nnz);
        }
    }
    
    /**
     * \brief Computes a matrix-vector product with the preconditioner.
     * \param[in] n the dimension of the preconditioner
     * \param[out] y the result
     * \param[in] x the input vector
     * \param[in] P a pointer to the preconditioner, i.e. an array of 
     *  n doubles with the inverse of the diagonal coefficients of the matrix
     */
    void precondMV(index_t n, const double* x, double* y, const double* P) {
        for(index_t i=0; i<n; ++i) {
            y[i] = P[i]*x[i];
        }
        flops_ += Numeric::uint64(n);
    }


#ifdef USE_OPENNL
    
    /**
     * \brief Solves a linear system with OpenNL
     * \details Just used for tests and debugging (quite unefficient, it copies
     *  the matrix to OpenNL).
     * \param[in] n dimension of the system
     * \param[in] rowptr, colind, val the matrix of the system in the 
     *   compressed row storage format
     * \param[in] rhs the right-hand side of the system
     * \param[out] x the solution of the linear system
     * \param[in] max_iter maximum number of iterations
     * \param[in] threshold the maximum value of 
     *    \f$ \| Ax - b \| / \| b \| \f$ before the iterative solver is stopped.
     * \param[in] symmetric_storage true if the matrix is stored with symmetric
     *    storage (i.e. only lower triangular part is stored), false otherwise
     *    (i.e. then all the coefficients of the matrix are stored).
     */
    void solve_with_OpenNL(
        index_t  n,
        index_t* rowptr,
        index_t* colind,
        double*  val,
        double*  rhs,
        double*  x,
        index_t  max_iter,
        double   threshold,
        bool     symmetric_storage
    ) {
        nlNewContext();
        nlSolverParameteri(NL_NB_VARIABLES, NLint(n));
        nlSolverParameteri(NL_SOLVER, NL_CG);
        nlSolverParameteri(NL_PRECONDITIONER, NL_PRECOND_JACOBI);
        nlSolverParameteri(NL_SYMMETRIC, symmetric_storage ? NL_TRUE: NL_FALSE);
        nlSolverParameterd(NL_THRESHOLD, threshold);
        nlSolverParameteri(NL_MAX_ITERATIONS, NLint(max_iter));                
        nlBegin(NL_SYSTEM);
        nlBegin(NL_MATRIX);

        for(index_t i=0; i<n; ++i) {
            nlAddIRightHandSide(i,rhs[i]);
            for(index_t jj=rowptr[i]; jj<rowptr[i+1]; ++jj) {
                nlAddIJCoefficient(i, colind[jj], val[jj]);
            }
        }
        
        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);
        nlSolve();
        
        {
            int used_iters;
            double elapsed_time;
            double gflops;
            double error;
            nlGetIntegerv(NL_USED_ITERATIONS, &used_iters);
            nlGetDoublev(NL_ELAPSED_TIME, &elapsed_time);
            nlGetDoublev(NL_GFLOPS, &gflops);
            nlGetDoublev(NL_ERROR, &error);                
            Logger::out("OpenNL")
                << "   "
                << used_iters << " iters in "
                << elapsed_time << " seconds "
                << gflops << " GFlop/s"
                << "  ||Ax-b||/||b||="
                << error
                << std::endl;
        }
        for(index_t i=0; i<n; ++i) {
            x[i] = nlGetVariable(NLuint(i));
        }
        nlDeleteContext(nlGetCurrent());
    }
#endif
    
}



namespace GEO {

    void solve_conjugate_gradient(
        index_t  n,
        index_t* rowptr,
        index_t* colind,
        double*  val,
        double*  rhs,
        double*  x,
        index_t  max_iter,
        double   threshold,
        bool     symmetric_storage
    ) {


	// This one is just for debugging, it is less efficient
	// (copies everything to OpenNL, solves, copies back)
#ifdef USE_OPENNL	
	solve_with_OpenNL(
	    n, rowptr, colind, val, rhs, x, max_iter, threshold,
	    symmetric_storage
	);
	return;            
#endif

        flops_ = 0;
        
        // step 1: compute the inverse of the diagonal (Jacobi
        // preconditioner).
        
        double* P = new double[n];
        for(index_t i=0; i<n; ++i) {
            for(index_t jj=rowptr[i]; jj<rowptr[i+1]; ++jj) {
                if(colind[jj] == i) {
                    P[i] = 1.0 / val[jj];
                }
            }
        }

        double start_time = SystemStopwatch::now();

        // step 2: conjugate gradient loop *******************

        double eps = threshold;
        const double* b = rhs;
        double* r = new double[n];
        double* d = new double[n];
        double* h = new double[n];
        double* Ad = h;
        index_t its = 0;
        double rh, alpha, beta;
        double b_square = nrm2(n,b);
        double err=eps*eps*b_square;
        index_t i;
        double* Ax = new double[n];
        double accu =0.0;
        double cur_err;

        spMV(n,x,r,rowptr,colind, val, symmetric_storage);
        axpy(n,-1.0,b,r);
        precondMV(n,r,d,P);
        copy(n,d,h);
        rh = dot(n,r,h);
        cur_err = nrm2(n,r);

        while ( cur_err >err && its < max_iter) {
            spMV(n,d,Ad,rowptr,colind, val, symmetric_storage);
            alpha=rh/dot(n,d,Ad);
            axpy(n,-alpha,d,x);
            axpy(n,-alpha,Ad,r);
            precondMV(n,r,h,P);
            beta=1./rh; rh=dot(n,r,h); beta*=rh;
            scal(n,beta,d);
            axpy(n,1.,h,d);
            ++its;
            cur_err = nrm2(n,r);
        }
        spMV(n,x,Ax,rowptr,colind, val, symmetric_storage);
        for(i = 0; i < n; ++i) {
            accu+=(Ax[i]-b[i])*(Ax[i]-b[i]);
        }
        if(b_square == 0.0) {
            double error = sqrt(accu);
            Logger::out("OTM linsolve")
                << "||Ax-b|| = " << error << std::endl;
        } else {
            double error = sqrt(accu/b_square);
            Logger::out("OTM linsolve")
                << "||Ax-b||/||b|| = " << error << std::endl;
        }

        delete[] Ax;
        delete[] r;
        delete[] d;
        delete[] h;
        delete[] P;

        double elapsed_time = SystemStopwatch::now()-start_time;
        double GFlops = double(flops_) / (elapsed_time * 1e9);

        Logger::out("OTM linsolve")
            << its << " iterations in " << elapsed_time << " seconds"
            << std::endl;

        Logger::out("OTM linsolve")
            << GFlops << " GFlops"
            << std::endl;
        
    }
}

