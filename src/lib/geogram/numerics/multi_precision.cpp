/*
 *  Copyright (c) 2004-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include <geogram/numerics/multi_precision.h>

namespace {

    using namespace GEO;


    /************************************************************************/

    class ExpansionAutoInitialize {
    public:
        ExpansionAutoInitialize() {
            expansion::initialize();
        }
    };

    ExpansionAutoInitialize init_;

    /************************************************************************/

    /**
     * \brief Computes the sum of two doubles into a length 2 expansion.
     * \details By Jonathan Shewchuk.
     * \param[in] a first argument
     * \param[in] b second argument
     * \param[out] x high-magnitude component of the result
     * \param[out] y low-magnitude component of the result
     * \pre |\p a| > |\p b|
     */
    inline void fast_two_sum(double a, double b, double& x, double& y) {
        x = a + b;
        double bvirt = x - a;
        y = b - bvirt;
    }

    /**
     * \brief Computes the difference of two doubles into a length 2 expansion.
     * \details By Jonathan Shewchuk.
     * \param[in] a first argument
     * \param[in] b second argument
     * \param[out] x high-magnitude component of the result
     * \param[out] y low-magnitude component of the result
     * \pre | \p a| > | \p b |
     */
    inline void fast_two_diff(double a, double b, double& x, double& y) {
        x = a - b;
        double bvirt = a - x;
        y = bvirt - b;
    }

    /**
     * \brief Computes the sum of a length 2 expansion and a double
     *  into a length 3 expansion.
     * \param[in] a1 high-magnitude component of first argument
     * \param[in] a0 low-magnitude component of first argument
     * \param[in] b second argument
     * \param[in] x2 high-magnitude component of the result
     * \param[in] x1 component of the result
     * \param[in] x0 low-magnitude component of the result
     * \details By Jonathan Shewchuk.
     */
    inline void two_one_sum(
        double a1, double a0, double b, double& x2, double& x1, double& x0
    ) {
        double _i;
        two_sum(a0, b, _i, x0);
        two_sum(a1, _i, x2, x1);
    }

    /**
     * \brief Computes the sum of a length 2 expansion and a double
     *  into a length 3 expansion.
     * \param[in] a1 high-magnitude component of first argument
     * \param[in] a0 low-magnitude component of first argument
     * \param[in] b1 high-magnitude component of second argument
     * \param[in] b0 high-magnitude component of second argument
     * \param[in] x3 high-magnitude component of the result
     * \param[in] x2 component of the result
     * \param[in] x1 component of the result
     * \param[in] x0 low-magnitude component of the result
     * \details By Jonathan Shewchuk.
     */
    inline void two_two_sum(
        double a1, double a0, double b1, double b0,
        double& x3, double& x2, double& x1, double& x0
    ) {
        double _j, _0;
        two_one_sum(a1, a0, b0, _j, _0, x0);
        two_one_sum(_j, _0, b1, x3, x2, x1);
    }

    /**
     * \brief Computes the product between two doubles where
     *  the second one have already been split.
     * \param[in] a first argument
     * \param[in] b second argument
     * \param[in] bhi high-magnitude part of second argument
     * \param[in] blo low-magnitude part of second argument
     * \param[out] x high-magnitude component of the result
     * \param[out] y low-magnitude component of the result
     * \details By Jonathan Shewchuk.
     */
    inline void two_product_presplit(
        double a, double b, double bhi, double blo, double& x, double& y
    ) {
        x = a * b;
        double ahi;
        double alo;
        split(a, ahi, alo);
        double err1 = x - (ahi * bhi);
        double err2 = err1 - (alo * bhi);
        double err3 = err2 - (ahi * blo);
        y = (alo * blo) - err3;
    }

    /**
     * \brief Computes the product between two doubles
     *  where both have already been split.
     * \param[in] a first argument
     * \param[in] ahi high-magnitude part of first argument
     * \param[in] alo low-magnitude part of first argument
     * \param[in] b second argument
     * \param[in] bhi high-magnitude part of second argument
     * \param[in] blo low-magnitude part of second argument
     * \param[out] x high-magnitude component of the result
     * \param[out] y low-magnitude component of the result
     * \details By Jonathan Shewchuk.
     */
    inline void two_product_2presplit(
        double a, double ahi, double alo,
        double b, double bhi, double blo,
        double& x, double& y
    ) {
        x = a * b;
        double err1 = x - (ahi * bhi);
        double err2 = err1 - (alo * bhi);
        double err3 = err2 - (ahi * blo);
        y = (alo * blo) - err3;
    }

    /**
     * \brief Computes the square of an expansion of length 2.
     * \param[in] a1 high-magnitude component of the argument
     * \param[in] a0 low-magnitude component of the argument
     * \param[out] x an array of six doubles to store the result.
     * \details By Jonathan Shewchuk.
     * An expansion of length two can be squared more quickly than finding the
     *  product of two different expansions of length two, and the result is
     *  guaranteed to have no more than six (rather than eight) components.
     */
    inline void two_square(
        double a1, double a0,
        double* x
    ) {
        double _0, _1, _2;
        double _j, _k, _l;
        square(a0, _j, x[0]);
        _0 = a0 + a0;
        two_product(a1, _0, _k, _1);
        two_one_sum(_k, _1, _j, _l, _2, x[1]);
        square(a1, _j, _1);
        two_two_sum(_j, _1, _l, _2, x[5], x[4], x[3], x[2]);
    }

    /**
     * \brief Computes the product of two expansions of length 2.
     * \param[in] a first argument (array of 2 doubles)
     * \param[in] b second argument (array of 2 doubles)
     * \param[out] x an array of 8 doubles to store the result
     * \details By Jonathan Shewchuk.
     */
    void two_two_product(
        const double* a,
        const double* b,
        double* x
    ) {

        double _0, _1, _2;
        double _i, _j, _k, _l, _m, _n;

        double a0hi, a0lo;
        split(a[0], a0hi, a0lo);
        double bhi, blo;
        split(b[0], bhi, blo);
        two_product_2presplit(
            a[0], a0hi, a0lo, b[0], bhi, blo, _i, x[0]
        );
        double a1hi, a1lo;
        split(a[1], a1hi, a1lo);
        two_product_2presplit(
            a[1], a1hi, a1lo, b[0], bhi, blo, _j, _0
        );
        two_sum(_i, _0, _k, _1);
        fast_two_sum(_j, _k, _l, _2);
        split(b[1], bhi, blo);
        two_product_2presplit(
            a[0], a0hi, a0lo, b[1], bhi, blo, _i, _0
        );
        two_sum(_1, _0, _k, x[1]);
        two_sum(_2, _k, _j, _1);
        two_sum(_l, _j, _m, _2);
        two_product_2presplit(
            a[1], a1hi, a1lo, b[1], bhi, blo, _j, _0
        );
        two_sum(_i, _0, _n, _0);
        two_sum(_1, _0, _i, x[2]);
        two_sum(_2, _i, _k, _1);
        two_sum(_m, _k, _l, _2);
        two_sum(_j, _n, _k, _0);
        two_sum(_1, _0, _j, x[3]);
        two_sum(_2, _j, _i, _1);
        two_sum(_l, _i, _m, _2);
        two_sum(_1, _k, _i, x[4]);
        two_sum(_2, _i, _k, x[5]);
        two_sum(_m, _k, x[7], x[6]);
    }

    /**
     * \brief Adds a scalar to an expansion, eliminating zero components
     *  from the output expansion.
     * \param[in] e first expansion
     * \param[in] b double to be added to \p e
     * \param[out] h the result \p e + \p b
     * \details Sets \p h = (\p e + \p b). \p e and \p h can be the same.
     *  This function is adapted from Jonathan Shewchuk's code.
     *  See the long version of his paper for details.
     *  Maintains the nonoverlapping property.  If round-to-even is used (as
     *  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent
     *  properties as well.  (That is, if e has one of these properties, so
     *  will h.)
     */
    void grow_expansion_zeroelim(
        const expansion& e, double b, expansion& h
    ) {
        double Q, hh;
        double Qnew;
        index_t eindex, hindex;
        index_t elen = e.length();

        hindex = 0;
        Q = b;
        for(eindex = 0; eindex < elen; eindex++) {
            double enow = e[eindex];
            two_sum(Q, enow, Qnew, hh);
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }

    /**
     * \brief Multiplies an expansion by a scalar,
     *  eliminating zero components from the
     *  output expansion.
     * \param[in] e an expansion
     * \param[in] b the double to be multiplied by \p e
     * \param[out] h the result \p b * \p e
     * \details (sets \p h = \p b * \p e). \p e and \p h cannot be the same.
     *  This function is adapted from Jonathan Shewchuk's code.
     *  See either version of his paper for details.
     *  Maintains the nonoverlapping property.  If round-to-even is used (as
     *  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent
     *  properties as well.  (That is, if e has one of these properties, so
     *  will h.)
     */
    void scale_expansion_zeroelim(
        const expansion& e, double b, expansion& h
    ) {
        double Q, sum;
        double hh;
        double product1;
        double product0;
        index_t eindex, hindex;
        double bhi, blo;
        index_t elen = e.length();

        // Sanity check: e and h cannot be the same.
        geo_debug_assert(&e != &h);

        split(b, bhi, blo);
        two_product_presplit(e[0], b, bhi, blo, Q, hh);
        hindex = 0;
        if(hh != 0) {
            h[hindex++] = hh;
        }
        for(eindex = 1; eindex < elen; eindex++) {
            double enow = e[eindex];
            two_product_presplit(enow, b, bhi, blo, product1, product0);
            two_sum(Q, product0, sum, hh);
            if(hh != 0) {
                h[hindex++] = hh;
            }
            fast_two_sum(product1, sum, Q, hh);
            if(hh != 0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }

    /**
     * \brief Sums two expansions, eliminating zero
     *  components from the output expansion (sets \p h = \p e + \p f).
     * \param[in] e the first expansion
     * \param[in] f the second expansion
     * \param[out] h the result \p e + \p f
     * \details h cannot be e or f.
     *  This function is adapted from Jonathan Shewchuk's code.
     *  See the long version of his paper for details.
     *  If round-to-even is used (as with IEEE 754), maintains the strongly
     *  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h
     *  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent
     *  properties.
     *
     */
    void fast_expansion_sum_zeroelim(
        const expansion& e, const expansion& f, expansion& h
    ) {
        double Q;
        double Qnew;
        double hh;
        index_t eindex, findex, hindex;
        double enow, fnow;
        index_t elen = e.length();
        index_t flen = f.length();

        // sanity check: h cannot be e or f
        geo_debug_assert(&h != &e);
        geo_debug_assert(&h != &f);

        enow = e[0];
        fnow = f[0];
        eindex = findex = 0;
        if((fnow > enow) == (fnow > -enow)) {
            Q = enow;
            enow = e[++eindex];
        } else {
            Q = fnow;
            fnow = f[++findex];
        }
        hindex = 0;
        if((eindex < elen) && (findex < flen)) {
            if((fnow > enow) == (fnow > -enow)) {
                fast_two_sum(enow, Q, Qnew, hh);
                enow = e[++eindex];
            } else {
                fast_two_sum(fnow, Q, Qnew, hh);
                fnow = f[++findex];
            }
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
            while((eindex < elen) && (findex < flen)) {
                if((fnow > enow) == (fnow > -enow)) {
                    two_sum(Q, enow, Qnew, hh);
                    enow = e[++eindex];
                } else {
                    two_sum(Q, fnow, Qnew, hh);
                    fnow = f[++findex];
                }
                Q = Qnew;
                if(hh != 0.0) {
                    h[hindex++] = hh;
                }
            }
        }
        while(eindex < elen) {
            two_sum(Q, enow, Qnew, hh);
            enow = e[++eindex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        while(findex < flen) {
            two_sum(Q, fnow, Qnew, hh);
            fnow = f[++findex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }

    /**
     * \brief Computes the difference of two expansions, eliminating zero
     *  components from the output expansion
     * \param[in] e first expansion
     * \param[in] f second expansion to be substracted from e
     * \param[out] h the result \p e - \p f
     * \details Sets \p h = (\p e - \p f). \p h cannot be \p e or \p f.
     *  This function is adapted from Jonathan Shewchuk's code.
     *  See the long version of his paper for details.
     *  If round-to-even is used (as with IEEE 754), maintains the strongly
     *  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h
     *  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent
     *  properties.
     */
    void fast_expansion_diff_zeroelim(
        const expansion& e, const expansion& f, expansion& h
    ) {
        double Q;
        double Qnew;
        double hh;
        index_t eindex, findex, hindex;
        double enow, fnow;
        index_t elen = e.length();
        index_t flen = f.length();

        // sanity check: h cannot be e or f
        geo_debug_assert(&h != &e);
        geo_debug_assert(&h != &f);

        enow = e[0];
        fnow = -f[0];
        eindex = findex = 0;
        if((fnow > enow) == (fnow > -enow)) {
            Q = enow;
            enow = e[++eindex];
        } else {
            Q = fnow;
            fnow = -f[++findex];
        }
        hindex = 0;
        if((eindex < elen) && (findex < flen)) {
            if((fnow > enow) == (fnow > -enow)) {
                fast_two_sum(enow, Q, Qnew, hh);
                enow = e[++eindex];
            } else {
                fast_two_sum(fnow, Q, Qnew, hh);
                fnow = -f[++findex];
            }
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
            while((eindex < elen) && (findex < flen)) {
                if((fnow > enow) == (fnow > -enow)) {
                    two_sum(Q, enow, Qnew, hh);
                    enow = e[++eindex];
                } else {
                    two_sum(Q, fnow, Qnew, hh);
                    fnow = -f[++findex];
                }
                Q = Qnew;
                if(hh != 0.0) {
                    h[hindex++] = hh;
                }
            }
        }
        while(eindex < elen) {
            two_sum(Q, enow, Qnew, hh);
            enow = e[++eindex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        while(findex < flen) {
            two_sum(Q, fnow, Qnew, hh);
            fnow = -f[++findex];
            Q = Qnew;
            if(hh != 0.0) {
                h[hindex++] = hh;
            }
        }
        if((Q != 0.0) || (hindex == 0)) {
            h[hindex++] = Q;
        }
        h.set_length(hindex);
    }
}

/****************************************************************************/

namespace GEO {

    double expansion_splitter_;
    double expansion_epsilon_;

    void expansion::initialize() {
        // Taken from Jonathan Shewchuk's exactinit.
        double half;
        double check, lastcheck;
        int every_other;

        every_other = 1;
        half = 0.5;
        expansion_epsilon_ = 1.0;
        expansion_splitter_ = 1.0;
        check = 1.0;
        // Repeatedly divide `epsilon' by two until it is too small to add to
        // one without causing roundoff.  (Also check if the sum is equal to
        // the previous sum, for machines that round up instead of using exact
        // rounding.  Not that this library will work on such machines anyway.
        do {
            lastcheck = check;
            expansion_epsilon_ *= half;
            if(every_other) {
                expansion_splitter_ *= 2.0;
            }
            every_other = !every_other;
            check = 1.0 + expansion_epsilon_;
        } while((check != 1.0) && (check != lastcheck));
        expansion_splitter_ += 1.0;
    }

    expansion* expansion::new_expansion_on_heap(index_t capa) {
        Memory::pointer addr = Memory::pointer(
            malloc(expansion::bytes(capa) + sizeof(index_t))
        );
        expansion* result = new(addr + sizeof(index_t))expansion(capa);
        expansion_refcount(result) = 0;
        return result;
    }

    void expansion::delete_expansion_on_heap(expansion* e) {
        free(expansion_baddr(e));
    }

    // ====== Initialization from expansion and double ===============

    expansion& expansion::assign_sum(const expansion& a, double b) {
        geo_debug_assert(capacity() >= sum_capacity(a, b));
        grow_expansion_zeroelim(a, b, *this);
        return *this;
    }

    expansion& expansion::assign_diff(const expansion& a, double b) {
        geo_debug_assert(capacity() >= diff_capacity(a, b));
        grow_expansion_zeroelim(a, -b, *this);
        return *this;
    }

    expansion& expansion::assign_product(const expansion& a, double b) {
        // TODO: implement special case where the double argument
        // is a power of two.
        geo_debug_assert(capacity() >= product_capacity(a, b));
        scale_expansion_zeroelim(a, b, *this);
        return *this;
    }

    // =============  expansion sum and difference =========================

    expansion& expansion::assign_sum(
        const expansion& a, const expansion& b
    ) {
        geo_debug_assert(capacity() >= sum_capacity(a, b));
        fast_expansion_sum_zeroelim(a, b, *this);
        return *this;
    }

    expansion& expansion::assign_sum(
        const expansion& a, const expansion& b, const expansion& c
    ) {
        geo_debug_assert(capacity() >= sum_capacity(a, b, c));
        expansion& ab = expansion_sum(a, b);
        this->assign_sum(ab, c);
        return *this;
    }

    expansion& expansion::assign_sum(
        const expansion& a, const expansion& b,
        const expansion& c, const expansion& d
    ) {
        geo_debug_assert(capacity() >= sum_capacity(a, b, c));
        expansion& ab = expansion_sum(a, b);
        expansion& cd = expansion_sum(c, d);
        this->assign_sum(ab, cd);
        return *this;
    }

    expansion& expansion::assign_diff(const expansion& a, const expansion& b) {
        geo_debug_assert(capacity() >= diff_capacity(a, b));
        fast_expansion_diff_zeroelim(a, b, *this);
        return *this;
    }

    // =============  expansion product ==================================

    // Recursive helper function for product implementation
    expansion& expansion::assign_sub_product(
        const double* a, index_t a_length, const expansion& b
    ) {
        geo_debug_assert(
            capacity() >= sub_product_capacity(a_length, b.length())
        );
        if(a_length == 1) {
            scale_expansion_zeroelim(b, a[0], *this);
        } else {
            // "Distillation" (see Shewchuk's paper) is computed recursively,
            // by splitting the list of expansions to sum into two halves.
            const double* a1 = a;
            index_t a1_length = a_length / 2;
            const double* a2 = a1 + a1_length;
            index_t a2_length = a_length - a1_length;
            expansion& a1b = expansion_sub_product(a1, a1_length, b);
            expansion& a2b = expansion_sub_product(a2, a2_length, b);
            this->assign_sum(a1b, a2b);
        }
        return *this;
    }

    expansion& expansion::assign_product(
        const expansion& a, const expansion& b
    ) {
        geo_debug_assert(capacity() >= product_capacity(a, b));
        if(a.length() == 0 || b.length() == 0) {
            x_[0] = 0.0;
            set_length(0);
        } else if(a.length() == 1 && b.length() == 1) {
            two_product(a[0], b[0], x_[1], x_[0]);
            set_length(2);
        } else if(a.length() == 1) {
            scale_expansion_zeroelim(b, a[0], *this);
        } else if(b.length() == 1) {
            scale_expansion_zeroelim(a, b[0], *this);
        } else if(a.length() == 2 && b.length() == 2) {
            two_two_product(a.data(), b.data(), x_);
            set_length(8);
        } else {
            // TODO: should we split the shortest or the longest
            // expansion ? 
            // Recursive distillation: the shortest expansion
            // is split into two parts.
            if(a.length() < b.length()) {
                const double* a1 = a.data();
                index_t a1_length = a.length() / 2;
                const double* a2 = a1 + a1_length;
                index_t a2_length = a.length() - a1_length;
                expansion& a1b = expansion_sub_product(a1, a1_length, b);
                expansion& a2b = expansion_sub_product(a2, a2_length, b);
                this->assign_sum(a1b, a2b);
            } else {
                const double* b1 = b.data();
                index_t b1_length = b.length() / 2;
                const double* b2 = b1 + b1_length;
                index_t b2_length = b.length() - b1_length;
                expansion& ab1 = expansion_sub_product(b1, b1_length, a);
                expansion& ab2 = expansion_sub_product(b2, b2_length, a);
                this->assign_sum(ab1, ab2);
            }
        }
        return *this;
    }

    expansion& expansion::assign_product(
        const expansion& a, const expansion& b, const expansion& c
    ) {
        const expansion& bc = expansion_product(b, c);
        this->assign_product(a, bc);
        return *this;
    }

    expansion& expansion::assign_square(const expansion& a) {
        geo_debug_assert(capacity() >= square_capacity(a));
        if(a.length() == 1) {
            square(a[0], x_[1], x_[0]);
            set_length(2);
        } else if(a.length() == 2) {
            two_square(a[1], a[0], x_);
            set_length(6);
        } else {
            this->assign_product(a, a);
        }
        return *this;
    }

    // =============  determinants ==========================================

    expansion& expansion::assign_det2x2(
        const expansion& a11, const expansion& a12,
        const expansion& a21, const expansion& a22
    ) {
        const expansion& a11a22 = expansion_product(a11, a22);
        const expansion& a12a21 = expansion_product(a12, a21);
        return this->assign_diff(a11a22, a12a21);
    }

    expansion& expansion::assign_det3x3(
        const expansion& a11, const expansion& a12, const expansion& a13,
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    ) {
        // Development w.r.t. first row
        const expansion& c11 = expansion_det2x2(a22, a23, a32, a33);
        const expansion& c12 = expansion_det2x2(a23, a21, a33, a31);
        const expansion& c13 = expansion_det2x2(a21, a22, a31, a32);
        const expansion& a11c11 = expansion_product(a11, c11);
        const expansion& a12c12 = expansion_product(a12, c12);
        const expansion& a13c13 = expansion_product(a13, c13);
        return this->assign_sum(a11c11, a12c12, a13c13);
    }

    expansion& expansion::assign_det_111_2x3(
        const expansion& a21, const expansion& a22, const expansion& a23,
        const expansion& a31, const expansion& a32, const expansion& a33
    ) {
        const expansion& c11 = expansion_det2x2(a22, a23, a32, a33);
        const expansion& c12 = expansion_det2x2(a23, a21, a33, a31);
        const expansion& c13 = expansion_det2x2(a21, a22, a31, a32);
        return this->assign_sum(c11, c12, c13);
    }

    // =============  geometric operations ==================================

    expansion& expansion::assign_sq_dist(
        const double* p1, const double* p2, coord_index_t dim
    ) {
        geo_debug_assert(capacity() >= sq_dist_capacity(dim));
        if(dim == 1) {
            double d0, d1;
            two_diff(p1[0], p2[0], d1, d0);
            two_square(d1, d0, x_);
            set_length(6);
        } else {
            // "Distillation" (see Shewchuk's paper) is computed recursively,
            // by splitting the list of expansions to sum into two halves.
            coord_index_t dim1 = dim / 2;
            coord_index_t dim2 = coord_index_t(dim - dim1);
            const double* p1_2 = p1 + dim1;
            const double* p2_2 = p2 + dim1;
            expansion& d1 = expansion_sq_dist(p1, p2, dim1);
            expansion& d2 = expansion_sq_dist(p1_2, p2_2, dim2);
            this->assign_sum(d1, d2);
        }
        return *this;
    }

    expansion& expansion::assign_dot_at(
        const double* p1, const double* p2, const double* p0,
        coord_index_t dim
    ) {
        geo_debug_assert(capacity() >= dot_at_capacity(dim));
        if(dim == 1) {
            double v[2];
            two_diff(p1[0], p0[0], v[1], v[0]);
            double w[2];
            two_diff(p2[0], p0[0], w[1], w[0]);
            two_two_product(v, w, x_);
            set_length(8);
        } else {
            // "Distillation" (see Shewchuk's paper) is computed recursively,
            // by splitting the list of expansions to sum into two halves.
            coord_index_t dim1 = dim / 2;
            coord_index_t dim2 = coord_index_t(dim - dim1);
            const double* p1_2 = p1 + dim1;
            const double* p2_2 = p2 + dim1;
            const double* p0_2 = p0 + dim1;
            expansion& d1 = expansion_dot_at(p1, p2, p0, dim1);
            expansion& d2 = expansion_dot_at(p1_2, p2_2, p0_2, dim2);
            this->assign_sum(d1, d2);
        }
        return *this;
    }

    /************************************************************************/

    expansion_nt& expansion_nt::operator+= (const expansion_nt& rhs) {
        index_t e_capa = expansion::sum_capacity(rep(), rhs.rep());
        expansion* e = expansion::new_expansion_on_heap(e_capa);
        e->assign_sum(rep(), rhs.rep());
        expansion::unref_expansion(rep_);
        rep_ = e;
        expansion::ref_expansion(rep_);
        return *this;
    }

    expansion_nt& expansion_nt::operator+= (double rhs) {
        index_t e_capa = expansion::sum_capacity(rep(), rhs);

        // TODO: optimized in-place version to be used
        //   if(!shared() && e_capa < rep().capacity())

        expansion* e = expansion::new_expansion_on_heap(e_capa);
        e->assign_sum(rep(), rhs);
        expansion::unref_expansion(rep_);
        rep_ = e;
        expansion::ref_expansion(rep_);

        return *this;
    }

    expansion_nt& expansion_nt::operator-= (const expansion_nt& rhs) {
        index_t e_capa = expansion::diff_capacity(rep(), rhs.rep());
        expansion* e = expansion::new_expansion_on_heap(e_capa);
        e->assign_diff(rep(), rhs.rep());
        expansion::unref_expansion(rep_);
        rep_ = e;
        expansion::ref_expansion(rep_);
        return *this;
    }

    expansion_nt& expansion_nt::operator-= (double rhs) {
        index_t e_capa = expansion::diff_capacity(rep(), rhs);

        // TODO: optimized in-place version to be used
        //   if(!shared() && e_capa < rep().capacity())

        expansion* e = expansion::new_expansion_on_heap(e_capa);
        e->assign_diff(rep(), rhs);
        expansion::unref_expansion(rep_);
        rep_ = e;
        expansion::ref_expansion(rep_);

        return *this;
    }

    expansion_nt& expansion_nt::operator*= (const expansion_nt& rhs) {
        index_t e_capa = expansion::product_capacity(rep(), rhs.rep());
        expansion* e = expansion::new_expansion_on_heap(e_capa);
        e->assign_product(rep(), rhs.rep());
        expansion::unref_expansion(rep_);
        rep_ = e;
        expansion::ref_expansion(rep_);
        return *this;
    }

    expansion_nt& expansion_nt::operator*= (double rhs) {
        index_t e_capa = expansion::product_capacity(rep(), rhs);

        // TODO: optimized in-place version to be used
        //   if(!shared() && e_capa < rep().capacity())

        expansion* e = expansion::new_expansion_on_heap(e_capa);
        e->assign_product(rep(), rhs);
        expansion::unref_expansion(rep_);
        rep_ = e;
        expansion::ref_expansion(rep_);
        return *this;
    }

    /************************************************************************/

    expansion_nt expansion_nt::operator+ (const expansion_nt& rhs) const {
        expansion* e = expansion::new_expansion_on_heap(
            expansion::sum_capacity(rep(), rhs.rep())
        );
        e->assign_sum(rep(), rhs.rep());
        return expansion_nt(e);
    }

    expansion_nt expansion_nt::operator- (const expansion_nt& rhs) const {
        expansion* e = expansion::new_expansion_on_heap(
            expansion::diff_capacity(rep(), rhs.rep())
        );
        e->assign_diff(rep(), rhs.rep());
        return expansion_nt(e);
    }

    expansion_nt expansion_nt::operator* (const expansion_nt& rhs) const {
        expansion* e = expansion::new_expansion_on_heap(
            expansion::product_capacity(rep(), rhs.rep())
        );
        e->assign_product(rep(), rhs.rep());
        return expansion_nt(e);
    }

    expansion_nt expansion_nt::operator+ (double rhs) const {
        expansion* e = expansion::new_expansion_on_heap(
            expansion::sum_capacity(rep(), rhs)
        );
        e->assign_sum(rep(), rhs);
        return expansion_nt(e);
    }

    expansion_nt expansion_nt::operator- (double rhs) const {
        expansion* e = expansion::new_expansion_on_heap(
            expansion::diff_capacity(rep(), rhs)
        );
        e->assign_diff(rep(), rhs);
        return expansion_nt(e);
    }

    expansion_nt expansion_nt::operator* (double rhs) const {
        expansion* e = expansion::new_expansion_on_heap(
            expansion::product_capacity(rep(), rhs)
        );
        e->assign_product(rep(), rhs);
        return expansion_nt(e);
    }

    /************************************************************************/

    expansion_nt expansion_nt::operator- () const {
        expansion_nt result(*this);
        result.rep().negate();
        return result;
    }

    /************************************************************************/

    bool expansion_nt::operator> (const expansion_nt& rhs) const {
        const expansion& diff = expansion_diff(rep(), rhs.rep());
        return diff.sign() > 0;
    }

    bool expansion_nt::operator>= (const expansion_nt& rhs) const {
        const expansion& diff = expansion_diff(rep(), rhs.rep());
        return diff.sign() >= 0;
    }

    bool expansion_nt::operator< (const expansion_nt& rhs) const {
        const expansion& diff = expansion_diff(rep(), rhs.rep());
        return diff.sign() < 0;
    }

    bool expansion_nt::operator<= (const expansion_nt& rhs) const {
        const expansion& diff = expansion_diff(rep(), rhs.rep());
        return diff.sign() <= 0;
    }

    bool expansion_nt::operator> (double rhs) const {
        if(rhs == 0.0) {
            return rep().sign() > 0;
        }
        const expansion& diff = expansion_diff(rep(), rhs);
        return diff.sign() > 0;
    }

    bool expansion_nt::operator>= (double rhs) const {
        if(rhs == 0.0) {
            return rep().sign() >= 0;
        }
        const expansion& diff = expansion_diff(rep(), rhs);
        return diff.sign() >= 0;
    }

    bool expansion_nt::operator< (double rhs) const {
        if(rhs == 0.0) {
            return rep().sign() < 0;
        }
        const expansion& diff = expansion_diff(rep(), rhs);
        return diff.sign() < 0;
    }

    bool expansion_nt::operator<= (double rhs) const {
        if(rhs == 0.0) {
            return rep().sign() <= 0;
        }
        const expansion& diff = expansion_diff(rep(), rhs);
        return diff.sign() <= 0;
    }

    /************************************************************************/
}

