//  Copyright (c) 2016 Semyon Kolganov
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FERMI_DIRAK_K1_HPP
#define BOOST_MATH_FERMI_DIRAK_K1_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/rational.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/assert.hpp>

// New quick metod of calculation Fermi-Dirak function 
// with high precision

namespace boost { namespace math { namespace detail {
	
// unified function for calculation fermi-dirak function of integer index
// at negative x < 0
template <typename T>
static T evaluate_negative_part(const T x, const T k,
                            const T& a, const T& b)
{
    const size_t N_base = 4;
    T S1 = 0, S2 = 0;
    T y = log(1 + exp(x));

    for (size_t n = 0; n < N_base + 1; n++) {
        S1 += a[n] * pow(y, n + 1);
    }

    for (size_t m = 0; m < N_base; m++) {
        S2 += b[m] * pow(y, m + 1);
    }

    T I_negative = (1 + S1) / (1 + S2);

    // return value of approximation according (14) in 2016, N. Kalitkin, S. Kolganov
    return boost::math::tgamma(k + 1) * y * pow(I_negative, k);
}
	
template <typename T>
T fermi_dirak_k1(const T x)
{
    const int k = 1;
    static const T a[] = { 
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.2715113138214362780964488)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0562661238060587633169245)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0067420740469345689743873)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0005169505155333205859985)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0000194771836765773190602)) 
	};

    static const T b[] = { 
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0215113138214352840651584)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0231105175729721417901084)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0003669081577365413477999)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0000610424408732720110769)) 
	};

    if (x < 0) {
        return evaluate_negative_part(x, k, a, b);
    }
    else {
        return x * x / 2 + 2 * I_1_0 - evaluate_negative_part(-x, k, a, b);
    }
}

}}} // namespaces

#endif // BOOST_MATH_FERMI_DIRAK_K1_HPP