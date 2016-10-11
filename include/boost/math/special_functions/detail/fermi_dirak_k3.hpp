//  Copyright (c) 2016 Semyon Kolganov
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FERMI_DIRAK_K3_HPP
#define BOOST_MATH_FERMI_DIRAK_K3_HPP

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

template <typename T>
T fermi_dirak_k3(const T x)
{
    const int k = 3;
    const T a[] = { 
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.1583482145380455955096383)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0460645149909308107878344)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0048861379108841469134267)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0004336733305971515517559)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0000173435613795895152436)) 
	};

    const T b[] = { 
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0125148812047107612191739)),
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0266693407000929631393759)),
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0003285431094547362504004)),
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0000820910787890062715299)) 
	};

    if (x < 0) {
        return get_negative_value(x, k, a, b);
    }
    else {
        return pow(x, 4) / 4 + 6 * x * x * I_1_0 + 2 * I_3_0 - 
               get_negative_value(-x, k, a, b);
    }
    
}

}}} // namespaces

#endif // BOOST_MATH_FERMI_DIRAK_K3_HPP