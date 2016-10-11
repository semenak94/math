//  Copyright (c) 2016 Semyon Kolganov
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FERMI_DIRAK_K2_HPP
#define BOOST_MATH_FERMI_DIRAK_K2_HPP

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
T fermi_dirak_k2(const T x)
{
    const int k = 2;
    const T a[] = { 
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.2263816364340698560028783)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0533684335574798857246766)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0062904756340795211604491)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0005023228274452983506998)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0000189379675088061004880))
	};

    const T b[] = { 
		static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0388816364340691133155655)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0243043998742774445085992)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0006290985326433190105734)),
        static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0000657018161945458806177)) 
	};

    if (x < 0) {
        return get_negative_value(x, k, a, b);
    }
    else {
        return pow(x, 3) / 3 + 4 * x * I_1_0 + 
               get_negative_value(-x, k, a, b);
    }
}

}}} // namespaces

#endif // BOOST_MATH_FERMI_DIRAK_K2_HPP