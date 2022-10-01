#ifndef SIMPLE_POLY_HPP 
#define SIMPLE_POLY_HPP 

#include <vector> 
#include <iostream>

template<typename _Scalar, int deg> 
struct Polynomial {
    Polynomial() : coef(deg + 1, _Scalar(0)) {} 
    Polynomial(_Scalar const (&coef)[deg+1]) : coef(deg + 1, _Scalar(0))  {
        for(int i = 0; i <= deg; ++i) {
            this -> coef[i] = coef[i];
        }
    }  

    _Scalar& operator[] (size_t idx) {
        return this -> coef[idx]; 
    }

    _Scalar const & operator[] (size_t idx) const {
        return this -> coef[idx]; 
    }

    _Scalar evaluate(_Scalar const &x) const { 
        _Scalar ret = this -> coef[deg];
        for(int i = deg - 1; i >= 0; --i) {
            ret = x * ret + this -> coef[i];
        }
        return ret; 
    }

    friend std::ostream& operator<<(std::ostream& os, Polynomial const &p) {
        for(auto c: p.coef) { 
            os << c << ", "; 
        }
        return os;
    }

    std::vector<_Scalar> coef; 
};

template<typename _Scalar, int deg1, int deg2> 
Polynomial<_Scalar,deg1 + deg2> operator* (Polynomial<_Scalar, deg1> const &p1, Polynomial<_Scalar, deg2> const &p2) {
    Polynomial<_Scalar,deg1 + deg2> ret; 
    for(int i = 0; i <= deg1; ++i) {
        for(int j = 0; j <= deg2; ++j) {
            ret[i + j] += p1[i] * p2[j]; 
        }
    }
    return ret; 
}

template<typename _Scalar, int deg> 
Polynomial<_Scalar,deg> operator* (_Scalar c, Polynomial<_Scalar, deg> const &p1) {
    Polynomial<_Scalar,deg> ret; 
    for(int i = 0; i <= deg; ++i) {
        ret[i] += c * p1[i]; 
    }
    return ret; 
}

template<typename _Scalar, int deg1, int deg2> 
Polynomial<_Scalar, std::max(deg1, deg2)> operator+(Polynomial<_Scalar, deg1> const &p1, Polynomial<_Scalar, deg2> const &p2) {
    if(deg1 > deg2) {
        return p2 + p1;
    }
    Polynomial<_Scalar, std::max(deg1, deg2)> ret; 
    int i; 
    for(i = 0; i <= deg1; ++i) 
        ret[i] = p1[i] + p2[i];
    for(; i <= deg2; ++i) 
        ret[i] = p2[i];
    return ret; 
}

template<typename _Scalar, int deg> 
Polynomial<_Scalar, deg> make_polynomial(_Scalar const (&coef)[deg+1]) {
    return Polynomial<_Scalar,deg>(coef);
}
#endif 