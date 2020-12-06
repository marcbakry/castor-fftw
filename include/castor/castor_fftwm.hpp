
#ifndef CASTOR_FFTWM_HPP
#define CASTOR_FFTWM_HPP

#include "castor/matrix.hpp"

namespace castor{
namespace fftw{

template<typename T>
matrix<T> fftfreq(std::size_t n, T d)
{
    matrix<T> freq = zeros(1,n);
    if(n%2 == 0)
    {
        for(std::size_t i=0; i<n/2; ++i) freq(i)     = i/(d*n);
        for(std::size_t i=0; i<n/2; ++i) freq(n-1-i) = -(1.+i)/(d*n);
    }
    else
    {
        for(std::size_t i=0; i<(n-1)/2+1; ++i) 
        {
            freq(i) = i/(d*n);
        }
        for(std::size_t i=0; i<(n-1)/2; ++i)   
        {
            freq(n-1-i) = -(1.+i)/(d*n);
        }
    }
    return freq;
}


template<typename T>
matrix<T> fftshift(matrix<T> const &A, int dim = 0)
{
    std::size_t m=size(A,1), n=size(A,2);
    matrix<T> As = zeros(m,n);

    //
    if(m == 1 || n == 1) // one-dimensional case
    {
        std::size_t mn  = m*n;
        std::size_t cnt = 0;
        if(mn%2 == 0)
        {
            for(std::size_t i=mn/2; i<mn; ++i) As(cnt++) = A(i);
            for(std::size_t i=0; i<mn/2; ++i)  As(cnt++) = A(i);
        }
        else
        {
            for(std::size_t i=(mn-1)/2+1; i<mn; ++i) As(cnt++) = A(i);
            for(std::size_t i=0; i<(mn-1)/2+1; ++i)  As(cnt++) = A(i);
        }
    }
    else
    {
        //
    }
    //
    return As;
}


template<typename T>
matrix<T> ifftshift(matrix<T> const &As, int dim = 0)
{
    std::size_t m=size(As,1), n=size(As,2);
    matrix<T> A= zeros(m,n);

    //
    if(m == 1 || n == 1)
    {
        std::size_t mn  = m*n;
        std::size_t cnt = 0;
        if(mn%2 == 0)
        {
            for(std::size_t i=mn/2; i<mn; ++i) A(cnt++) = As(i);
            for(std::size_t i=0; i<mn/2; ++i) A(cnt++) = As(i);
        }
        else
        {
            for(std::size_t i=(mn-1)/2; i<mn; ++i) A(cnt++) = As(i);
            for(std::size_t i=0; i<(mn-1)/2; ++i) A(cnt++) = As(i);
        }
    }
    else
    {
        /* code */
    }
    return A;
}

        
}
}

#endif