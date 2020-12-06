
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
matrix<T> fftshift_1d(matrix<T> const &A)
{
    std::size_t m=size(A,1), n=size(A,2);
    matrix<T> As = zeros(m,n);
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
    return As;
}

template<typename T>
matrix<T> fftshift_2d1d(matrix<T> const &A, int dim=1)
{
    std::size_t m = size(A,1), n = size(A,2);
    matrix<T> As = zeros(m,n);
    if(dim == 1)
    {
        if(m%2 == 0)
        {
            for(std::size_t j=0; j<n; ++j)
            {
                std::size_t cnt = 0;
                for(std::size_t i=m/2; i<m; ++i) As(cnt++,j) = A(i,j);
                for(std::size_t i=0; i<m/2; ++i) As(cnt++,j) = A(i,j);
            }
        }
        else
        {
            for(std::size_t j=0; j<n; ++j)
            {
                std::size_t cnt = 0;
                for(std::size_t i=(m-1)/2+1; i<m; ++i) As(cnt++,j) = A(i,j);
                for(std::size_t i=0; i<(m-1)/2+1; ++i) As(cnt++,j) = A(i,j);
            }
        }
    }
    else if(dim == 2)
    {
        if(n%2 == 0)
        {
            for(std::size_t j=0; j<m; ++j)
            {
                std::size_t cnt = 0;
                for(std::size_t i=n/2; i<n; ++i) As(j,cnt++) = A(j,i);
                for(std::size_t i=0; i<n/2; ++i) As(j,cnt++) = A(j,i);
            }
        }
        else
        {
            for(std::size_t j=0; j<m; ++j)
            {
                std::size_t cnt = 0;
                for(std::size_t i=(n-1)/2+1; i<n; ++i) As(j,cnt++) = A(j,i);
                for(std::size_t i=0; i<(n-1)/2+1; ++i) As(j,cnt++) = A(j,i);
            }
        }
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"'dim' parameter should be 1 or 2.");
    }
    return As;
}


template<typename T>
matrix<T> ifftshift_1d(matrix<T> const &As)
{
    std::size_t m=size(As,1), n=size(As,2);
    matrix<T> A = zeros(m,n);
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
    return A;
}

template<typename T>
matrix<T> ifftshift_2d1d(matrix<T> const &As, int dim = 0)
{
    std::size_t m = size(As,1), n = size(As,2);
    matrix<T> A = zeros(m,n);
    if(dim == 1)
    {
        if(m%2 == 0)
        {
            for(std::size_t j=0; j<n; ++j)
            {
                std::size_t cnt = 0;
                for(std::size_t i=m/2; i<m; ++i) A(cnt++,j) = As(i,j);
                for(std::size_t i=0; i<m/2; ++i) A(cnt++,j) = As(i,j);
            }
        }
        else
        {
            for(std::size_t j=0; j<n; ++j)
            {
                std::size_t cnt = 0;
                for(std::size_t i=(m-1)/2; i<m; ++i) A(cnt++,j) = As(i,j);
                for(std::size_t i=0; i<(m-1)/2; ++i) A(cnt++,j) = As(i,j);
            }
        }
    }
    else if(dim == 2)
    {
        if(n%2 == 0)
        {
            for(std::size_t j=0; j<m; ++j)
            {
                std::size_t cnt = 0;
                for(std::size_t i=n/2; i<n; ++i) A(j,cnt++) = As(j,i);
                for(std::size_t i=0; i<n/2; ++i) A(j,cnt++) = As(j,i);
            }
        }
        else
        {
            for(std::size_t j=0; j<m; ++j)
            {
                std::size_t cnt = 0;
                for(std::size_t i=(n-1)/2; i<n; ++i) A(j,cnt++) = As(j,i);
                for(std::size_t i=0; i<(n-1)/2; ++i) A(j,cnt++) = As(j,i);
            }
        }
    }
    else
    {
        error(__FILE__, __LINE__, __FUNCTION__,"'dim' parameter should be 1 or 2.");
    }
    return A;
}





// PUBLIC INTERFACE
template<typename T>
matrix<T> fftshift(matrix<T> const &A, int dim = 0)
{
    std::size_t m=size(A,1), n=size(A,2);
    //
    if(m == 1 || n == 1) // one-dimensional case
    {
        return fftshift_1d(A);
    }
    else if(dim != 0)
    {
        return fftshift_2d1d(A,dim);
    }
    //
    return fftshift_2d1d(fftshift_2d1d(A,1),2);
}


template<typename T>
matrix<T> ifftshift(matrix<T> const &As, int dim = 0)
{
    std::size_t m=size(As,1), n=size(As,2);
    //
    if(m == 1 || n == 1)
    {
        return ifftshift_1d(As);
    }
    else if(dim != 0)
    {
        return ifftshift_2d1d(As,dim);
    }
    return ifftshift_2d1d(ifftshift_2d1d(As,2),1);
}


// END OF NAMESPACE        
}
}

#endif