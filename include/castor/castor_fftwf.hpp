/*

+-----------------------------------------------------------+
|                                                           |
|  (c) 2020 - PROPERTY OF MARC BAKRY - MIT license          |
|                                                           |
+-----------------------------------------------------------+
|                                                           |
|       FILE       : castor_fftwf.hpp                       |
|       VERSION    : 0.1.0                                  |
|       AUTHOR     : Marc Bakry                             |
|       CREATION   : 29.11.2020                             |
|       LAST MODIF : 29.11.2020                             |
|       SYNOPSIS   : Wraps the well-known FFTW3 C library   |
|                    within the 'castor' framework for      |
|                    single precision accuracy.             |
|                                                           |
+-----------------------------------------------------------+

 */
#ifndef CASTOR_FFTWF_HPP
#define CASTOR_FFTWF_HPP

#include <algorithm>
extern "C"{
    #include <fftw3.h>
}

#include "castor/matrix.hpp"


namespace castor{


namespace fftw{

//////////////////////////////////////////
// DECLARATIONS OF THE PUBLIC INTERFACE //
//////////////////////////////////////////
matrix<std::complex<float>> fft(matrix<float> &X, int dim=1);
matrix<std::complex<float>> fft(matrix<std::complex<float>> &X, int dim=1);
matrix<std::complex<float>> ifft(matrix<float> &X, int dim=1);
matrix<std::complex<float>> ifft(matrix<std::complex<float>> &X, int dim=1);
matrix<std::complex<float>> fft2(matrix<float> &X);
matrix<std::complex<float>> fft2(matrix<std::complex<float>> &X);
matrix<std::complex<float>> fft2(std::size_t m, std::size_t n, matrix<float> &X);
matrix<std::complex<float>> fft2(std::size_t m, std::size_t n, matrix<std::complex<float>> &X);
matrix<std::complex<float>> ifft2(matrix<float> &X);
matrix<std::complex<float>> ifft2(matrix<std::complex<float>> &X);
matrix<std::complex<float>> ifft2(std::size_t m, std::size_t n, matrix<float> &X);
matrix<std::complex<float>> ifft2(std::size_t m, std::size_t n, matrix<std::complex<float>> &X);

matrix<float> fftfreq(std::size_t n, float d = 1.);

matrix<std::complex<float>> fftshift(std::complex<float> const &A, int dim=-1);
matrix<std::complex<float>> ifftshift(std::complex<float> const &A, int dim=-1);


//////////////////////////////
// NON-DOCUMENTED INTERFACE //
//////////////////////////////

// interface to forward/backward 1d FFT depending on sign
matrix<std::complex<float>> xfft(matrix<std::complex<float>> &X, int sign, int dim=1)
{
    auto m = size(X,1); // nb. of lines
    auto n = size(X,2); // nb. of columns
    bool isMatrix = (m != 1 && n != 1);
    matrix<std::complex<float>> Y = zeros<std::complex<float>>(m,n);
    fftwf_complex *in = reinterpret_cast<fftwf_complex*>(&X(0));
    fftwf_complex *out = reinterpret_cast<fftwf_complex*>(&Y(0));
    if(!isMatrix) // vector case
    {
        fftwf_plan plan = fftwf_plan_dft_1d(std::max(m,n), in, out, sign, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
    }
    else
    {   
        // default values : FFT column by column
        int rank    = 1; 
        int howmany = n;
        int istride = n;
        int idist   = 1;
        auto N = std::vector<int>();
        if(dim == 1)
        {
            N = std::vector<int>(howmany,m);
        }
        else if(dim == 2) // if FFT line by line
        {
            howmany = m; istride = 1; idist = n; N = std::vector<int>(howmany,n);
        }
        else
        {
            throw;
        }
        fftwf_plan plan = fftwf_plan_many_dft(rank,&N[0],howmany,in,nullptr,istride,idist,out,nullptr,istride,idist,sign,FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
    }
    return Y;
}


// interface to forward/backward 2d FFT depending on sign
matrix<std::complex<float>> xfft2(matrix<std::complex<float>> &X, int sign)
{
    auto m = size(X,1); // nb. of lines
    auto n = size(X,2); // nb. of columns
    matrix<std::complex<float>> Y = zeros<std::complex<float>>(m,n);
    fftwf_complex *in  = reinterpret_cast<fftwf_complex*>(&X(0));
    fftwf_complex *out = reinterpret_cast<fftwf_complex*>(&Y(0));
    //
    fftwf_plan plan  = fftwf_plan_dft_2d(m,n,in,out,sign,FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    return Y;
}

// interface to forward/backward 2d FFT depending on sign
matrix<std::complex<float>> xfft2(std::size_t m, std::size_t n, matrix<std::complex<float>> &X, int sign)
{
    // check consistency
    if(m*n != size(X,1)*size(X,2)) error(__FILE__, __LINE__, __FUNCTION__,"Size of input given as parameter is not consistent with the number of elements in the input matrix.");
    matrix<std::complex<float>> Y = zeros<std::complex<float>>(size(X,1),size(X,2));
    fftwf_complex *in  = reinterpret_cast<fftwf_complex*>(&X(0));
    fftwf_complex *out = reinterpret_cast<fftwf_complex*>(&Y(0));
    //
    fftwf_plan plan  = fftwf_plan_dft_2d(static_cast<int>(m),static_cast<int>(n),in,out,sign,FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    return Y;
}

}

matrix<std::complex<float>> fftw::fft(matrix<std::complex<float>> &X, int dim)
{
    return fftw::xfft(X,FFTW_FORWARD,dim);
}
matrix<std::complex<float>> fftw::fft(matrix<float> &X, int dim)
{
    matrix<std::complex<float>> XX = cast<std::complex<float>>(X);
    return fftw::xfft(XX,FFTW_FORWARD,dim);
}
matrix<std::complex<float>> fftw::ifft(matrix<std::complex<float>> &X, int dim)
{
    return fftw::xfft(X,FFTW_BACKWARD,dim);
}
matrix<std::complex<float>> fftw::ifft(matrix<float> &X, int dim)
{
    matrix<std::complex<float>> XX = cast<std::complex<float>>(X);
    return fftw::xfft(XX,FFTW_BACKWARD,dim);
}
matrix<std::complex<float>> fftw::fft2(matrix<std::complex<float>> &X)
{
    return fftw::xfft2(X,FFTW_FORWARD);
}
matrix<std::complex<float>> fftw::fft2(matrix<float> &X)
{
    matrix<std::complex<float>> XX = cast<std::complex<float>>(X);
    return fftw::xfft2(XX,FFTW_FORWARD);
}
matrix<std::complex<float>> fftw::fft2(std::size_t m, std::size_t n, matrix<std::complex<float>> &X)
{
    return fftw::xfft2(m,n,X,FFTW_FORWARD);
}
matrix<std::complex<float>> fftw::fft2(std::size_t m, std::size_t n, matrix<float> &X)
{
    matrix<std::complex<float>> XX = cast<std::complex<float>>(X);
    return fftw::xfft2(m,n,XX,FFTW_FORWARD);
}
matrix<std::complex<float>> fftw::ifft2(matrix<std::complex<float>> &X)
{
    return fftw::xfft2(X,FFTW_BACKWARD);
}
matrix<std::complex<float>> fftw::ifft2(matrix<float> &X)
{
    matrix<std::complex<float>> XX = cast<std::complex<float>>(X);
    return fftw::xfft2(XX,FFTW_BACKWARD);
}
matrix<std::complex<float>> fftw::ifft2(std::size_t m, std::size_t n, matrix<std::complex<float>> &X)
{
    return fftw::xfft2(m,n,X,FFTW_BACKWARD);
}
matrix<std::complex<float>> fftw::ifft2(std::size_t m, std::size_t n, matrix<float> &X)
{
    matrix<std::complex<float>> XX = cast<std::complex<float>>(X);
    return fftw::xfft2(m,n,XX,FFTW_BACKWARD);
}



matrix<float> fftw::fftfreq(std::size_t n, float d)
{
    matrix<float> freq = zeros(1,n);
    if(n%2 == 0)
    {
        for(std::size_t i=0; i<n/2; ++i) freq(i)     = i/(d*n);
        for(std::size_t i=0; i<n/2; ++i) freq(n-1-i) = -(1+i)/(d*n);
    }
    else
    {
        for(std::size_t i=0; i<(n-1)/2+1; ++i) freq(i) = i/(d*n);
        for(std::size_t i=0; i<(n-1)/2; ++i)   freq(i) = -(1+i)/(d*n);
    }
    return freq;
}


matrix<std::complex<float>> fftshift(matrix<std::complex<float>> const &A, int dim)
{
    std::size_t m=size(A,1), n=size(A,2);
    matrix<std::complex<float>> As = zeros(m,n);

    //
    if(m == 1 || n == 1) // one-dimensional case
    {
        std::size_t mn  = m*n;
        std::size_t cnt = 0;
        if(mn%2 == 0)
        {
            for(std::size_t i=n/2; i<n; ++i) As(cnt++) = A(i);
            for(std::size_t i=0; i<n/2; ++i) As(cnt++) = A(i);
        }
        else
        {
            // 
        }
    }
    else
    {
        //
    }
    //
    return As;
}

// END OF NAMESPACE
}

#endif // CASTOR_FFTWF_HPP