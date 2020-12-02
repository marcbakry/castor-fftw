/*

+-----------------------------------------------------------+
|                                                           |
|  (c) 2020 - PROPERTY OF MARC BAKRY - MIT license          |
|                                                           |
+-----------------------------------------------------------+
|                                                           |
|       FILE       : castor_fftw.hpp                        |
|       VERSION    : 0.1.0                                  |
|       AUTHOR     : Marc Bakry                             |
|       CREATION   : 29.11.2020                             |
|       LAST MODIF : 29.11.2020                             |
|       SYNOPSIS   : Wraps the well-known FFTW3 C library   |
|                    within the 'castor' framework.         |
|                                                           |
+-----------------------------------------------------------+

 */
#ifndef CASTOR_FFTW_HPP
#define CASTOR_FFTW_HPP

#include <algorithm>
extern "C"{
    #include <fftw3.h>
}

#include "castor/matrix.hpp"


namespace castor{


namespace fftw{

/////////////////////////////////
// SINGLE PRECISION TRANSFORMS //
/////////////////////////////////
matrix<std::complex<float>> fft(matrix<float> &X, int dim=1);
matrix<std::complex<float>> fft(matrix<std::complex<float>> &X, int dim=1);
matrix<std::complex<float>> ifft(matrix<float> &X, int dim=1);
matrix<std::complex<float>> ifft(matrix<std::complex<float>> &X, int dim=1);

// interface to forward/backward 1d FFT depending on sign
matrix<std::complex<float>> xfft(matrix<std::complex<float>> &X, int sign, int dim=1)
{
    auto m = size(X,1); // nb. of lines
    auto n = size(X,2); // nb. of columns
    bool isMatrix = m == 1 || n == 1;
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

// // interface to forward/backward 2d FFT depending on sign
// matrix<std::complex<float>> fftw::xfft2(matrix<std::complex<float>> &X, int sign)
// {
//     auto m = size(X,1); // nb. of lines
//     auto n = size(X,2); // nb. of columns
//     matrix<std::complex<float>> Y = zeros<std::complex<float>>(m,n);
//     return Y;
// }


// END OF NAMESPACE
}

#endif // CASTOR_FFTW_HPP