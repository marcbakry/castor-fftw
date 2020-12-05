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
|                    within the 'castor' framework for      |
|                    double precision accuracy.             |
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

//////////////////////////////////////////
// DECLARATIONS OF THE PUBLIC INTERFACE //
//////////////////////////////////////////
matrix<std::complex<double>> fft(matrix<double> &X, int dim=1);
matrix<std::complex<double>> fft(matrix<std::complex<double>> &X, int dim=1);
matrix<std::complex<double>> ifft(matrix<double> &X, int dim=1);
matrix<std::complex<double>> ifft(matrix<std::complex<double>> &X, int dim=1);
matrix<std::complex<double>> fft2(matrix<double> &X);
matrix<std::complex<double>> fft2(matrix<std::complex<double>> &X);
matrix<std::complex<double>> fft2(std::size_t m, std::size_t n, matrix<double> &X);
matrix<std::complex<double>> fft2(std::size_t m, std::size_t n, matrix<std::complex<double>> &X);
matrix<std::complex<double>> ifft2(matrix<double> &X);
matrix<std::complex<double>> ifft2(matrix<std::complex<double>> &X);
matrix<std::complex<double>> ifft2(std::size_t m, std::size_t n, matrix<double> &X);
matrix<std::complex<double>> ifft2(std::size_t m, std::size_t n, matrix<std::complex<double>> &X);

matrix<double> fftfreq(std::size_t n, double d = 1.);


//////////////////////////////
// NON-DOCUMENTED INTERFACE //
//////////////////////////////

// interface to forward/backward 1d FFT depending on sign
matrix<std::complex<double>> xfft(matrix<std::complex<double>> &X, int sign, int dim=1)
{
    auto m = size(X,1); // nb. of lines
    auto n = size(X,2); // nb. of columns
    bool isMatrix = (m != 1 && n != 1);
    matrix<std::complex<double>> Y = zeros<std::complex<double>>(m,n);
    fftw_complex *in = reinterpret_cast<fftw_complex*>(&X(0));
    fftw_complex *out = reinterpret_cast<fftw_complex*>(&Y(0));
    if(!isMatrix) // vector case
    {
        fftw_plan plan = fftw_plan_dft_1d(std::max(m,n), in, out, sign, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
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
        fftw_plan plan = fftw_plan_many_dft(rank,&N[0],howmany,in,nullptr,istride,idist,out,nullptr,istride,idist,sign,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    }
    return Y;
}


// interface to forward/backward 2d FFT depending on sign
matrix<std::complex<double>> xfft2(matrix<std::complex<double>> &X, int sign)
{
    auto m = size(X,1); // nb. of lines
    auto n = size(X,2); // nb. of columns
    matrix<std::complex<double>> Y = zeros<std::complex<double>>(m,n);
    fftw_complex *in  = reinterpret_cast<fftw_complex*>(&X(0));
    fftw_complex *out = reinterpret_cast<fftw_complex*>(&Y(0));
    //
    fftw_plan plan  = fftw_plan_dft_2d(m,n,in,out,sign,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return Y;
}

// interface to forward/backward 2d FFT depending on sign
matrix<std::complex<double>> xfft2(std::size_t m, std::size_t n, matrix<std::complex<double>> &X, int sign)
{
    // check consistency
    if(m*n != size(X,1)*size(X,2)) error(__FILE__, __LINE__, __FUNCTION__,"Size of input given as parameter is not consistent with the number of elements in the input matrix.");
    matrix<std::complex<double>> Y = zeros<std::complex<double>>(size(X,1),size(X,2));
    fftw_complex *in  = reinterpret_cast<fftw_complex*>(&X(0));
    fftw_complex *out = reinterpret_cast<fftw_complex*>(&Y(0));
    //
    fftw_plan plan  = fftw_plan_dft_2d(static_cast<int>(m),static_cast<int>(n),in,out,sign,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    return Y;
}

}

matrix<std::complex<double>> fftw::fft(matrix<std::complex<double>> &X, int dim)
{
    return fftw::xfft(X,FFTW_FORWARD,dim);
}
matrix<std::complex<double>> fftw::fft(matrix<double> &X, int dim)
{
    matrix<std::complex<double>> XX = cast<std::complex<double>>(X);
    return fftw::xfft(XX,FFTW_FORWARD,dim);
}
matrix<std::complex<double>> fftw::ifft(matrix<std::complex<double>> &X, int dim)
{
    return fftw::xfft(X,FFTW_BACKWARD,dim);
}
matrix<std::complex<double>> fftw::ifft(matrix<double> &X, int dim)
{
    matrix<std::complex<double>> XX = cast<std::complex<double>>(X);
    return fftw::xfft(XX,FFTW_BACKWARD,dim);
}
matrix<std::complex<double>> fftw::fft2(matrix<std::complex<double>> &X)
{
    return fftw::xfft2(X,FFTW_FORWARD);
}
matrix<std::complex<double>> fftw::fft2(matrix<double> &X)
{
    matrix<std::complex<double>> XX = cast<std::complex<double>>(X);
    return fftw::xfft2(XX,FFTW_FORWARD);
}
matrix<std::complex<double>> fftw::fft2(std::size_t m, std::size_t n, matrix<std::complex<double>> &X)
{
    return fftw::xfft2(m,n,X,FFTW_FORWARD);
}
matrix<std::complex<double>> fftw::fft2(std::size_t m, std::size_t n, matrix<double> &X)
{
    matrix<std::complex<double>> XX = cast<std::complex<double>>(X);
    return fftw::xfft2(m,n,XX,FFTW_FORWARD);
}
matrix<std::complex<double>> fftw::ifft2(matrix<std::complex<double>> &X)
{
    return fftw::xfft2(X,FFTW_BACKWARD);
}
matrix<std::complex<double>> fftw::ifft2(matrix<double> &X)
{
    matrix<std::complex<double>> XX = cast<std::complex<double>>(X);
    return fftw::xfft2(XX,FFTW_BACKWARD);
}
matrix<std::complex<double>> fftw::ifft2(std::size_t m, std::size_t n, matrix<std::complex<double>> &X)
{
    return fftw::xfft2(m,n,X,FFTW_BACKWARD);
}
matrix<std::complex<double>> fftw::ifft2(std::size_t m, std::size_t n, matrix<double> &X)
{
    matrix<std::complex<double>> XX = cast<std::complex<double>>(X);
    return fftw::xfft2(m,n,XX,FFTW_BACKWARD);
}


matrix<double> fftw::fftfreq(std::size_t n, double d)
{
    matrix<double> freq;
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

// END OF NAMESPACE
}

#endif // CASTOR_FFTW_HPP