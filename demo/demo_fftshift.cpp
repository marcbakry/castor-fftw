#include <iostream>

#include "castor/matrix.hpp"
#include "castor/castor_fftwf.hpp"
#include "castor/castor_fftwm.hpp"

using namespace castor;

int main()
{
    std::cout << "+------------------------------------------+" << std::endl;
    std::cout << "| TESTING CASTOR-FFTW: FFTFREQ/(I)FFTSHIFT |" << std::endl;
    std::cout << "+------------------------------------------+" << std::endl;
    //
    std::cout << " - 1D version" << std::endl;
    std::cout << " -- N odd" << std::endl;
    std::size_t N = 7;
    matrix<float> freq = fftw::fftfreq(N,1./N);
    disp(freq,1);
    disp(fftw::fftshift(freq),1);
    disp(norm(fftw::ifftshift(fftw::fftshift(freq)) - freq,"inf"),1);

    std::cout << std::endl;
    std::cout << " -- N even" << std::endl;
    N = 6;
    freq = fftw::fftfreq(N,1./N);
    disp(freq,1);
    disp(fftw::fftshift(freq),1);
    disp(norm(fftw::ifftshift(fftw::fftshift(freq)) - freq,"inf"),1);

    //
    std::cout << std::endl << " - 2D version: with dimension parameter" << std::endl;
    std::size_t M = 7;
    matrix<float> freq2 = zeros(M,N);
    for(std::size_t i=0;i<M;++i)
    {
        for(std::size_t j=0; j<N; ++j) freq2(i,j) = freq(j);
    }
    disp(freq2,1);
    disp(fftw::fftshift(freq2,2),1);
    disp(norm(fftw::ifftshift(fftw::fftshift(freq2,2))-freq2,"inf"),1);

    freq2 = transpose(freq2);
    disp(fftw::fftshift(freq2,1),1);
    disp(norm(fftw::ifftshift(fftw::fftshift(freq2,1))-freq2,"inf"),1);

    //
    std::cout << std::endl << " - 2D version" << std::endl;
    freq2 = matrix<float>({{0,1,2},{3,4,-4},{-3,-2,-1}});
    disp(fftw::fftshift(freq2),1);
    disp(norm(fftw::ifftshift(fftw::fftshift(freq2))-freq2,"inf"),1);
    // the end
    return EXIT_SUCCESS;
}