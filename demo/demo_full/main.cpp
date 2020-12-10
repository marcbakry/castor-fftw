#include "castor/matrix.hpp"
#include "castor/graphics.hpp"
#include "castor/castor_fftw.hpp"
#include "castor/castor_fftwm.hpp"

using namespace castor;

int main()
{
    //
    std::size_t N    = 1000;
    double      tmin = 0.;
    double      tmax = 1.;
    double      dt   = (tmax-tmin)/N;
    // initialize figures
    figure fig;
    figure fighat;
    //
    matrix<double> t = zeros(1,N);
    for(auto i=0; i<N; ++i) t(i) = i*dt;

    matrix<double> f = 2.*sin(2*M_PI*3.*t) + 10.*sin(2*M_PI*10.*t);
    // add 'f' to plot
    plot(fig,t,f,{"r-+"},{"'f'"});

    matrix<std::complex<double>> fhat = fftw::fftshift(fftw::fft(f))/(double)N;
    matrix<double>               freq = fftw::fftshift(fftw::fftfreq(N,dt));
    // add fhat to plot
    plot(fighat,freq,abs(fhat),{"r-+"},{"F(f)"});

    // now we filter the data
    auto fhat_filtered = fhat;
    for(auto i=0; i<N; ++i)
    {
        if(std::abs(freq(i)) > 5.) fhat_filtered(i) = 0;
    }
    plot(fighat,freq,abs(fhat_filtered),{"b-"},{"F(f) filtered"});

    // and we compute the backward transform
    fhat_filtered = fftw::ifftshift(fhat_filtered);
    auto f_filtered = fftw::ifft(fhat_filtered);

    // plot filtered signal
    plot(fig,t,real(f_filtered),{"b-"},{"'f' filtered"});

    // draw all figures
    drawnow(fig);

    // the end
    return EXIT_SUCCESS;
}