This section describes the API of **castor-fftw**.

Main API
++++++++

.. _label-fft:

fft
---
.. doxygenfunction:: fft(matrix<float> &X, int dim = 1)
    :project: castorfftw
.. doxygenfunction:: fft(matrix<std::complex<float>> &X, int dim = 1)
    :project: castorfftw
.. doxygenfunction:: fft(matrix<double> &X, int dim = 1)
    :project: castorfftw
.. doxygenfunction:: fft(matrix<std::complex<double>> &X, int dim = 1)
    :project: castorfftw

See also :ref:`label-ifft`.


.. _label-ifft:

ifft
----
.. doxygenfunction:: ifft(matrix<float> &X, int dim = 1)
    :project: castorfftw
.. doxygenfunction:: ifft(matrix<std::complex<float>> &X, int dim = 1)
    :project: castorfftw
.. doxygenfunction:: ifft(matrix<double> &X, int dim = 1)
    :project: castorfftw
.. doxygenfunction:: ifft(matrix<std::complex<double>> &X, int dim = 1)
    :project: castorfftw

See also :ref:`label-fft`.


.. _label-fft2:

fft2
----
.. doxygenfunction:: fft2(std::size_t m, std::size_t n, matrix<float> &X)
    :project: castorfftw
.. doxygenfunction:: fft2(std::size_t m, std::size_t n, matrix<std::complex<float>> &X)
    :project: castorfftw
.. doxygenfunction:: fft2(matrix<float> &X)
    :project: castorfftw
.. doxygenfunction:: fft2(matrix<std::complex<float>> &X)
    :project: castorfftw
.. doxygenfunction:: fft2(std::size_t m, std::size_t n, matrix<double> &X)
    :project: castorfftw
.. doxygenfunction:: fft2(std::size_t m, std::size_t n, matrix<std::complex<double>> &X)
    :project: castorfftw
.. doxygenfunction:: fft2(matrix<double> &X)
    :project: castorfftw
.. doxygenfunction:: fft2(matrix<std::complex<double>> &X)
    :project: castorfftw

See also :ref:`label-ifft2`.

.. _label-ifft2:

ifft2
-----
.. doxygenfunction:: ifft2(std::size_t m, std::size_t n, matrix<float> &X)
    :project: castorfftw
.. doxygenfunction:: ifft2(std::size_t m, std::size_t n, matrix<std::complex<float>> &X)
    :project: castorfftw
.. doxygenfunction:: ifft2(matrix<float> &X)
    :project: castorfftw
.. doxygenfunction:: ifft2(matrix<std::complex<float>> &X)
    :project: castorfftw
.. doxygenfunction:: ifft2(std::size_t m, std::size_t n, matrix<double> &X)
    :project: castorfftw
.. doxygenfunction:: ifft2(std::size_t m, std::size_t n, matrix<std::complex<double>> &X)
    :project: castorfftw
.. doxygenfunction:: ifft2(matrix<double> &X)
    :project: castorfftw
.. doxygenfunction:: ifft2(matrix<std::complex<double>> &X)
    :project: castorfftw

See also :ref:`label-fft2`.

Helper functions
++++++++++++++++

fftfreq
-------
.. doxygenfunction:: fftfreq(std::size_t n, T d)
    :project: castorfftw


.. _label-fftshift:

fftshift
--------
.. doxygenfunction:: fftshift(matrix<T> const &A, int dim=0)
    :project: castorfftw

See :ref:`label-ifftshift`, :ref:`label-fft`, :ref:`label-fft2`.


.. _label-ifftshift:

ifftshift
---------
.. doxygenfunction:: ifftshift(matrix<T> const &A, int dim=0)
    :project: castorfftw

See :ref:`label-fftshift`, :ref:`label-ifft`, :ref:`label-ifft2`.