Clara 2
=======

Clara2 - a parallel classical radiation calculator based on Liénard-Wiechert potentials


Introduction
------------

Clara stands for CLAssical RAdiation and is a redevelopment from scratch of [Clara1](https://github.com/ComputationalRadiationPhysics/clara1).
It has been developed as part of a [diploma thesis](http://www.hzdr.de/db/Cms?pOid=38997) in 2012.
In contrast to [Clara1](https://github.com/ComputationalRadiationPhysics/clara1) it is parallelized using MPI and OpenMP to run efficiently on large CPU clusters. 


Referencing
-----------

Clara2 is a *scientific project*. If you **present and/or publish** scientific
results that used Clara2, you should set this as a **reference**.

Our according **up-to-date publication** at **the time of your publication**
should be inquired from:
- [REFERENCE.md](REFERENCE.md)



Software License
----------------

*Clara2* is licensed under the **GPLv3+**. You can use any of our *libraries* with
**GPLv3+ or LGPLv3+** (they are *dual-licensed*).
Please refer to our [LICENSE](LICENSE)


Dependency
----------

*Clara2* uses the FFTW library when used with the (faster) fft detector.
If you install *Clara2*, you need to install FFTW as well. You can finde 
compiled code at: http://www.fftw.org/ or the source code at 
https://github.com/FFTW/fftw3. FFTW3 is under GNU General Public License v2.0.
