# Verification and benchmark of time domain code module of Clara2.0
> Author: Hao PENG, email: penghao1028311@gmail.com Data：7/6/2021
> 
> This module is based on [Clara2.0](https://github.com/ComputationalRadiationPhysics/clara2) (authorized by Richard Pausch, email: r.pausch@hzdr.de) and his diploma [thesis](http://www.hzdr.de/db/Cms?pOid=38997).

## Introduction to the module

This module is basically devoted to calculating the electrical field (including near field and radiation field) of an electron bunch using their trajectories as input, based on the Lienard-Wichert Potentials: 

$\vec{E}(\vec{r}, t)=\frac{q}{4\pi\epsilon_0}\left[ \frac{\vec{n}-\vec{\beta}}{R^2\gamma^2(1-\vec{\beta}\cdot\vec{n})^3} + \frac{\vec{n}\times[(\vec{n}-\vec{\beta})\times\dot{\vec{\beta}}]}{cR(1-\vec{\beta}\cdot\vec{n})^3}\right]$,

at advanced time: $t_{adv}=t+\frac{|\vec{R}|}{c}$. Where $\vec{R}(t)$ is the vector from the particle to the observation point at $t$: $\vec{R}(t)=\vec{r}_{obs}-\vec{r}_e(t)$, and $\vec{n}$ is the unit vector: $\vec{n}=\frac{\vec{R}}{|\vec{R}|}$.

In the *setting.cpp* of the source code, users can turn on this module by setting *USE_uop=true*, and define the position and and the observation time domain of the observation plane by following parameters: *X_uop*, *Y_max*, *Y_min*, *Z_min*, *Z_max*, *Ny*, *Nz*, *t_start*, *dt_obs*, *Nt_obs*. The scheme is as shown in  [scheme](figures/scheme.jpg). The time domain is from $t_{start}$
to $t_{start}+Nt_{obs}*dt_{obs}$. Other settings are the same as the original module of Clara2.0 to calculate the electron radiation spectrum based on electron trajectories, please see [how to set up a first Clara2 simulation](https://github.com/ComputationalRadiationPhysics/clara2/issues/96) and [Clara2 wiki](https://github.com/ComputationalRadiationPhysics/clara2/wiki/Users).

Beside setting up the input parameters in *setting.cpp*, users need to put the electron traces in the *src* directories in the form of *trace_%04d.txt*. Once the code is successfully run and ended, files in the form of *my_eField_trace_%06d.dat* will be dumped as output, in which electric fields are stored with the nested order: directions -> timesteps -> coordinations(separated with newlines, tabs and spaces, respectively), as shown in [figure](figures/structureOfOutputFile.jpg). Then by running *./process_data_uop*, the electric field emitted by different electrons(calculated based on different traces) will be added coherently. At last files *my_efield_all_%03d.dat* will be dumped, each with the electric fields in one *y* direction. For example, *my_efield_all_000.dat* stores the electric fields with *y=y_min*. (If users want to input a little bit more traces, say 1000 traces, running *./process_data_uop* will be super slow and the dumped files will cost huge memory. I am now modifying the code, so that one of the processor will receive the computed electric fields from other processors and add them together. Then no more separate electric field files will be dumped, only one last file with the coherent added electric fields will be dumped.)

The module is benchmarked against three different simulations, which are basically the same as those in Section 5.1 of Richard Pausch's [thesis](https://zenodo.org/record/843510#.YMBfKb7iuUk), except in his thesis, the observation plane is in the same plane of the electron movement, while for the benchmarks of the first two simulations here the observation plane is not. 

## Benchmark 1: static electron

For the first benchmark simulation, a static electron is put at the origin, while the observation plane is put at $x_{uop}=-13.6~{\rm \mu m}$, $y_{max}=z_{max}=2.1~{\rm \mu m}$, $y_{min}=z_{min}=-2~{\rm \mu m}$, $N_y=N_z=41$, $t_{start}=4.0\times10^{-14} s$, $dt_{obs}=1.0\times10^{-15} s$ and $Nt_{obs}=20$.The trace file is [staticElectron_trace_0000.txt](staticElectron_trace_0000.txt) (note if one wants to use this as input file, one needs change the name to *trace_0000.txt*). Theoretically, only near field (i.e. the Coulomb field) is emitted in this case, which can be calculated as:

$\vec{E}=\frac{q}{4\pi\epsilon_0}\frac{\vec{e}_r}{r^2}$.

The difference between the theoretical value and the simulated electric field is shown in [figure](./figures/deltaAbsEInDetectorPlane.pdf). Good agreement is found and the relative error is on the order of $10^{-12}$.

## Benchmark 2: electron at constant speed

For the second benchmark simulation, the trace of an electron moving at a constant speed($\beta=0.9$ towards $+x$) is used as input file. The trace file is [constantSpeed_trace_0000.txt](constantSpeed_trace_0000.txt). The electron passes the origin at $t=0$.  The observation plane is put at $x_{uop}=-13.6~{\rm \mu m}$, $y_{max}=z_{max}=2.1~{\rm \mu m}$, $y_{min}=z_{min}=-2~{\rm \mu m}$, $N_y=N_z=41$, $t_{start}=-2.0\times10^{-14} s$, $dt_{obs}=1.0\times10^{-15} s$ and $Nt_{obs}=80$. In this case, the electron also has no accelerations, only near field is generated. Theoretically, the electrical field at $t=0$ should be:

$\vec{E}=\frac{q}{4\pi\epsilon_0\sqrt{1-\beta^2}}\frac{\vec{r}}{\left[\frac{x^2}{1-\beta^2}+y^2+z^2\right]^{3/2}}$.

The difference between the theoretical value and the simulated electric field is shown in [figure](./figures/deltaAbsEInDetectorPlane_beta_0d9.pdf). Still, good agreement is found and the relative error is on the order of $10^{-4}$ to $10^{-3}$.

## Benchmark 3: synchrotron radiation

For the last benchmark simulation, the electron is moving in a circle, similar to the electron motion in a synchrotron. It moves clockwise with a radius of $\rho=1 m$ and a velocity $\beta=0.9$, in the $y-z$ plane with $x=0$. The trace file is [synchrotron_trace_0000.txt](synchrotron_trace_0000.txt). The observation plane is in the same plane of the electron motion, i.e. $x_{uop}=0~{\rm \mu m}$. Other input parameters are $y_{max}=z_{max}=30.01~{\rm m}$, $y_{min}=z_{min}=-30~{\rm m}$, $N_y=N_z=601$, $t_{start}=-1.5\times10^{-7} s$, $dt_{obs}=\times10^{-10} s$ and $Nt_{obs}=3000$. At $t=0$, the electron moves to $(y, z)=(0, 1m)$, and the electrical field (including near field and radiation field) at $t=0$ is shown in [figure](./figures/eFieldInDetectorPlane_synchrotron.pdf), which is quite similar to Fig 5.3 in R. Pausch's [thesis](http://www.hzdr.de/db/Cms?pOid=38997). 

Theoretically, because the electron trajectory is regular, so the "emitting time"(or retarded time) of electric field at any positions at a particular time can be determined numerically as explained in [Ref](https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.16.010701). To check the coherence between the theory and the simulation following R. Pausch, the one dimensional electric field at $z=0$, from $y=0$ to $y=30~m$ at $t=0$ is simulated. The input parameters are $x_{uop}=0~{\rm \mu m}$, $y_{max}=30.000625~{\rm m}$, $y_{min}=0~{\rm m}$, $N_y=4801$, $y_{max}=0~{\rm m}$, $y_{min}=0~{\rm m}$, $N_z=1$, $t_{start}=-1.5\times10^{-7} s$, $dt_{obs}=\times10^{-10} s$ and $Nt_{obs}=3000$. The near field(or the Coulomb field, the velocity field) and the radiation field are dumped and compared to the theoretical value separately, as shown in [figure](./figures/Evel.pdf) and [figure](./figures/Erad.pdf), respectively. They agree very well.
