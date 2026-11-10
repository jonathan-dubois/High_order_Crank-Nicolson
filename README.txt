
# High-order Crank-Nicolson method via diagonal Padé factorization

MATLAB code accompanying
> L. Freschi, C. Lévêque, J. Caillat, R. Taïeb and J. Dubois,
> High-order Crank-Nicolson method via diagonal Padé factorization,
> submitted to Computer Physics Communications (2026).

The code solves the time-dependent Schrödinger equation with a propagator of arbitrary order 2q in the time step dt. 
One time step costs q tridiagonal matrix-vector products and q Thomas solves, i.e. O(qN) operations. 
For q=1 the scheme reduces to the standard Crank-Nicolson method.

For Hamiltonians that are sums of tridiagonal and diagonal operators, such as multi-dimensional grids, the propagator is combined with the symmetric high-order splitting schemes of McLachlan (Verlet, BM4, BM6). 
The convergence order is then set by the splitting, for every q.

## Contents

| Script          | Reproduces    | Description                                                             
|-----------------|---------------|--------------------------------------------------------------------
| 'Quantum1D_N.m' | Figs. 1 and 2 | 1D free particle; Padé orders q=1,...,5; error, norm and timing vs. dt 
| 'Quantum2D_N.m' | Fig. 3        | 2D particle in a non-separable Gaussian potential; Verlet, BM4 or BM6 splitting combined with q=1,...,5 

Each script is self-contained.
It includes the tabulated coefficients for q=1,...,9 (Table 1 of the paper) and the splitting coefficients (Table 2).

## Requirements

- MATLAB R2023b or later, for the built-in 'expmv', which computes the reference solution.
- No toolbox is required.

## Usage

Run a script from the MATLAB command window:

>> matlab
Quantum1D_N                 % Figs. 1-2
Quantum2D_N                 % Fig. 3, one splitting scheme per run

For Fig. 3, set the splitting scheme at the top of 'Quantum2D_N.m' and run the script once per scheme:

>> matlab
P.integrator = 'BM4';       % 'Verlet', 'BM4' or 'BM6'

All physical and numerical parameters are defined in the 'Parameters' block at the top of each script. 
The defaults are those used in the paper.

| Parameter              | 1D                             | 2D                                   
|------------------------|--------------------------------|--------------------------------------
| Grid                   | 512 points, dx=0.1             | 512x512, dx=dy=0.1
| Initial wavepacket     | x0=-5, sigma=2, k0=2           | (x0,y0)=(-5,-5), sigma=2, k_0=2
| Potential              | none (free particle)           | V=exp(-0.5x^2-0.5y^2+0.3xy)
| Final time             | T=4                            | T=4                              
| Padé orders            | q=1,...,5                      | q=1,...,5                      
| Time steps             | dt=4/N_t, N_t = 2^0,...,2^15   | dt=4/N_t, N_t=2^0,...,2^15 

All quantities are in atomic units.

## Output

Each script plots, as a function of dt:

- the L^2 error at time T with respect to the 'expmv' reference, and
- the norm deviation |1 - ||psi(T)||^2|.

The results are saved to a '.mat' file ('data2D_<scheme>.mat' for the 2D script) containing the parameters, errors, norm deviations and computation times.

## Runtime

The 1D script runs in a few minutes. A full 2D run down to dt=4/2^15, with all five Padé orders, takes several hours in MATLAB, especially for BM6. 
To test the code quickly, reduce the range of time steps, e.g. 'P.Nt = 2^(0:8);'.

## Notes

- The accuracy of the 'expmv' reference is a few 10^{-12} for these parameters. The errors therefore saturate at this level for small dt (see Sec. 2.4 of the paper).
- In 2D, the wavefunction is stored as an Ny x Nx array. The full Hamiltonian used for the reference is assembled consistently with MATLAB's column-major vectorization, so non-square grids, dx!=dy and asymmetric potentials are all handled correctly.
- The tridiagonal solves are vectorized over grid lines. Each grid line is an independent problem, which makes the method straightforward to parallelize in compiled languages.

## Citation

If you use this code, please cite:

>> bibtex
@article{Freschi2026,
  author  = {Freschi, Lou and L{\'e}v{\^e}que, Camille and Caillat, J{\'e}r{\'e}mie and Ta{\"i}eb, Richard and Dubois, Jonathan},
  title   = {High-order Crank--Nicolson method via diagonal Pad{\'e} factorization},
  journal = {Computer Physics Communications},
  year    = {2026},
  note    = {submitted}
}

## License

Released under the MIT License (see 'LICENSE').

## Contact

Jonathan Dubois, Sorbonne Université, CNRS, Laboratoire de Chimie Physique - Matière et Rayonnement (LCPMR), Paris, France
jonathan.dubois@sorbonne-universite.fr
