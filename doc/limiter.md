# Slope Limiters

[TOC]

The slope limiter plays an important role in suppressing spurious numerical oscillations.
They can be written in different equivalent forms and may be equivalent to the corresponding flux limiters.
Most of them are derived strictly for a single nonlinear conservation equation or constant coefficient equations in 1D space.
While for nonlinear conservative equations, slope limiters lack strict theoretical proof.
That might be why sometimes they are not so reliable.

AMRVAC provides many limiters to reconstruct data at cell face.
Most of them are applied for finite volume methods.
Therefore, take care when using them in combination with finite difference schemes, especially for asymmetric ones.
Most TVD limiters in mod_limiter.t can be found in [Yee (1989)](https://ntrs.nasa.gov/citations/19890016281).

In the current version, limiters can be chosen from the following list.

Limiter Type | Limiter | Order | Ghost cells | Reference
---|---|---|---|---
TVD | 'minmod' | 2 | 2 | e.g., Roe (1985), Yee (1989), LeVeque (2002), Toro (2009)
TVD | 'superbee' | 2 | 2 | Sweby (1984), Roe(1985)
TVD | 'vanleer' | 2 | 2 | van Leer (1974)
TVD | 'woodward' | 2 | 2 | van Leer (1977), Woodward et al. (1984)
TVD | 'mcbeta' | 2 | 2 | van Leer (1977) 
TVD | 'albada' | 2 | 2 | van Albada et al. (1974)
TVD | 'koren' | 3 | 2 | Koren (1993)
TVD | 'ppm' | 3 | 3 | Colella and Woodward (1984)
beyond TVD | 'cada' | 2 | 2 | Cada et al. (2009)
beyond TVD | 'cada3' | 3 | 2 | Cada et al. (2009)
beyond TVD | 'schmid1' | 3 | 2 | Schmidtmann et al. (2016)
beyond TVD | 'schmid2' | 3 | 2 | Schmidtmann et al. (2016)
beyond TVD | 'venk' | 2 | 2 | Venkatakrishnan (1995)
ENO-based | 'weno3' | 3 | 2 | Jiang et al. (1996)
ENO-based | 'wenoyc3' | 3 | 2 | Yamaleev et al. (2009), Arandiga et al. (2014)
ENO-based | 'weno5' | 5 | 3 | Jiang et al. (1996), Shu (2009)
ENO-based | 'wenoz5' | 5 | 3 | Borges et al. (2008)
ENO-based | 'wenozp5' | 5 | 3 | Acker et al. (2016)
ENO-based | 'weno5nm' | 5 | 3 | Huang et al. (2018)
ENO-based | 'wenoz5nm' | 5 | 3 | Huang et al. (2018)
ENO-based | 'wenozp5nm' | 5 | 3 | Huang et al. (2018)
ENO-based | 'teno5ad' | 5 | 3 | Peng et al. (2021)
ENO-based | 'weno5cu6' | 6 | 3 | Huang et al. (2018)
ENO-based | 'weno7' | 7 | 4 | Balsara et al. (2000)
ENO-based | 'mpweno7' | 7 | 4 | Balsara et al. (2000)
              
'minmod', 'superbee': Both are classic second-order symmetric TVD limiters. They are found in most textbooks or review papers. as listed in the table above. 'minmod', a.k.a. minbee or mina, might be the most diffusive (non-linear) limiter and thus suitable for testing codes or modules. While 'superbee', or sometimes known as supera, is sharper compared with other second-order limiters. It was designed by Roe (1985), devised by the "discontinuity-sharpning" suggestion from Colella (1984). But it first appear in Sweby (1984) already. The name of A or B (or nickname bee) comes from the A-function or B-function used in the early literatures to quantify the slope.

'vanleer': A classic and widely used second-order symmetric TVD limiter named after Prof. van Leer. It is robust and less diffusive than 'minmod'. Recommended when using second-order limiters.

'woodward': a.k.a. double minmod limiter or monotonized central (MC) limiter, is also a second-order symmetric TVD limiter. Should behave better than 'minmod' and 'superbee' and similar with 'vanleer', also works well in practical simulations.

'mcbeta': A second-order symmetric TVD limiter with a parameter \f$\beta\f$ whose default value is 1.4. Some other usually used values could be 1.5 (see the H-AMR code) or 1.9 (see the KORAL code). In principle, this \f$\beta\f$ could be changed from 1 to 2. If \f$\beta\f$ = 1, it will degenerate into the 'minmod' limiter and if \f$\beta\f$ = 2, it will become the 'woodward' limiter.

'albada': Also a classic second-order symmetric TVD limiter. Similar with the 'vanleer' limiter but smoother. Note that \f$\Delta^2\f$ is set to be 1.d-12 in our code to avoid clipping smooth extrema. However, a better idea to set \f$\Delta^2\f$ is to make it proportional to \f$\Delta x^3\f$, as we can see in the settings of the 'venkatakrishnan' limiter or 'wenoyc3' limiter.

'koren': A third-order asymmetric TVD limiter. Uses less resources than 'cada3' but is more diffusive as a third-order limiter. Maybe as robust as 'wenoyc3'.

'ppm': A three-point centered stencil third-order TVD limiter, the acronym stands for the Piecewise Parabolic Method (Colella and Woodward 1984). See also the references in mod_ppm.t for the implementation of this method.

'cada', 'cada3': Second-order and third-order asymmetric limiter, a.k.a. LIMO or LIMO3. They are designed to simplify the complex TVD limiter created by Artebrant et al. (2005).  However, the simplication makes it not strictly TVD any more. As a third-order limiter, 'cada3' performs better than 'koren' limiter. Thus, although it has some drawbacks like being asymmetric or relying on some artificial parameters, this limiter is still recommended when choosing third-order limiters.

'schmid1', 'schmid2': An improved version 'cada3' limiter. It changed some parameters of the 'cada3' limiter so that it is now a symmetric limiter. Besides, it changes the criteria to switch on/off the TVD function in 'cada3' from a quite artificial approach to a more reliable method. According to Schmidtmann (2016), it should be related to "the max value of the second derivative of all the variables", which is calculated automatically in 'schmid1' settings. But considering the time cost of communication, 'schmid2' is recommended. The users are expected to give "the max value of the second derivative of all the variables at t=0" manually. See examples of the Shu-Osher test in the hd test folder.

'venk': The 'venkatakrishnan' limiter is popular in unstructured meshes. It is a little bit similar with the 'van albada' limiter. This limiter is claimed to keep a better balance between accuracy and oscillations, and thus can provide some more details than other second-order limiters. The parameter K, as mentioned in their paper, is set to be 0.3. Other values like 1 or 10 are also widely used in literature.

'mp5': MP5 is a high order (fifth-order) five-point Monotonicity Preserving limiter, which preserves the so-called MP property. MP5 can provide more details than other limiters provided in AMRVAC, but since this MP property cannot prevent small oscillations, even with small CFL number (0.4 or 0.2), it may crash in many cases. Note that the parameter \f$\alpha\f$ is chosen as 4, which is recommended Suresh et al. (1997). Therefore, in principle, CFL number should be less than 1 / (1 + \f$\alpha\f$) = 0.2. But in practice, 0.4 still yields nonoscillatory results.

'weno3', 'weno5', 'weno7': The classic WENO or Weighted Essentially Non(-)oscillatory Scheme, a popular variation of ENO scheme, which means that we allow some kind of oscillations that is not allowed in the TVD scheme, but these oscillations could not be too large. This idea is similar with the so-called TVB scheme, where B stands for bounded. Check Liu et al. (1994)  for the original idea of WENO scheme. The so-called classic WENO or WENO-JS scheme is from Jiang et al. (1996), who modified it into a more practical way. High-order WENO schemes are especially suitable for problems containing both shocks and a large number of complex smooth structures. 'weno3', 'weno5', 'weno7' stand for third-order, fifth-order and seventh-order WENO-JS scheme, respectively. As mentioned in Shu (2009), the coefficients named 'reconstruction' are used here. And in 'weno5', we also provide another set of coefficients named 'interpolation', as applied in the ECHO code (Del Zanna et al. 2007). Check mod_weno.t if it is necessary to use this set of coefficients. The parameter \f$\Delta^2\f$, which has a similar usage with the van Albada and Venkatakrishnan limiter, is set to be \f$ 10^{-18}\f$ by default. Values as low as \f$ 10^{-40}\f$ is also possible. Change it if your grid size is very small. For the performance, although the WENO scheme will allow some oscillations, it seems to be more diffusive than other third-order or fifth-order schemes, and thus robust in many cases. Therefore, as a compromise between robutstness and accuracy (and maybe also efficiency), 'weno3' is note so recommended compared with other third-order limiters while 'weno5' should be a good choice when using fifth-order limiters. Notice should be taken that 'weno7' perfers to using characteristic variables instead of primitive or conservative variables (actually, this is true for all the limiters, but hige-order limiters are more sensitive). Thus, in the current version of AMRVAC, 'weno7' is not so robust and not recommended.

'wenoz5': Fifth-order WENO-Z scheme. A popular variation of the classic WENO scheme which is innovated by WENO-M but more efficient than WENO-M. It claims to perform better at critical points than WENO-JS, but is not as robust as WENO-JS in some extreme cases. Although the power index is set to be 1 in Borges et al. (2008), we still set it to be 2 so that it works in a similar way with 'weno5' and 'wenozp5'. Change it in mod_weno.t if necessary. Apart from the widely used value 2, values from 1 to 6 could also be found in previous papers.

'wenoyc3': It could be considered as a thrid-order version of 'wenoz5', see Yamaleev et al. (2009). But the version installed in the code is actually from Arandiga et al. (2014).

'wenozp5': Fifth-order WENO-Z+ scheme. An improved version of WENO-Z. More precision in the sacrifice of robustness.

'weno5nm', 'wenoz5nm', 'wenozp5nm': In principle, all the limiters mentioned above are designed for uniform grids. While in stretched grids, they are not strictly correct so that it cannot be as accurate as what they claimed to be. While these -NM limiters are specially designed to deal with stretched Cartesian grids. The precision will improve in the sacrifice of a little bit efficiency.

'mpweno7': Seventh-order MP-WENO scheme. It is basically a combination of MP5 and WENO scheme to preserve the MP property in the results of especially high-order WENO schemes, e.g. seventh-order of even higher. The advantage might be that it is more stable than 'weno7'. However, it shows two drawbacks in the tests, one is this scheme is time-comsuming, and another is the result is more diffusive, which is a straighforward result from the design of this scheme. If 'weno5' or 'weno7' could already be used very well in your case, 'mpweno7' is not recommeneded.

'weno5cu6': CU means central-upwind here, so actually this limiter should be called WENO6. It is based on two ideas. First, in recent desigments for WENO-type schemes, instead of the traditional way of making combination of low-order linear reconstructions, the combination of high-order and low-order linear reconstructions is preferred. Second, in the smooth region, central scheme can always give less dissipative solution than the traditional upwind schemes. Huang et al. (2018) gives sucn an example to construct sixth-order scheme based on the idea above. In our tests, however, the robustness of this 'weno5cu6' is sort of disappointing. 

'teno5ad': Recently, T(argeted) ENO or TENO schemes are popular in ENO-based scheme designments. WENO-type schemes will always be non-linear so that diffusive even in the smooth region. While TENO-type schemes aim to be totally linear in the regions they recognize as smooth region. Here we select a recent one in Peng et al. (2021), where AD here means adaptive dissipation. However, in the current stage, TENO-type schemes are still sort of in experiment, tends to crash in many practical applications.

References
1. Acker, F., 2016, An Improved WENO-Z Scheme.
2. Arandiga, F., 2014, Weights Design for Maximal Order WENO Schemes.
3. Artebrant, R., 2005, Conservative Logarithmic Reconstructions and Finite Volume Methods.
4. Balsara, D., 2000, Monotonicity Preserving Weighted Essentially Non-Oscillatory Scheme with Increasingly High Order of Accuracy.
5. Borges, R., 2008, An Improved Weighted Essentially Non-Oscillatory Scheme for Hyperbolic Conservation Laws.
6. Cada, M., 2009, Compact Third-Order Limiter Functions for Finite Volume Methods.
7. Colella, P., 1984, The Piecewise Parabolic Method (PPM) for Gas-Dynamical Simulations.
8. Del Zanna, L., 2007, ECHO: a Eulerian Conservative High-Order Scheme for General Relativistic Magnetohydrodynamics and Magnetodynamics.
9. Huang, C., 2018, A new adaptively central-upwind sixth-order WENO scheme.
10. Huang, W., 2018, A Simple Algorithm to Improve the Performance of the WENO Scheme on Non-uniform Grids.
11. Jiang, G., 1996, Efficient Implementation of Weighted ENO Schemes.
12. Koren, B., 1993, A Robust Upwind Discretisation Method for Advection, Diffusion and Source Terms.
13. LeVeque, R., 2002, Finite Volume Methods for Hyperbolic Problems.
14. Liu, X., 1994, Weighted Essentially Non-Oscillatory Schemes.
15. Peng, J., 2021, An efficient targeted ENO scheme with local adaptive dissipation for compressible flow simulation.
16. Roe, P., 1985, Some Contributions to the Modelling of Discontinuous Flows.
17. Schmidtmann, B., 2016, Relations Between WENO3 and Third-Order Limiting in Finite Volume Methods.
18. Shu, C., 2009, High Order Weighted Essentially Nonoscillatory Schemes for Convection Dominated Problems.
19. Suresh, A., 1997, Accurate Monotonicity-Preserving Schemes with Runge–Kutta Time Stepping.
20. Sweby, P, 1984, High Resolution Schemes Using Flux Limiters for Hyperbolic Conservation Laws.
21. Toro, E, 2009, Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction.
22. van Albada, G., 1982, A Comparative Study of Computional Methods in Cosmic Gas Dynamic.
23. van Leer, B., 1974. Towards the Ultimate Conservative Difference Scheme, II: Monotonicity and Conservation Combined in a Second-Order Scheme.
24. van Leer, B., 1977, Towards the Ultimate Conservative Difference Scheme. III. Upstream-Centered Finite-Difference Schemes for Ideal Compressible Flow
25. Venkatakrishnan, V., 1993, On the Accuracy of Limiters and Convergence to Steady State Solutions.
26. Woodward, P., 1984, The Numerical Simulation of Two-Dimensional Fluid Flow with Strong Shock.
27. Yamaleev, N., 2009, A systematic methodology for constructing high-order energy stable WENO schemes
28. Yee, H., 1989, A Class of High-Resolution Explicit and Implicit Shock-Capturing Methods.
