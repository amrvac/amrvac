# Slope Limiters

The slope limiter plays an importrant role in suppressing the spurious numerical oscillations.
They can be written in different equivalent forms and could be equivalent to the correspoding flux limiters.
Most of them are derived strictly for a single nonlinear conservation equation or constant coefficient equations in 1D space.
While for nonlinear conservative equations, slope limiters lack strict theoretical proof.
That might be why sometimes they are not so reliable.

AMRVAC provide tens of limiters to reconstruct data at cell face.
Most of them are applied for finite volume methods.
Therefore, take care when using them in combination with finite difference schemes, especially for asymmetric ones.
And, be careful when choosing limiters when you are using stretching grids because some of the limiters only suitable for uniform grids.
Most limiters in mod_limiter.t could be found in Reference [1].

In the current version, limiters could be chosen from the following 19 presets.

'minmod', 'superbee': MINMOD and SUPERBEE are both classic second-order symmetric TVD limiters. They could be found in most text books or review papers, for example, Reference [1-3]. MINMOD, a.k.a. MINBEE or MINA, might be the most diffusive limiter and thus suitable for testing codes or modules. While SUPERBEE, or sometimes known as SUPERA, is a little bit sharp compared with other second-order limiters. Another similar limiter in the same family is ULTRABEE, which is not implemneted in our code.

'woodward': WOODWARD, a.k.a. double MINMOD limiter or MONOTONIZED CENTRAL(MC) limiter, is also a second-order symmetric TVD limiter which could be found in Reference [4-5]. Should play a better balance than MINMOD and SUPERBEE.

'mcbeta': MC_beta limiter is a second-order symmetric TVD limiter with a default beta = 1.4. Some other usually used values could be 1.5 (in the H-AMR code) or 1.9 (in the KORAL code). In principle, this beta could be changed from 1 to 2. If beta = 1, it will reduce to the MINMOD limiter and if beta = 2, it will reduce to the WOODWARD limiter. It is also proposed in Reference [4].

'vanleer': A classic and widely used second-order symmetric TVD limiter named after Prof. van Leer. Its performance should be comparable with the WOODWARD limiter. First propsed in Reference [6].

'albada': Also a classic second-order symmetric TVD limiter. Similar with the VANLEER limiter but smoother. See Reference [7] for details if you have access to this book. Note that delta^2 is set to be 1.d-12 in our code to avoid clipping smooth extrema. You can lower this value to, for example, 1.d-42 for current double precision limit to reduce the diffusion. However, a better idea to set delta^2 is to make it proportional to (delta x)^3, as we can see in the settings of VENKATAKRISHNAN limiter.

'koren': A third-order asymmetric TVD limiter. See Reference [8] for details. Cost less resource than CADA3 but more diffusive as a third-order limiter.

'cada', 'cada3': Second-order and third-order asymmetric limiter proposed in Reference [9], a.k.a. LIMO or LIMO3. They are designed to simplify the complex TVD limiter created by Artebrant and Schroll. However, itself is not strictly TVD any more. Related parameters are chosen as -0.5, 2.0, 1.6 as recommended in their paper. As a third-order limiter, it perform better than KOREN limiter as shown in Reference [9]. Thus, this limiter is recommended for when choosing third-order limiters.

'venk': The famous VENKATAKRISHNAN limiter is popular in unstructured meshes. It is a little bit similar with the VANALBADA limiter, see Reference [10] for details. This limiter is claimed to keep a better balance between accuracy and oscillations, and thus can provide some more details than other second-order limtiters. The parameter K, as mentioned in their paper, is set to be 0.3. Other values like 1 or 10 are also widely used in other literatures.

'ppm': A four-point centered stencil third-order TVD limiter stands for the Piecewise Parabolic Method. See Reference [11] for details and see the references in mod_ppm.t for the implement of this method. Since 4 ghostcells is needed, it would be time-consuming as a third-order limiter.

'mp5': Limiters afterwards are not TVD any more. MP5 is a high order (fifth-order) five-point Monolicity Preserving limiter, which preserves the so-called MP property. see Reference [12] for details. MP5 can provide more details than other limiters provided in AMRVAC, but since this MP property cannot prevent small oscillations, even with small CFL number(0.4 or 0.2), it will crash in many cases. Note that the parameter alpha is chosen as 4, which is recommended in Reference [12]. Therefore, in principle, CFL number should be less than 1 / (1 + alpha) = 0.2. But in practice, 0.4 still yields nonoscillatory results.

'weno3', 'weno5', 'weno7': The classic WENO or Weighted Essentially Non(-)oscillatory Scheme, a popular variation of ENO scheme, which means that we allow kind of oscillations that is not allowed in the TVD scheme, but these oscillations could not be too large. This idea is similar with the so-called TVB scheme, where B stands for bounded. Check Reference [13] for the original idea of WENO scheme. The so-called classic WENO or WENO-JS scheme is from Reference [14], who modified it into a more practical way. High-order WENO scheme are especially suitable for problems containing both shocks and a large number of complex smooth structures. 'weno3', 'weno5', 'weno7' stand for third-order, fifth-order and seventh-order WENO-JS scheme, respectively. As mentioned in Reference [15], the coefficients named 'reconstruction' are used here. And in 'weno5', we also provide annother set of coefficients named 'interpolation', as applied in the ECHO code (see Reference [16]). Check mod_weno.t if it is necessary for you to use this set of coefficients. The parameter delta^2, which has a similar usage with van Albada and Venkatakrishnan limiter, is set to be 10^-18 by default. Values as low as 10^-40 is also possible Change it if your grid size is very small. For the performance, although the WENO scheme will allow sort of oscillation, it performs to be more diffusive than other third-order or fifth-order schemes, and thus robust in many cases.

'wenoz5': fifth-order WENO-Z scheme, see Reference [17] for details. One variation of the classic WENO scheme which is innovated by WENO-M but more efficient than WENO-M. It claims to perform better at critical points than WENO-JS, but is not as robust as WENO-JS in some extrme cases. Although the power index is set to be 1 in Reference [17], we still set it to be 2 so that it works in a similar way with WENO5 and WENOZ+5. You can change it in mod_weno.t if necessary. Apart from the common used value 2, value from 1 to 6 could all be found in previous papers.

'wenozp5': fifth-order WENO-Z+ scheme, see Reference [18] for details. An improved version of WENO-Z.

'mpweno7': seventh-order MP-WENO scheme, see Reference [19]. It is basically a combination of MP5 and WENO scheme to preserve the MP property in the results of, especially high-order WENO schemes, namely seventh-order of even higher. The advantage might be that it is more stable than WENO7. However, it shows two drawbacks in the tests, one is this shceme is time-comsuming, and another is the result is more diffusive, which is a straighforward result from the design of this scheme. Maybe efficiency could be improved here, but anyway, if WENO5 or WENO7 could be used properly in your case, MPWENO7 is not recommeneded.

'exeno7': seventh-order extended ENO scheme, also a variation of ENO scheme. It is based on the idea of Targeted ENO or TENO scheme which uses the smoothness measurement in the WENO scheme as a shock detector. But right now the this scheme is only used for tests because the choice of the crucial parameter C_T is still under testing. See Reference [20] for details.

References
1. Yee, H., 1989, A Class of High-Resolution Explicit and Implicit Shock-Capturing Methods.
2. Roe, P., 1983, Some Contributions to the Modelling of Discontinuous Flows.
3. Roe, P., 1986, Characteristic-Based schemes for the Euler Equations.
4. van Leer, B., 1977, Towards the Ultimate Conservative Difference Scheme. IV: A New Approach to Numerical Convection.
5. Woodward, P., 1984, The Numerical Simulation of Two-Dimensional Fluid Flow with Strong Shock.
6. van Leer, B., 1974. Towards the Ultimate Conservative Difference Scheme, II: Monotonicity and Conservation Combined in a Second-Order Scheme.
7. van Albada, G., 1982, A Comparative Study of Computional Methods in Cosmic Gas Dynamic.
8. Koren, B., 1993, A Robust Upwind Discretisation Method for Advection, Diffusion and Source Terms.
9. Cada, M., 2009, Compact Third-Order Limiter Functions for Finite Volume Methods.
10. Venkatakrishnan, V., 1993, On the Accuracy of Limiters and Convergence to Steady State Solutions.
11. Colella, P., 1984, The Piecewise Parabolic Method (PPM) for Gas-Dynamical Simulations.
12. Suresh, A., 1997, Accurate Monotonicity-Preserving Schemes with Runge–Kutta Time Stepping.
13. Liu, X., 1994, Weighted Essentially Non-Oscillatory Schemes.
14. Jiang, G., 1996, Efficient Implementation of Weighted ENO Schemes.
15. Shu, C., 2009, High Order Weighted Essentially Nonoscillatory Schemes for Convection Dominated Problems.
16. Del Zanna, L., 2007, ECHO: a Eulerian Conservative High-Order Scheme for General Relativistic Magnetohydrodynamics and Magnetodynamics.
17. Borges, R., 2008, An Improved Weighted Essentially Non-Oscillatory Scheme for Hyperbolic Conservation Laws.
18. Acker, F., 2016, An Improved WENO-Z Scheme.
19. Balsara, D., 2000, Monotonicity Preserving Weighted Essentially Non-Oscillatory Scheme with Increasingly High Order of Accuracy.
20. Xu, C., 2019, Arbitrary high-order extended ENO schemes for hyperbolic conservation laws.
