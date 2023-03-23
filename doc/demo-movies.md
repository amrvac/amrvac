# MPI-AMRVAC 3.0 Demo simulations

# Demos {#demo_top}

A set of Demo simulations, all found in our tests/demo folder

* @ref demo_vac 2D Advection and particle sampling
* @ref demo_rd 2D Reaction-Diffusion
* @ref demo_1d 1D hydro 
* @ref demo_kh 2D hydro Kelvin-Helmholtz
* @ref demo_khdust 2D hydro Kelvin-Helmholtz with gas-dust coupling
* @ref demo_ti 2D hydro thermal runaway
* @ref demo_alfven 2D MHD shock-cloud: Alfven hits Alfven
* @ref demo_tilt 2D resistive MHD: ideal tilt leads to chaotic reconnection
* @ref demo_part Charged particle motion in MHD fields

# 2D linear advection of the VAC-logo, and particle sampling {#demo_vac}

Our favorite density pattern in the VAC-logo shape gets diagonally advected, while we sample the solution in three locations.

\htmlonly
<video width="1444" height="916" controls>
  <source src="Fig19-VACpart.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
\endhtmlonly

@ref demo_top To top of page

# 2D Reaction-Diffusion Gray-Scott evolution {#demo_rd}

A 2D reaction-diffusion PDE system shows self-replicating spot formation.

\htmlonly
<video width="800" height="750" controls>
  <source src="Fig11_GrayScott.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
\endhtmlonly

@ref demo_top To top of page

# 1D hydro Shu-Osher test {#demo_1d}

This 1D hydro test shows a shock impinging on a sinusoidally varying density region. We demonstrate how varying reconstructions/limiters handle especially the compressed density variation trailing the right-ward moving shock.

\htmlonly
<video width="800" height="600" controls>
  <source src="Fig1-ShuOsher.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
\endhtmlonly

@ref demo_top To top of page

# 2D hydro Kelvin-Helmholtz evolution {#demo_kh}

A 2D pure hydro evolution of an initially sharp-interfaced shear layer: with the identical initial setup and only changing the reconstruction/limiting strategy, qualitative differences appear in resolving density fine-structure.

\htmlonly
<video width="875" height="750" controls>
  <source src="khi.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
\endhtmlonly

@ref demo_top To top of page

# 2D gas-dust coupled Kelvin-Helmholtz simulation {#demo_khdust}

A similar setup, but now with dynamically coupled gas-dust species, where the coupling strength parameter varies (top row to bottom row, increasing the coupling). Left panels show the gas evolutions, right panels the dust species.

\htmlonly
<video width="1450" height="950" controls>
  <source src="khi-dust.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
\endhtmlonly

@ref demo_top To top of page

# 2D hydro Thermal instability and runaway {#demo_ti}

This test shows how thermal instability driven by radiative losses in a 2D hydro setup can evolve towards a runaway condensation followed by a fragmentary, dynamically adjusting multiphase medium. Note how the slowly evolving gradual thermal imbalance suddenly changes into collapse and further erratic behavior.

\htmlonly
<video width="1086" height="920" controls>
  <source src="Fig5-TI.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
\endhtmlonly

@ref demo_top To top of page

# 2D MHD Alfven shock {#demo_alfven}

Here we let an Alfven (intermediate) shock hit a density variation containing Alfven's image. Three different flux schemes produce overall consistent end results, also involving a corrugated shock front.

\htmlonly
<video width="1600" height="550" controls>
  <source src="alfven-lr.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
\endhtmlonly

@ref demo_top To top of page

# 2D MHD Tilt instability {#demo_tilt}

Only varying the numeric magnetic monopole control in a 2D MHD test: we find consistent evolutions of an ideal MHD tilt instability, where the finite resistivity causes chaotic islands to form in the final reconnecting stages.

\htmlonly
<video width="1500" height="1250" controls>
  <source src="tilt.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
\endhtmlonly

@ref demo_top To top of page

# Charged particle motion {#demo_part}

In the central region of the tilt setup from above, the current density variation shows lots of fine-structure, and this could be an ideal site for particle acceleration. Here we show how we trace particle dynamics in a fixed snapshot from that simulation.

\htmlonly
<video width="1480" height="866" controls>
  <source src="Fig20-charged_part.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
\endhtmlonly

@ref demo_top To top of page
