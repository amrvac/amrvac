title: Coordinate systems

# Available geometries {: #geom_overview }

AGILE supports five coordinate settings: Cartesian (x, y, z), spherical
(r, θ, φ), cylindrical (r, z, φ), and the two logarithmic-radius variants
`logSpherical` and `logCylindrical` (the same two curvilinear systems with a
radial mesh that is uniform in ln(r) rather than in r -- see [Logarithmically
stretched radius](#geom_log_radius) below). All are always 3D, so a
curvilinear run is the full 3D system, never an axisymmetric or 1D/2D
reduction.

The coordinate system is a **compile-time** choice: the finite-volume
kernels are specialized for it, so that the (much cheaper) Cartesian kernel
stays free of curvilinear bookkeeping. It is selected via `geometry` in
`&meshlist` of the parameter file:

```fortran
&meshlist
   geometry = 'spherical'   ! or 'cylindrical', 'logSpherical',
                            ! 'logCylindrical', or 'Cartesian' (the default)
/
```

`logSpherical`/`logCylindrical` are a spherical/cylindrical build (`GEOM`,
and the runtime `coordinate`, stay `'spherical'`/`'cylindrical'`) plus a
`LOG_RADIUS` flag; every geometry-specific branch stays correct untouched and
only the radial spacing differs.

Conventions to know when setting up a curvilinear case:

* **Angular domain bounds are given in units of 2π.** So for spherical,
  `xprobmin2 = 0.125d0`, `xprobmax2 = 0.375d0` places θ between π/4 and
  3π/4. For cylindrical, only the φ bound (`xprobmin3`/`xprobmax3`) gets
  this treatment; the z bound (`xprobmin2`/`xprobmax2`) is a plain length,
  exactly like r.
* **Axis order** is (r, θ, φ) for spherical and (r, z, φ) for cylindrical --
  in both cases φ is the third coordinate direction.
* Under `logSpherical`/`logCylindrical`, `xprobmin1`/`xprobmax1` are
  **physical radii** in the parameter file and are converted to the logical
  ln-radius coordinate by `read_par_files`. A case's `mod_usr.fpp` that reads
  those bounds for a physical radius must map them back with `r_of_s` from
  `mod_geometry` (not by hand -- the map is only a plain exponential when
  `log_r0 = 0`).
* **Cell positions `ps(igrid)%x` are the volume barycentre of the cell, not
  the midpoint of its faces** -- so `x ± dx/2` is not a cell face. See [Cell
  positions are volume barycentres](#geom_barycentre) below.

## Cell positions are volume barycentres {: #geom_barycentre }

`ps(igrid)%x` -- which every case's `mod_usr.fpp` reads -- holds the volume
barycentre of the cell. A cell average equals the point value at the
barycentre to second order (the linear term of the Taylor expansion
integrates to zero only about that point), so the barycentre is where the
state, the initial condition and any position-dependent source term belong.
Only the directions with a non-uniform volume weight are affected: r in both
systems, and θ in spherical. φ in both, and z in cylindrical, carry uniform
weight, so their barycentres are their midpoints. The offset is about 1% of a
cell width radially on the test grids and grows to dθ/6 -- a sixth of a cell,
pushed away from the axis -- for the first polar cell off θ = 0.

<figure class="doc-figure">
<svg viewBox="0 0 720 288" role="img" aria-label="A radial cell drawn unrolled as a wedge that is narrow at the inner face r_L and wide at the outer face r_R; the face midpoint sits at the centre in r, but the volume barycentre sits further out because more of the cell's volume lies on the outer side">
<defs>
<marker id="ah2" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="8" markerHeight="8" orient="auto-start-reverse">
<path d="M0 0L10 5L0 10z" fill="var(--fig-accent,#B0501F)"/>
</marker>
</defs>
<path d="M 150 155 L 560 96 L 560 244 L 150 191 Z" fill="var(--fig-diagram-bg,#DDEBEC)" stroke="var(--fig-diagram,#2B6C7C)" stroke-width="1.6"/>
<line x1="150" y1="155" x2="150" y2="191" stroke="var(--fig-diagram,#2B6C7C)" stroke-width="3.5"/>
<line x1="560" y1="96" x2="560" y2="244" stroke="var(--fig-diagram,#2B6C7C)" stroke-width="3.5"/>
<text x="150" y="214" text-anchor="middle" font-family="monospace" font-size="12" fill="currentColor">r_L</text>
<text x="560" y="262" text-anchor="middle" font-family="monospace" font-size="12" fill="currentColor">r_R</text>
<text x="120" y="140" font-family="monospace" font-size="9.5" fill="var(--fig-faint,#8B98A1)">less volume</text>
<text x="520" y="80" font-family="monospace" font-size="9.5" fill="var(--fig-faint,#8B98A1)">more volume</text>
<line x1="150" y1="173" x2="560" y2="170" stroke="currentColor" stroke-width="0.75" stroke-dasharray="2 3"/>
<circle cx="355" cy="171.5" r="6" fill="var(--fig-surface,#fff)" stroke="currentColor" stroke-width="1.6"/>
<text x="355" y="120" text-anchor="middle" font-family="monospace" font-size="11" fill="var(--fig-soft,#566470)">r_c = 1/2(r_L + r_R)</text>
<text x="355" y="105" text-anchor="middle" font-family="monospace" font-size="9.5" fill="var(--fig-faint,#8B98A1)">face midpoint — the old stored x</text>
<line x1="355" y1="126" x2="355" y2="164" stroke="currentColor" stroke-width="0.75"/>
<circle cx="404" cy="170.5" r="6.5" fill="var(--fig-accent,#B0501F)"/>
<path d="M 355 194 Q 380 205 404 194" fill="none" stroke="var(--fig-accent,#B0501F)" stroke-width="1.5" marker-end="url(#ah2)"/>
<text x="380" y="230" text-anchor="middle" font-family="monospace" font-size="11" fill="var(--fig-accent,#B0501F)">r-bar — volume barycentre, +O(h²) outward</text>
<text x="380" y="246" text-anchor="middle" font-family="monospace" font-size="10" fill="var(--fig-soft,#566470)">spherical: r-bar = 3/4·(r_R⁴ − r_L⁴)/(r_R³ − r_L³)</text>
<line x1="404" y1="177" x2="404" y2="188" stroke="var(--fig-accent,#B0501F)" stroke-width="1"/>
<text x="360" y="279" text-anchor="middle" font-family="monospace" font-size="10.5" fill="var(--fig-soft,#566470)">shift ≈ 1% of dr radially  ·  grows to dθ/6 for the first cell off θ = 0</text>
</svg>
<figcaption>
<b>The point moves outward to the centre of volume.</b> Only directions with a non-uniform volume weight are affected: the radius in both systems, and <code>θ</code> in spherical. <code>φ</code> (both) and cylindrical <code>z</code> carry uniform weight, so their barycentres already were their midpoints. <code>fill_geometry_device</code> evaluates the offset forms with the common factor cancelled analytically, because computing an <code>O(h²)</code> quantity as the difference of two <code>O(1)</code> ones would throw away most of the mantissa on a fine grid.
</figcaption>
</figure>

`fill_geometry_device` is written face-first: the radial-cell macro produces
the two faces, the midpoint and the extent, and every area and volume is then
written in terms of those so the expressions stay algebraically exact for any
face pair. Three places in the tree that assumed `x ± dx/2` was a face -- the
face coordinates handed to the Riemann solver in `mod_finite_volume`, the VTU
corner grid in `calc_x`, and `fill_geometry_device` itself -- rebuild the face
from `rnode` instead. `bgeo%ds` deliberately stays on the cell midpoint: it is
a geometric extent that `setdt` turns into the CFL length, and following the
barycentre there would quietly relax the timestep.

## Why curvilinear coordinates need extra source terms {: #geom_why_sources }

AGILE discretizes the generic conservation law

\begin{equation}
    \partial_t \mathbf{U} + \nabla \cdot \mathbf{F}(\mathbf{U}) = \mathbf{S}(\mathbf{U},\nabla \mathbf{U},\mathbf{x},t)
\end{equation}

(see [Physics modules and equations](equations.html)) with a finite-volume
update that divides the net flux through a cell's faces by the cell volume.
In spherical or cylindrical coordinates the face areas and cell volume
depend on position (e.g. a spherical shell's outer face is larger than its
inner face), so this ratio is what actually implements the curvilinear
divergence for scalar quantities -- no extra term is needed for the mass,
energy, or (for `ffhd`) the field-aligned quantities, since these are
advected as genuine scalars and Gauss's theorem then guarantees the finite
volume update alone reproduces \( \nabla \cdot \mathbf{F} \) correctly in
any coordinate system.

Momentum (and, for `mhd` and `srmhd`, the magnetic field) is different: it is a
*vector* quantity, and its flux is a tensor built from a position-dependent
(r, θ, φ) or (r, z, φ) basis. Writing the vector Laplacian/divergence out in
curvilinear coordinates leaves over extra terms beyond the flux-divergence
that the finite-volume update captures -- physically, these account for how
the basis vectors themselves change from one face of a cell to the other.
AGILE adds these leftover terms explicitly as a per-cell geometric source,
in the `addsource_geometry` routine of each physics module
(`src/{hd,mhd,srhd,srmhd,ffhd}/mod_*_templates.fpp`); a physics module that does
not implement one aborts at startup rather than silently omitting these
terms. The explicit expressions below are exactly what these routines
evaluate (rewritten into physical notation matching
[Physics modules and equations](equations.html); see the source files for
the exact discretization, in particular how the isotropic-pressure-like
prefactor is evaluated from the actual discrete face areas rather than
plugging in the continuous coordinate factor directly, so that a uniform
pressure produces no net force to the same precision the code otherwise
resolves -- for cylindrical this recovers the continuous prefactor exactly,
since the radial face area of a cylindrical cell is linear in r). Because
those discrete factors are the *exact* volume averages of the continuous
`2/r`, `cot(θ)/r` (spherical) or `1/r` (cylindrical) prefactors for any face
pair, no curvature term depends on where in the cell the stored position
sits, and `addsource_geometry` takes no position argument -- which is what
keeps it correct now that `ps%x` is the volume barycentre rather than the
face midpoint.

# Geometric source terms {: #geom_sources }

The equations below give the additional source term \( \mathbf{S}_{\mathrm{geom}} \)
that a curvilinear run adds to the corresponding equation from
[Physics modules and equations](equations.html); everywhere else the
equations are unchanged. `ffhd` needs none, for any geometry (see above), so
it is omitted below beyond noting \( \mathbf{S}_{\mathrm{geom}} = 0 \).

## Spherical (r, θ, φ) {: #geom_spherical_sources }

### Hydrodynamics: hd {: #geom_spherical_hd }

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{2p+\rho\left(v_\theta^2+v_\varphi^2\right)}{r} \,,\\
    S_{\mathrm{geom}}[m_\theta] &=& \frac{p\cot\theta+\rho\left(v_\varphi^2\cot\theta-v_r v_\theta\right)}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& -\frac{\rho v_\varphi\left(v_r+v_\theta\cot\theta\right)}{r} \,.
\end{eqnarray}

### Magnetohydrodynamics: mhd {: #geom_spherical_mhd }

With the total (thermal + magnetic) pressure \( p_{\mathrm{tot}} = p + \frac{1}{2}B^2 \):

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{2p_{\mathrm{tot}}+\rho\left(v_\theta^2+v_\varphi^2\right)-\left(B_\theta^2+B_\varphi^2\right)}{r} \,,\\
    S_{\mathrm{geom}}[m_\theta] &=& \frac{p_{\mathrm{tot}}\cot\theta+\cot\theta\left(\rho v_\varphi^2-B_\varphi^2\right)-\rho v_r v_\theta+B_r B_\theta}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& \frac{-\rho v_\varphi\left(v_r+v_\theta\cot\theta\right)+B_\varphi\left(B_r+B_\theta\cot\theta\right)}{r} \,,\\
    S_{\mathrm{geom}}[B_r] &=& \frac{2\psi}{r} \,,\\
    S_{\mathrm{geom}}[B_\theta] &=& \frac{v_r B_\theta-v_\theta B_r+\psi\cot\theta}{r} \,,\\
    S_{\mathrm{geom}}[B_\varphi] &=& \frac{v_r B_\varphi-v_\varphi B_r+\cot\theta\left(v_\theta B_\varphi-v_\varphi B_\theta\right)}{r} \,,
\end{eqnarray}

where ψ is the GLM field (see [Divergence B source
treatments](equations.html#eq_divB_fix)); ψ itself has no geometric source.

### Special relativistic hydrodynamics: srhd {: #geom_spherical_srhd }

Structurally identical to the `hd` case, with the momentum density
\( \mathbf{m} = \xi\mathbf{v} \) (rather than \( \rho\mathbf{v} \)) and
using the spatial four-velocity components \( \mathbf{u}\equiv\Gamma\mathbf{v} \)
(the primitive momentum-slot variables; \( \mathbf{v}=\mathbf{u}/\Gamma \)),
so that \( \xi v_iv_j = (\xi/\Gamma^2)\,u_iu_j \):

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{2p+(\xi/\Gamma^2)\left(u_\theta^2+u_\varphi^2\right)}{r} \,,\\
    S_{\mathrm{geom}}[m_\theta] &=& \frac{p\cot\theta+(\xi/\Gamma^2)\left(u_\varphi^2\cot\theta-u_r u_\theta\right)}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& -\frac{(\xi/\Gamma^2)\,u_\varphi\left(u_r+u_\theta\cot\theta\right)}{r} \,,
\end{eqnarray}

with \( \xi=\Gamma^2\rho h \) as in [Physics modules and
equations](equations.html#eq_srhd).

### Special relativistic magnetohydrodynamics: srmhd {: #geom_spherical_srmhd }

These are the curvature terms of the [srmhd](equations.html#eq_srmhd)
momentum flux, written with its stress tensor (pressure excluded),

\begin{equation}
    T_{ij} = W v_iv_j - (\mathbf{v}\cdot\mathbf{B})\left(B_iv_j+v_iB_j\right) - \frac{B_iB_j}{\Gamma^2} \,,
    \qquad W = \Gamma^2\rho h + B^2 \,,
\end{equation}

and the total pressure \( p_{\mathrm{tot}} = p + \frac{1}{2}\left(B^2/\Gamma^2+(\mathbf{v}\cdot\mathbf{B})^2\right) \).
Here \( \mathbf{v}=\mathbf{u}/\Gamma \) comes from the primitive four-velocity,
as in `srhd`. The terms are

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{2p_{\mathrm{tot}}+T_{\theta\theta}+T_{\varphi\varphi}}{r} \,,\\
    S_{\mathrm{geom}}[m_\theta] &=& \frac{\left(p_{\mathrm{tot}}+T_{\varphi\varphi}\right)\cot\theta-T_{r\theta}}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& -\frac{T_{r\varphi}+T_{\theta\varphi}\cot\theta}{r} \,,
\end{eqnarray}

The pressure gradient is part of the flux divergence, as for the other
modules. For **B** = 0 they reduce to the `srhd` terms above. For
Γ = 1, dropping the (**v**·**B**) terms, they reduce to the `mhd` terms. The induction and
ψ terms are the `mhd` ones evaluated with **v** = **u**/Γ:

\begin{eqnarray}
    S_{\mathrm{geom}}[B_r] &=& \frac{2\psi}{r} \,,\\
    S_{\mathrm{geom}}[B_\theta] &=& \frac{v_r B_\theta-v_\theta B_r+\psi\cot\theta}{r} \,,\\
    S_{\mathrm{geom}}[B_\varphi] &=& \frac{v_r B_\varphi-v_\varphi B_r+\cot\theta\left(v_\theta B_\varphi-v_\varphi B_\theta\right)}{r} \,.
\end{eqnarray}

## Cylindrical (r, z, φ) {: #geom_cylindrical_sources }

Only the \(m_r\)/\(m_\varphi\) (and, for `mhd` and `srmhd`, \(B_r\)/\(B_\varphi\))
components pick up a geometric source; \(m_z\)/\(B_z\) do not, since a
cylindrical volume element's z-extent does not depend on r.

### Hydrodynamics: hd {: #geom_cylindrical_hd }

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{p+\rho v_\varphi^2}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& -\frac{\rho v_\varphi v_r}{r} \,.
\end{eqnarray}

### Magnetohydrodynamics: mhd {: #geom_cylindrical_mhd }

With \( p_{\mathrm{tot}} = p + \frac{1}{2}B^2 \):

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{p_{\mathrm{tot}}-B_\varphi^2+\rho v_\varphi^2}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& \frac{-\rho v_\varphi v_r+B_\varphi B_r}{r} \,,\\
    S_{\mathrm{geom}}[B_r] &=& \frac{\psi}{r} \,,\\
    S_{\mathrm{geom}}[B_\varphi] &=& \frac{B_\varphi v_r-B_r v_\varphi}{r} \,.
\end{eqnarray}

### Special relativistic hydrodynamics: srhd {: #geom_cylindrical_srhd }

As for spherical, with \( \mathbf{u}\equiv\Gamma\mathbf{v} \) and
\( \xi=\Gamma^2\rho h \):

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{p+(\xi/\Gamma^2)\,u_\varphi^2}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& -\frac{(\xi/\Gamma^2)\,u_\varphi u_r}{r} \,.
\end{eqnarray}

### Special relativistic magnetohydrodynamics: srmhd {: #geom_cylindrical_srmhd }

With the same \( T_{ij} \) and \( p_{\mathrm{tot}} \) as for spherical:

\begin{eqnarray}
    S_{\mathrm{geom}}[m_r] &=& \frac{p_{\mathrm{tot}}+T_{\varphi\varphi}}{r} \,,\\
    S_{\mathrm{geom}}[m_\varphi] &=& -\frac{T_{r\varphi}}{r} \,,\\
    S_{\mathrm{geom}}[B_r] &=& \frac{\psi}{r} \,,\\
    S_{\mathrm{geom}}[B_\varphi] &=& \frac{B_\varphi v_r-B_r v_\varphi}{r} \,.
\end{eqnarray}

# Logarithmically stretched radius {: #geom_log_radius }

`geometry = 'logSpherical'` and `geometry = 'logCylindrical'` place the
radial mesh uniformly in ln(r + r0) rather than in r. The offset `r0` is the
`log_r0` entry of `&meshlist` and defaults to zero, the plain ln(r): the cell
width then grows in proportion to r, so the relative resolution dr/r is
constant -- the one stretching astrophysics usually wants, and the only one
supported.

The design principle is that taking the *logical* radial coordinate to be
ξ = ln(r + r0) leaves the mesh uniform in ξ, which is exactly what the
geometry kernel already assumes -- an analytic function of the block corner
and a constant spacing. The tree, the AMR level ladder, `rnode`, prolongation,
coarsening and the ghost exchange all keep working unchanged in ξ; only the
step from logical faces to physical ones takes an `exp`. It needs none of
upstream MPI-AMRVAC's general stretched-grid machinery (`qstretch`, `dxfirst`,
per-level stretch tables) -- which AGILE has removed -- and refinement and log
stretching commute for free.

A positive `log_r0` moves the singularity of the map from r = 0 out to
r = −r0, making the mesh genuinely uniform (cell width r0·dξ) for r ≪ r0 and
logarithmic for r ≫ r0. That is what lets a stretched radial grid **reach
r = 0**: with the default `log_r0 = 0` a log-radial grid can never reach the
axis, so `xprobmin1 ≤ 0` is rejected unless `log_r0 > 0`. `logCylindrical`
with `log_r0 > 0` therefore supports the cylindrical axis; `logSpherical`
still has no usable origin (nothing fills the ghost cells across r = 0).

<figure class="doc-figure">
<svg viewBox="0 0 720 330" role="img" aria-label="Three radial meshes on the same axis: uniform in r; uniform in ln(r), which cannot reach r=0; and the offset map uniform-then-logarithmic, whose pivot at s=0 is exactly r=0 and whose ghost cells mirror the interior">
<defs>
<marker id="ah3" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
<path d="M0 0L10 5L0 10z" fill="currentColor"/>
</marker>
</defs>
<text x="20" y="40" font-family="monospace" font-size="11" fill="var(--fig-soft,#566470)">uniform in r</text>
<line x1="150" y1="46" x2="690" y2="46" stroke="currentColor" stroke-width="1"/>
<g stroke="var(--fig-diagram,#2B6C7C)" stroke-width="1.4">
<line x1="150" y1="38" x2="150" y2="54"/><line x1="204" y1="38" x2="204" y2="54"/>
<line x1="258" y1="38" x2="258" y2="54"/><line x1="312" y1="38" x2="312" y2="54"/>
<line x1="366" y1="38" x2="366" y2="54"/><line x1="420" y1="38" x2="420" y2="54"/>
<line x1="474" y1="38" x2="474" y2="54"/><line x1="528" y1="38" x2="528" y2="54"/>
<line x1="582" y1="38" x2="582" y2="54"/><line x1="636" y1="38" x2="636" y2="54"/>
<line x1="690" y1="38" x2="690" y2="54"/>
</g>
<text x="150" y="70" text-anchor="middle" font-family="monospace" font-size="10" fill="var(--fig-faint,#8B98A1)">r=0</text>
<text x="20" y="140" font-family="monospace" font-size="11" fill="var(--fig-soft,#566470)">uniform in ln(r)</text>
<line x1="150" y1="146" x2="690" y2="146" stroke="currentColor" stroke-width="1"/>
<g stroke="var(--fig-diagram,#2B6C7C)" stroke-width="1.4">
<line x1="182" y1="138" x2="182" y2="154"/><line x1="206" y1="138" x2="206" y2="154"/>
<line x1="238" y1="138" x2="238" y2="154"/><line x1="282" y1="138" x2="282" y2="154"/>
<line x1="342" y1="138" x2="342" y2="154"/><line x1="424" y1="138" x2="424" y2="154"/>
<line x1="536" y1="138" x2="536" y2="154"/><line x1="690" y1="138" x2="690" y2="154"/>
</g>
<line x1="150" y1="130" x2="150" y2="162" stroke="var(--fig-drop,#9A5350)" stroke-width="1.5" stroke-dasharray="3 3"/>
<text x="150" y="182" text-anchor="middle" font-family="monospace" font-size="10" fill="var(--fig-drop,#9A5350)">r=0 unreachable</text>
<text x="176" y="130" font-family="monospace" font-size="10" fill="var(--fig-faint,#8B98A1)">r_min &gt; 0</text>
<text x="20" y="242" font-family="monospace" font-size="11" fill="var(--fig-soft,#566470)">s = ln(1 + r/r0)</text>
<line x1="150" y1="248" x2="690" y2="248" stroke="currentColor" stroke-width="1"/>
<g stroke="var(--fig-accent,#B0501F)" stroke-width="1.4" stroke-dasharray="2 2">
<line x1="366" y1="240" x2="366" y2="256"/><line x1="338" y1="240" x2="338" y2="256"/>
<line x1="306" y1="240" x2="306" y2="256"/>
</g>
<circle cx="394" cy="248" r="4.5" fill="var(--fig-accent,#B0501F)"/>
<text x="394" y="278" text-anchor="middle" font-family="monospace" font-size="10" fill="var(--fig-accent,#B0501F)">s = 0  =  r = 0  (exactly)</text>
<g stroke="var(--fig-diagram,#2B6C7C)" stroke-width="1.4">
<line x1="394" y1="240" x2="394" y2="256"/><line x1="422" y1="240" x2="422" y2="256"/>
<line x1="452" y1="240" x2="452" y2="256"/><line x1="486" y1="240" x2="486" y2="256"/>
<line x1="528" y1="240" x2="528" y2="256"/><line x1="582" y1="240" x2="582" y2="256"/>
<line x1="652" y1="240" x2="652" y2="256"/>
</g>
<path d="M 366 224 A 40 40 0 0 1 422 224" fill="none" stroke="var(--fig-accent,#B0501F)" stroke-width="1.3" marker-end="url(#ah3)" marker-start="url(#ah3)"/>
<text x="394" y="214" text-anchor="middle" font-family="monospace" font-size="9.5" fill="var(--fig-accent,#B0501F)">mirror: odd about s = 0</text>
<text x="470" y="304" text-anchor="middle" font-family="monospace" font-size="10" fill="var(--fig-soft,#566470)">uniform for r ≪ r0   ·   logarithmic for r ≫ r0</text>
</svg>
<figcaption>
<b>Reaching the axis.</b> With the default <code>log_r0 = 0</code> the map is plain <code>ln(r)</code> and cannot reach <code>r = 0</code>. A positive <code>log_r0</code> moves the singularity out to <code>r = −r0</code>, so the mesh is genuinely uniform near the axis and logarithmic far out. Three properties of that form are deliberate: the <b>pivot is the axis</b> (<code>s = 0</code> is <code>r = 0</code> exactly, so <code>rnode</code>'s corner is exactly zero); the <b>map is odd</b>, so the pole copy fills a ghost cell from an interior cell of <em>mirrored volume</em>; and the <b>odd branch is confined to the ghost layer</b>, so the physical domain still uses the cheaper <code>log_ra·exp(s) + log_rb</code> and a plain <code>log_r0 = 0</code> build is bit-for-bit unchanged.
</figcaption>
</figure>

Two limitations of the stretching:

* **The reconstruction stays uniform-stencil.** The MUSCL reconstruction and
  the limiters assume equal cells, so the second-order constant degrades as
  the stretch ratio grows. Keep it modest -- the test cases run at ratios
  around 1.04 to 1.09 per cell. A positive `log_r0` helps here, since the mesh
  is genuinely uniform in the region next to the axis where the stencil
  assumption would otherwise be worst.
* **The timestep tightens.** dr = (r + r0)·dξ, so a small `log_r0` chosen for
  near-axis resolution buys a CFL restriction on top of the one the polar axis
  already imposes.

# The polar axis {: #geom_polar_axis }

A spherical domain may span the full θ from 0 to π, and a cylindrical one may
start at r = 0. The axis is not a boundary: a ghost cell just past the axis at
azimuth φ *is* the interior cell at φ + π, so it is filled from the block half
a revolution away (a "pi-periodic" treatment) rather than by a boundary
condition. Setting it up takes three things in `&meshlist`/`&boundlist`, and
`check_pole_setup` stops the run with a specific message if any is missing:

* `'pole'` on the axis face(s) -- `typeboundary_min2`/`max2` in spherical,
  `typeboundary_min1` in cylindrical; it is all-or-nothing per face.
* φ periodic over exactly one full turn (`xprobmin3 = 0`, `xprobmax3 = 1` in
  2π units, `periodic` at both φ faces).
* an even number of level-1 blocks across φ.

The sign rule is the same in both systems: the component along the mirrored
direction flips, the φ component flips, everything else is symmetric. Note
that this makes the **radial** component odd at the cylindrical axis
(m_r, B_r as well as m_φ, B_φ) -- upstream MPI-AMRVAC marks only the φ
component antisymmetric there, which is wrong.

The polar axis is supported for `hd`, `mhd`, `srhd`, `srmhd` and `ffhd`. Two things to
know:

* **The timestep collapses near the axis.** ds(φ) = r·sin(θ)·dφ → 0 there, so
  `setdt` gives a severe CFL restriction. No axis filtering or ring averaging
  is implemented; keep pole cases at modest resolution.
* **Volume averages cannot validate the pole copy** -- a cell on the axis has
  vanishing volume and a zero-area axis face, so a wrong ghost value barely
  moves `mean(w)`. The pole test cases therefore assert on the ghost cells
  directly against the analytic state at `it = 0`.

# Current limitations {: #geom_limitations }

* **AMR across curvilinear refinement levels is only lightly tested.**
  Prolongation and coarsening volume-weight correctly, but `fix_conserve` is
  disabled for all geometries in `src/mod_advance.fpp`.
* **General grid stretching is unavailable.** A `stretch_dim` that varied the
  cell spacing within a block could not be expressed by the analytic metric
  derivation, so upstream MPI-AMRVAC's stretched-grid machinery was removed and
  the `stretch_dim` / `qstretch_baselevel` / `nstretchedblocks_baselevel` /
  `stretch_uncentered` keys are no longer read from `&meshlist`. The
  logarithmic radius above is the one supported non-uniform mesh, precisely
  because it is still uniform in ξ = ln(r + r0).
* `divvector`/`curlvector` in `src/mod_geometry.fpp` still assume a face lies
  midway between two cell centres, which is untrue for the volume barycentre.
  Neither routine is on a path this fork currently exercises.
