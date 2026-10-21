# Dumping arbitrary fields for debugging {#debug_field_dump}

[TOC]

# Motivation {#debug_field_dump_motivation}

The standard output routine `saveamrfile` can only write the registered set of
conservative variables and the `nwauxio` auxiliary fields. Many quantities of
interest during debugging are computed *inside* a block loop and never stored:
the fluxes `fC`, the constrained-transport electric field `fE`, source terms, a
reconstructed primitive state, or an intermediate Runge-Kutta stage. These are
invisible to the normal IO.

The field-dump tool lets you stash **any array, from anywhere inside a block
loop**, into a per-block scratch buffer and then flush it to a standard
@ref fileformat.md ".dat" file that is read by `amrvac_pytools` and yt without
any special handling.

The tool is opt-in: it is disabled by default (`wdebug_on = .false.`) and has no
runtime cost until you enable it. In use, you add capture lines temporarily while
chasing a problem and remove those lines afterwards; the tool itself is a
permanent part of the IO layer.

# How it works {#debug_field_dump_design}

The tool deliberately separates two phases:

* **Capture** is per-block and local, so it needs no synchronisation. Each block
  carries its own buffer `ps(igrid)%wdebug(ixG^T, n_wdebug)`, which you fill in
  place wherever the quantity is live.
* **Flush** reuses the existing collective snapshot writer at an
  MPI-coordinated point. It walks the full forest, so the result is AMR-correct
  and is a genuine `.dat` file.

The practical rule that follows is: **write inside the loop, flush outside it.**

The pieces are:

| Symbol | Defined in | Role |
|---|---|---|
| `wdebug` (field of the `state` type) | `mod_physicaldata.t` | per-block buffer, lazily allocated |
| `wdebug_on`, `n_wdebug` | `mod_physicaldata.t` | global enable flag (default `.false.`) and slot count |
| `debug_alloc(n, names)` | `mod_input_output.t` | register `n` named slots and set `wdebug_on=.true.` |
| `save_wdebug(suffix)` | `mod_input_output.t` | collective flush of `ps(:)%wdebug` to a `.dat` |

# Usage {#debug_field_dump_usage}

## Step 1: enable and name the slots

Call this once, before the region you want to inspect (for example just before
the time loop in `amrvac.t`). The names become the variable names in the output
file.

```{fortran}
use mod_input_output, only: debug_alloc, save_wdebug
...
call debug_alloc(3, [character(len=name_len) :: 'fE1', 'fE2', 'fE3'])
```

## Step 2: capture, inside a block loop

Write the field into the current block's buffer where it is computed. The global
`block` pointer refers to the persistent `ps(igrid)`, so the data survives to the
flush. Gate the capture on `wdebug_on` so the line is a no-op in production.

```{fortran}
if (wdebug_on) then
   if (.not. allocated(block%wdebug)) then
      allocate(block%wdebug(ixG^T, n_wdebug))
      block%wdebug = 0.0d0
   end if
   block%wdebug(ixI^S, 1)   = my_field(ixI^S)   ! slot 1
   block%wdebug(ixI^S, 2:3) = fE(ixI^S, 1:2)    ! slots 2-3
end if
```

Face-centred quantities (fluxes, `fE`) are simply stored at the host-cell index;
cell-averaging is not needed for localisation work.

## Step 3: flush at a barrier-safe point

Call `save_wdebug` only where every rank arrives together (between substeps, at
the end of a step) and **never** inside a block loop.

```{fortran}
if (wdebug_on) call save_wdebug('_fE')
```

This writes `<base_filename>_fE NNNN.dat`, where `NNNN` is the current iteration
counter `it`.

# Reading the output {#debug_field_dump_reading}

The file is an ordinary cell-centred `.dat`:

```{python}
import amrvac_pytools as apt
from amrvac_pytools.datfiles.reading import datfile_utilities as D
r = apt.load_datfile("run_fE0000.dat")
u = D.get_uniform_data(r)        # shape [nx, ny, nz, n_wdebug]
print(r.header['w_names'])       # the names passed to debug_alloc
```

# Inspecting ghost cells {#debug_field_dump_ghosts}

`save_wdebug` writes only the **interior** of each block; `get_uniform_data`
never exposes ghost cells. To compare a boundary ghost layer (for example at the
lower z boundary), shift the source index down by `nghostcells` so the ghost
layers land in interior slots:

```{fortran}
ps(igrid)%wdebug(ixMlo1:ixMhi1, ixMlo2:ixMhi2, ixMlo3:ixMhi3, 1:nw) = &
   ps1(igrid)%w(ixMlo1:ixMhi1, ixMlo2:ixMhi2, &
                ixMlo3-nghostcells:ixMhi3-nghostcells, 1:nw)
```

The bottom rows of the resulting uniform array are then the physical ghost cells.

# Notes and pitfalls {#debug_field_dump_notes}

* The flush is collective; place `save_wdebug` where all ranks arrive together.
* Only interior cells are written; use the index shift above to inspect ghosts.
* Capture into the persistent `ps(igrid)` (which the global `block` points to in
  the advance path). Writing into a temporary state's buffer will not be seen by
  `save_wdebug`, which reads `ps(:)%wdebug`.
* The header is written by reusing `snapshot_write_header1` with a temporary
  `stagger_grid = .false.` toggle; do not hand-roll a header.
* A subroutine-local `use mod_input_output, only: saveamrfile` shadows the
  module-level import. Add `debug_alloc, save_wdebug` to that `only` list. The
  module `mod_advance` may `use mod_input_output` directly (no circular
  dependency).
* If `n_wdebug` is zero, `save_wdebug` returns immediately, but a capture line
  that writes slots `1:k` will go out of bounds. Always gate captures on
  `wdebug_on`.

# Production safety {#debug_field_dump_safety}

The tool ships in the build but is inert at runtime: `wdebug_on` defaults to
`.false.` and the per-block `wdebug` buffers are allocated only once `debug_alloc`
is called, so a run that never calls it pays no memory or time cost. Only the
temporary capture lines you add for a given investigation need to be removed
afterwards; the `debug_alloc` / `save_wdebug` infrastructure stays in place.
