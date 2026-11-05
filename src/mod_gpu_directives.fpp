!=====================================================================
! GPU offload directive macros.
!
! Every `!$acc`/`!$omp` directive in the codebase is written by calling
! one of the macros below, so that a single source line expands to the
! correct directive for whichever GPU backend is selected at build time
! (`make OPENACC=1` or `make OPENMP=1`.
!=====================================================================

#:if defined('OPENACC')

#:def GPU_DECLARE_CREATE(vars)
!$acc declare create(${vars}$)
#:enddef

#:def GPU_DECLARE_COPYIN(vars)
!$acc declare copyin(${vars}$)
#:enddef

#:def GPU_ROUTINE_SEQ()
!$acc routine seq
#:enddef

#:def GPU_ROUTINE_VECTOR()
!$acc routine vector
#:enddef

#:def GPU_ROUTINE()
!$acc routine
#:enddef

#:def GPU_PARALLEL_LOOP(clauses='')
!$acc parallel loop default(present) ${clauses}$
#:enddef

#:def GPU_PARALLEL_LOOP_GANG(clauses='')
!$acc parallel loop gang default(present) ${clauses}$
#:enddef

#:def GPU_LOOP_VECTOR(clauses='')
!$acc loop vector ${clauses}$
#:enddef

#:def GPU_LOOP_SEQ()
!$acc loop seq
#:enddef

#:def GPU_INDEPENDENT()
independent
#:enddef

#:def GPU_ENTER_DATA_COPYIN(vars)
!$acc enter data copyin(${vars}$)
#:enddef

#:def GPU_ENTER_DATA_CREATE(vars)
!$acc enter data create(${vars}$)
#:enddef

#:def GPU_EXIT_DATA_DELETE(vars)
!$acc exit data delete(${vars}$)
#:enddef

#:def GPU_EXIT_DATA_DELETE_IF_PRESENT(vars)
#:set depth = 0
#:set parts = ['']
#:for ch in vars
#:if ch == '('
#:set depth = depth + 1
#:endif
#:if ch == ')'
#:set depth = depth - 1
#:endif
#:if ch == ',' and depth == 0
#:set parts = parts + ['']
#:else
#:set parts = parts[:-1] + [parts[-1] + ch]
#:endif
#:endfor
#:set parts = [p.strip() for p in parts]
#:for part in parts
!$acc exit data delete(${part}$) if(acc_is_present(${part}$))
#:endfor
#:enddef

#:def GPU_USE_IS_PRESENT()
use openacc, only : acc_is_present
#:enddef

#:def GPU_UPDATE_DEVICE(vars)
!$acc update device(${vars}$)
#:enddef

#:def GPU_UPDATE_HOST(vars)
!$acc update host(${vars}$)
#:enddef

#:def GPU_HOST_DATA_USE_DEVICE(vars)
!$acc host_data use_device(${vars}$)
#:enddef

#:def GPU_END_HOST_DATA()
!$acc end host_data
#:enddef

#:elif defined('OPENMP')

#:def GPU_DECLARE_CREATE(vars)
!$omp declare target(${vars}$)
#:enddef

#:def GPU_DECLARE_COPYIN(vars)
!$omp declare target(${vars}$)
#:enddef

#:def GPU_ROUTINE_SEQ()
!$omp declare target
#:enddef

#:def GPU_ROUTINE_VECTOR()
!$omp declare target
#:enddef

#:def GPU_ROUTINE()
!$omp declare target
#:enddef

#:def GPU_PARALLEL_LOOP(clauses='')
#! NVHPC does not support defaultmap with present
#! (amd)flang does not support multiple defaultmap clauses on a single directive,
#! there we only insert the most important one: allocatables should be present
#! Empty line on purpose

#if defined(__NVCOMPILER)
!$omp target teams distribute parallel do ${clauses}$
#elif defined(__flang__)
!$omp target teams distribute parallel do ${clauses}$ defaultmap(present:allocatable)
#else
!$omp target teams distribute parallel do ${clauses}$ defaultmap(alloc:aggregate) defaultmap(present:allocatable) defaultmap(present:pointer)
#endif
#:enddef

#:def GPU_PARALLEL_LOOP_GANG(clauses='')
#! see comments on GPU_PARALLEL_LOOP
#! NVHPC compiles teams distribute with a nested parallel do in generic mode
#! (inner region outlined, ~255 registers per thread), so use loop there;
#! see GPU_LOOP_VECTOR. Not for other compilers: flang and CCE lower a
#! nested loop bind(parallel) to a worksharing do with no enclosing
#! parallel, which leaves one thread per team.
#! Empty line on purpose

#if defined(__NVCOMPILER)
!$omp target teams loop ${clauses}$
#elif defined(__flang__)
!$omp target teams distribute ${clauses}$ defaultmap(present:allocatable)
#else
!$omp target teams distribute ${clauses}$ defaultmap(alloc:aggregate) defaultmap(present:allocatable) defaultmap(present:pointer)
#endif
#:enddef

#:def GPU_LOOP_VECTOR(clauses='')
#! bind(parallel) is needed: without it NVHPC's dependence analysis can
#! serialise the loop onto one thread per team
#! Empty line on purpose

#if defined(__NVCOMPILER)
!$omp loop bind(parallel) ${clauses}$
#else
!$omp parallel do ${clauses}$
#endif
#:enddef

#:def GPU_LOOP_SEQ()
#:enddef

#:def GPU_INDEPENDENT()
#:enddef

#:def GPU_ENTER_DATA_COPYIN(vars)
!$omp target enter data map(to: ${vars}$)
#:enddef

#:def GPU_ENTER_DATA_CREATE(vars)
!$omp target enter data map(alloc: ${vars}$)
#:enddef

#:def GPU_EXIT_DATA_DELETE(vars)
!$omp target exit data map(delete: ${vars}$)
#:enddef

#:def GPU_EXIT_DATA_DELETE_IF_PRESENT(vars)
#:set depth = 0
#:set parts = ['']
#:for ch in vars
#:if ch == '('
#:set depth = depth + 1
#:endif
#:if ch == ')'
#:set depth = depth - 1
#:endif
#:if ch == ',' and depth == 0
#:set parts = parts + ['']
#:else
#:set parts = parts[:-1] + [parts[-1] + ch]
#:endif
#:endfor
#:set parts = [p.strip() for p in parts]
#:for part in parts
if (omp_target_is_present(c_loc(${part}$), omp_get_default_device()) /= 0) then
    !$omp target exit data map(delete: ${part}$)
end if
#:endfor
#:enddef

#:def GPU_USE_IS_PRESENT()
use omp_lib, only: omp_target_is_present, omp_get_default_device
use, intrinsic :: iso_c_binding, only: c_loc
#:enddef

#:def GPU_UPDATE_DEVICE(vars)
#! Ensure each referenced variable has device storage before refreshing its
#! value. Entering the *base* variable (not the specific, possibly
#! subscripted, expression) is essential: repeatedly calling `enter data`
#! with a different array-element expression each time (e.g. rnode(1,igrid)
#! for varying igrid) confuses nvfortran's OpenMP runtime and silently drops
#! the data motion for every element after the first. Entering the same base
#! variable repeatedly is a harmless no-op once it is present, and the
#! `update to` below (which always transfers, regardless of presence) then
#! correctly refreshes the specific element.
#:set depth = 0
#:set parts = ['']
#:for ch in vars
#:if ch == '('
#:set depth = depth + 1
#:endif
#:if ch == ')'
#:set depth = depth - 1
#:endif
#:if ch == ',' and depth == 0
#:set parts = parts + ['']
#:else
#:set parts = parts[:-1] + [parts[-1] + ch]
#:endif
#:endfor
#:set parts = [p.strip() for p in parts]
#:set bases = list(dict.fromkeys([p.split('(')[0].strip() for p in parts]))
!$omp target enter data map(to: ${', '.join(bases)}$)
!$omp target update to(${vars}$)
#:enddef

#:def GPU_UPDATE_HOST(vars)
!$omp target update from(${vars}$)
#:enddef

#:def GPU_HOST_DATA_USE_DEVICE(vars)
!$omp target data use_device_addr(${vars}$)
#:enddef

#:def GPU_END_HOST_DATA()
!$omp end target data
#:enddef

#:else
#! CPU build

#:def GPU_DECLARE_CREATE(vars)
#:enddef

#:def GPU_DECLARE_COPYIN(vars)
#:enddef

#:def GPU_ROUTINE_SEQ()
#:enddef

#:def GPU_ROUTINE_VECTOR()
#:enddef

#:def GPU_ROUTINE()
#:enddef

#:def GPU_PARALLEL_LOOP(clauses='')
#:enddef

#:def GPU_PARALLEL_LOOP_GANG(clauses='')
#:enddef

#:def GPU_LOOP_VECTOR(clauses='')
#:enddef

#:def GPU_LOOP_SEQ()
#:enddef

#:def GPU_INDEPENDENT()
#:enddef

#:def GPU_ENTER_DATA_COPYIN(vars)
#:enddef

#:def GPU_ENTER_DATA_CREATE(vars)
#:enddef

#:def GPU_EXIT_DATA_DELETE(vars)
#:enddef

#:def GPU_EXIT_DATA_DELETE_IF_PRESENT(vars)
#:enddef

#:def GPU_USE_IS_PRESENT()
#:enddef

#:def GPU_UPDATE_DEVICE(vars)
#:enddef

#:def GPU_UPDATE_HOST(vars)
#:enddef

#:def GPU_HOST_DATA_USE_DEVICE(vars)
#:enddef

#:def GPU_END_HOST_DATA()
#:enddef

#:endif
