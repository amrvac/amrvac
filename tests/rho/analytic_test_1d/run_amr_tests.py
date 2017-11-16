#!/usr/bin/env python

from subprocess import check_output
import re
import sys

limiters = ['minmod' , 'woodward', 'cada3' ]

flux_schemes = ['hll']

time_integrators = ['twostep', 'threestep', 'fourstep', 'ssprk54']

refine_max_levels = [ 1, 2, 3, 4, 5]

all_results = []

n_results = len(flux_schemes) * len(time_integrators) * len(limiters) * len(refine_max_levels)
n_done = 0

def progress_bar(pct):
    n = int(pct * 0.4)
    sys.stdout.write('\r')
    sys.stdout.write('[{0:40s}] {1:.1f}%'.format('=' * n, pct))
    sys.stdout.flush()

for fs in flux_schemes:
    for ti in time_integrators:
        for lim in limiters:
            for lev in refine_max_levels:
                with open('this_scheme.par', 'w') as f:
                    f.write('! Temporary par file\n')
                    f.write('&methodlist\n')
                    f.write(' time_integrator = "{}"\n'.format(ti))
                    f.write(' flux_scheme = "{}"\n'.format(fs))
                    f.write(' limiter = "{}"\n'.format(lim))
                    f.write('/\n')
                    f.write('&meshlist\n')
                    f.write(' refine_max_level = {}\n'.format(lev))
                    f.write('/\n')

                res = check_output(["./amrvac", "-i", "amrvac.par", "this_scheme.par"]).decode(encoding="utf-8")
                mobj = re.search(r'time -- RMSE:  0.10000000E\+01\s+(\S*)', res, re.MULTILINE)
                error = float(mobj.group(1))
                this_result = [fs, ti, lim, lev, error]
                all_results.append(this_result)
                n_done += 1
                progress_bar((100.0 * n_done) / n_results)

sorted_results = sorted(all_results, key=lambda x: x[4])

print("\n\nSorted results:\n")

best_err = sorted_results[0][4]
for i, res in enumerate(sorted_results):
    print("{:5d} {:12.2e} (+{:5.1f}%): {:8s} {:20s} {:12s} {:5d}".format(
        i+1, res[4], 100*(res[4]-best_err)/best_err,
        res[0], res[1], res[2], res[3]))
