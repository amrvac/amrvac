#!/usr/bin/env python

from subprocess import check_output
import re
import sys

limiters = ['minmod' , 'woodward', 'mcbeta', 'superbee', 'albada',
            'koren', 'vanleer', 'cada', 'cada3', 'mp5']

flux_schemes = ['hll']

time_integrators = ['twostep', 'twostep_trapezoidal', 'threestep',
                    'fourstep', 'ssprk54', 'ssprk43', 'rk4', 'jameson']

all_results = []

n_results = len(flux_schemes) * len(time_integrators) * len(limiters)
n_done = 0

def progress_bar(pct):
    n = int(pct * 0.4)
    sys.stdout.write('\r')
    sys.stdout.write('[{0:40s}] {1:.1f}%'.format('=' * n, pct))
    sys.stdout.flush()

for fs in flux_schemes:
    for ti in time_integrators:
        for lim in limiters:
            with open('this_scheme.par', 'w') as f:
                f.write('! Temporary par file\n')
                f.write('&methodlist\n')
                f.write(' time_integrator = "{}"\n'.format(ti))
                f.write(' flux_scheme = "{}"\n'.format(fs))
                f.write(' limiter = "{}"\n'.format(lim))
                f.write('/\n')

            res = check_output(["./amrvac", "-i", "amrvac.par", "this_scheme.par"]).decode(encoding="utf-8")
            mobj = re.search(r'time -- RMSE:  0.10000000E\+01\s+(\S*)', res, re.MULTILINE)
            error = float(mobj.group(1))
            this_result = [fs, ti, lim, error]
            all_results.append(this_result)
            n_done += 1
            progress_bar((100.0 * n_done) / n_results)

sorted_results = sorted(all_results, key=lambda x: x[3])

print("\n\nSorted results:\n")

best_err = sorted_results[0][3]
for i, res in enumerate(sorted_results):
    print("{:5d} {:12.2e} (+{:5.1f}%): {:8s} {:20s} {:12s}".format(
        i+1, res[3], 100*(res[3]-best_err)/best_err,
        res[0], res[1], res[2]))
