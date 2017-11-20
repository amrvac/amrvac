#!/usr/bin/env python

from subprocess import check_output
import re
import numpy as np
import sys
import argparse

def get_args():
    # Get and parse the command line arguments
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='''Tool to run convergence tests''')
    p.add_argument('-limiters', type=str, nargs='+', default=['koren'],
                        help='Which limiter(s) to include')
    p.add_argument('-schemes', type=str, nargs='+', default=['hll'],
                        help='Which flux scheme(s) to include')
    p.add_argument('-integrators', type=str, nargs='+', default=['threestep'],
                        help='Which time integrators(s) to include')

    p.add_argument('-nx', type=int, default=1024,
                        help='Maximum number of grid points')
    p.add_argument('-nxlow', type=int, default=16,
                        help='Minimum number of grid points')
    return p.parse_args()


def progress_bar(pct):
    n = int(pct * 0.4)
    sys.stdout.write('\r')
    sys.stdout.write('[{0:40s}] {1:.1f}%'.format('=' * n, pct))
    sys.stdout.flush()


if __name__ == '__main__':
    p = get_args()

    all_results = []
    ngrids = int(np.log2(float(p.nx)/p.nxlow) + 1.5)
    nx_list = p.nxlow * np.power(2 * np.ones(ngrids), np.arange(ngrids))
    nx_list = nx_list.astype(int)
    n_results = len(p.schemes) * len(p.integrators) * len(p.limiters) * len(nx_list)
    n_done = 0

    for fs in p.schemes:
        for ti in p.integrators:
            for lim in p.limiters:
                nx_results = np.zeros(ngrids,)

                for i, nx in enumerate(nx_list):
                    with open('this_scheme.par', 'w') as f:
                       f.write('! Temporary par file\n')
                       f.write('&methodlist\n')
                       f.write(' time_integrator = "{}"\n'.format(ti))
                       f.write(' flux_scheme = "{}"\n'.format(fs))
                       f.write(' limiter = "{}"\n'.format(lim))
                       f.write('/\n')
                       f.write('&meshlist\n')
                       f.write(' domain_nx1 = {}\n'.format(nx))
                       f.write('/\n')

                    res = check_output(["./amrvac", "-i", "amrvac.par", "this_scheme.par"]).decode(encoding="utf-8")
                    mobj = re.search(r'time -- RMSE:  0.10000000E\+01\s+(\S*)', res, re.MULTILINE)
                    error = float(mobj.group(1))
                    nx_results[i] = error
                    n_done += 1
                    progress_bar((100.0 * n_done) / n_results)

                all_results.append([fs, ti, lim, nx_results])


    print("\nResults for nx = {}\n".format(nx_list))
    for res in all_results:
        print("{:30s}: {}".format('-'.join(res[0:3]), res[3]))
