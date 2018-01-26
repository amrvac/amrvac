#!/usr/bin/env python

from subprocess import check_output, CalledProcessError
import re
import numpy as np
import sys
import argparse

def get_args():
    # Get and parse the command line arguments
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Tool to run convergence tests with MPI-AMRVAC, by increasing the grid
        resolution. Flux schemes, time integrators and limiters can be compared
        (simultaneously). Without the -fig argument, the program outputs the
        errors. The program requires a suitable par-file to be present.""",
        epilog='Example: ./check_convergence.py -limiters superbee koren -fig test.pdf')

    p.add_argument('-par', type=str, default='amrvac.par',
                   help='Which par files to use')

    p.add_argument('-limiters', type=str, nargs='+', default=['koren'],
                   help='Which limiter(s) to include')
    p.add_argument('-schemes', type=str, nargs='+', default=['hll'],
                   help='Which flux scheme(s) to include')
    p.add_argument('-integrators', type=str, nargs='+', default=['threestep'],
                   help='Which time integrators(s) to include')
    p.add_argument('-norm', type=str, default='2', choices=['1', '2', 'inf'],
                   help='Which norm to use (1, 2, inf)')
    p.add_argument('-powers', type=float, nargs='+', default=[],
                   help='Show these powers in the graph')
    p.add_argument('-iprob', type=int, default=1,
                   help='Which iprob value to use')
    p.add_argument('-dim', type=int, default=1,
                   help='Test problem dimension')
    p.add_argument('-nx', type=int, default=1024,
                   help='Maximum number of grid points')
    p.add_argument('-nxlow', type=int, default=16,
                   help='Minimum number of grid points')
    p.add_argument('-fig', type=str,
                   help='Store a figure with this name (e.g., test.pdf)')
    return p.parse_args()


def progress_bar(pct):
    n = int(pct * 0.4)
    sys.stdout.write('\r')
    sys.stdout.write('[{0:40s}] {1:.1f}%'.format('=' * n, pct))
    sys.stdout.flush()


if __name__ == '__main__':
    p = get_args()

    all_results = []

    # Number of grids to use
    ngrids = int(np.log2(float(p.nx)/p.nxlow) + 1.5)

    # Generate list with grid numbers (e.g., 16, 32, 64, ...)
    nx_list = p.nxlow * np.power(2 * np.ones(ngrids), np.arange(ngrids))
    nx_list = nx_list.astype(int)

    # Number of results to generate
    n_results = len(p.schemes) * len(p.integrators) * len(p.limiters) * ngrids
    n_done    = 0

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
                       f.write(' block_nx1 = {}\n'.format(p.nxlow))
                       if p.dim > 1:
                           f.write(' domain_nx2 = {}\n'.format(nx))
                           f.write(' block_nx2 = {}\n'.format(p.nxlow))
                       if p.dim > 2:
                           f.write(' domain_nx3 = {}\n'.format(nx))
                           f.write(' block_nx3 = {}\n'.format(p.nxlow))
                       f.write(' iprob = {}\n'.format(p.iprob))
                       f.write('/\n')

                    try:
                        res = check_output(["./amrvac", "-i", p.par,
                                            "this_scheme.par"])
                        res = res.decode(encoding="utf-8")
                    except CalledProcessError as e:
                        print("\nAMRVAC returned an error")
                        print(e.output)
                        sys.exit(1)

                    # Catch line before "Total timeloop took"
                    mobj = re.search(r'^(.*)$\n^\s*Total timeloop took',
                                     res, re.MULTILINE)

                    # Split the string which contains the errors
                    errors = mobj.group(1).split()

                    # Extract the appropriate norm
                    if p.norm == '1':
                        nx_results[i] = float(errors[-3])
                    elif p.norm == '2':
                        nx_results[i] = float(errors[-2])
                    else:
                        nx_results[i] = float(errors[-1])

                    n_done += 1
                    progress_bar((100.0 * n_done) / n_results)

                all_results.append([fs, ti, lim, nx_results])

    if p.fig:
        # Generate figure
        import matplotlib.pyplot as plt

        plt.figure(figsize=[10.,6.])
        for res in all_results:
            lbl = '-'.join(res[0:3])
            line, = plt.loglog(nx_list, res[3], linestyle='-',
                               marker='.', label=lbl)

        for pp in p.powers:
            plt.loglog(nx_list, (float(nx_list[0])/nx_list)**pp,
                       linestyle='--', label='order-{:.1f}'.format(pp))
        plt.ylabel('Error {}-norm'.format(p.norm))
        plt.xlabel('Grid points')
        plt.legend(loc='lower left')
        plt.xlim([nx_list[0]*0.9, nx_list[-1]*1.1])
        plt.savefig(p.fig)
        print('\nSaved figure {}'.format(p.fig))
    else:
        # Output text to screen
        print("\n{}-norm for nx = {}\n".format(p.norm, nx_list))
        for res in all_results:
            errs = ' '.join(['{:10.3e}'.format(err) for err in res[3]])
            print("{:s}: {}".format('-'.join(res[0:3]), errs))
