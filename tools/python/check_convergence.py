#!/usr/bin/env python

from subprocess import check_output, CalledProcessError
import re
import numpy as np
import sys
import argparse
import itertools

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
    p.add_argument('-var', type=str, default='err_1',
                   help='Which variable to consider')
    p.add_argument('-limiters', type=str, nargs='+', default=['koren'],
                   help='Which limiter(s) to include')
    p.add_argument('-schemes', type=str, nargs='+', default=['hll'],
                   help='Which flux scheme(s) to include')
    p.add_argument('-integrators', type=str, nargs='+', default=['threestep'],
                   help='Which time integrators(s) to include')
    p.add_argument('-powers', type=float, nargs='+', default=[],
                   help='Show these powers in the graph')
    p.add_argument('-iprob', type=int, default=1,
                   help='Which iprob value to use')
    p.add_argument('-dim', type=int, default=1,
                   help='Test problem dimension')
    p.add_argument('-nx', type=int, default=1024,
                   help='Maximum number of grid points')
    p.add_argument('-nxlow', type=int, default=32,
                   help='Minimum number of grid points')
    p.add_argument('-nxblock', type=int, default=16,
                   help='Size of grid blocks')
    p.add_argument('-fig', type=str,
                   help='Store a figure with this name (e.g., test.pdf)')
    p.add_argument('-np', type=int, default=4,
                   help='Number of CPUs to use (with mpirun)')
    return p.parse_args()


def progress_bar(pct):
    n = int(pct * 0.4)
    sys.stdout.write('\r')
    sys.stdout.write('[{0:40s}] {1:.1f}%'.format('=' * n, pct))
    sys.stdout.flush()


if __name__ == '__main__':
    p = get_args()

    all_results = []

    if p.nxlow % p.nxblock != 0:
        print("nxblock={} does not divide nxlow={}".format(
            p.nxblock, p.nxlow))
        sys.exit(1)

    if p.nx < p.nxlow:
        print("Error: nx={} smaller than nxlow={}".format(p.nx, p.nxlow))
        sys.exit(1)

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
                    with open('CONV_TEST.par', 'w') as f:
                       f.write('! Temporary par file\n')
                       f.write('&methodlist\n')
                       f.write(' time_integrator = "{}"\n'.format(ti))
                       f.write(' flux_scheme = "{}"\n'.format(fs))
                       f.write(' limiter = "{}"\n'.format(lim))
                       f.write('/\n')
                       f.write('&meshlist\n')
                       f.write(' domain_nx1 = {}\n'.format(nx))
                       f.write(' block_nx1 = {}\n'.format(p.nxblock))
                       if p.dim > 1:
                           f.write(' domain_nx2 = {}\n'.format(nx))
                           f.write(' block_nx2 = {}\n'.format(p.nxblock))
                       if p.dim > 2:
                           f.write(' domain_nx3 = {}\n'.format(nx))
                           f.write(' block_nx3 = {}\n'.format(p.nxblock))
                       f.write(' iprob = {}\n'.format(p.iprob))
                       f.write('/\n')

                    try:
                        res = check_output(["mpirun", "-np", str(p.np),
                                            "./amrvac", "-i", p.par,
                                            "CONV_TEST.par"])
                        res = res.decode(encoding="utf-8")
                    except CalledProcessError as e:
                        print("\nAMRVAC returned an error")
                        print(e.output)
                        sys.exit(1)

                    # Catch last line starting with CONVTEST
                    line = re.findall(r'^\s*CONVTEST.*$', res, re.MULTILINE)[-1]

                    # Format should be CONVTEST (var1,var2,...): val1 val2 val3
                    # Extract variable names and values
                    i0 = line.find('(')
                    i1 = line.find(')')
                    names = line[i0+1:i1].lstrip()
                    names = re.split(r'[ ,\t\n]+', names)
                    i0 = line.find(':')
                    vals = line[i0+1:].lstrip()
                    vals = re.split(r'[ ,\t\n]+', vals)
                    vals = map(float, vals)

                    try:
                        i0 = names.index(p.var)
                    except:
                        print("{} not in list of variables: {}".format(
                            p.var, ' '.join(names)))
                        sys.exit(1)

                    nx_results[i] = vals[i0]
                    n_done += 1
                    progress_bar((100.0 * n_done) / n_results)

                all_results.append([fs, ti, lim, nx_results])

    if p.fig:
        # Generate figure
        import matplotlib.pyplot as plt

        plt.figure(figsize=[10.,6.])
        marker = itertools.cycle(('+', '.', 'x', '*', 'v',
                                  '^', '<', '>', '8'))
        linestyle = itertools.cycle((':', '-.', '--', '-'))

        mean_val = 0.0
        for res in all_results:
            lbl = '-'.join(res[0:3])
            line, = plt.loglog(nx_list, res[3], linestyle=linestyle.next(),
                               marker=marker.next(), label=lbl)
            mean_val = mean_val + res[3][0] / len(all_results)

        for pp in p.powers:
            plt.loglog(nx_list, mean_val * (float(nx_list[0])/nx_list)**pp,
                       linestyle='--', label='order-{:.1f}'.format(pp))
        plt.ylabel('Error {}'.format(p.var))
        plt.xlabel('Grid points')
        leg = plt.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.5)
        plt.xlim([nx_list[0]*0.9, nx_list[-1]*1.1])
        plt.savefig(p.fig)
        print('\nSaved figure {}'.format(p.fig))
    else:
        # Output text to screen
        print("\n{} error for nx = {}\n".format(p.var, nx_list))
        for res in all_results:
            errs = ' '.join(['{:10.3e}'.format(err) for err in res[3]])
            print("{:s}: {}".format('-'.join(res[0:3]), errs))
