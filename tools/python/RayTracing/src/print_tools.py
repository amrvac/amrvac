"""
Module containing some methods to facilitate printing to console.

@author: Niels Claes
"""

import sys

def trim_filename(filename):
    """
    Returns name of snapshot, used to save and load interpolated data.
    Assumes filenames ending in MPI-AMRVAC format 0000.dat
    """
    try:
        filenumber = int(filename[-8:-4])
        for c in range(len(filename)-1, -1, -1):
            if filename[c] == '/':
                return filename[c+1:-4]
        return filename[::-4]
    except ValueError:
        return ""
    
def progress(count, total, status=''):
    bar_len = 50
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '#' * filled_len + ' ' * (bar_len - filled_len)

    sys.stdout.write('\r    [%s] %s%s %s' % (bar, percents, '%', status))
    sys.stdout.flush()

