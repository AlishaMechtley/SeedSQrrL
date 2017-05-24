from __future__ import division, print_function
import matplotlib.pyplot as pp
import numpy as np
from StringIO import StringIO


def parse_resource_log(filename):
    # preparse to remove header lines
    with open(filename) as fh:
        lines = fh.readlines()
    header = lines[0]
    lines = [line for line in lines if 'PID' not in line]
    lines.insert(0, header)
    contents = '\n'.join(lines)
    contents = StringIO(contents)

    arr = np.genfromtxt(contents, dtype=None, names=True, usecols=(0, 2, 3))
    # Every row matching the first PID is parent python process, shows up
    # first in each block. Sum resource stats over successive blocks.
    parent_rows = np.where(arr['PID'] == arr[0]['PID'])
    cpu_arr = np.add.reduceat(arr['C'], parent_rows[0])
    mem_arr = np.add.reduceat(arr['RSS'], parent_rows[0])
    return cpu_arr, mem_arr


def cpu_histogram(cpu_arr, ncpus=8, clip=True):
    # ps has some funky sampling such that it can report over 100% cpu usage.
    # clip values to 100%
    cpu_unnorm = cpu_arr * ncpus
    if clip:
        cpu_unnorm = np.clip(cpu_unnorm, 0, 100 * ncpus)
    pp.hist(cpu_unnorm, color='RoyalBlue')
    stats = 'median: {:0.0f}%\nmax: {:0.0f}%'.format(
        np.median(cpu_unnorm), np.max(cpu_unnorm))
    pp.text(0.95, 0.95, stats, transform=pp.gca().transAxes,
            fontsize=16, va='top', ha='right')
    pp.xlabel('CPU Usage (1 processor = 100%)')
    pp.show()
    return


def memory_histogram(mem_arr):
    # convert from kB to GB
    mem_gb = mem_arr.astype('float64') / 1024 ** 2
    pp.hist(mem_gb, color='RoyalBlue')
    stats = 'median: {:0.1f}GB\nmax: {:0.1f}GB'.format(
        np.median(mem_gb), np.max(mem_gb))
    pp.text(0.95, 0.95, stats, transform=pp.gca().transAxes,
            fontsize=16, va='top', ha='right')
    pp.xlabel('Memory Usage (GB)')
    pp.show()
    return


def make_histograms():
    import sys
    filename = sys.argv[1]
    cpu_arr, mem_arr = parse_resource_log(filename)

    cpu_histogram(cpu_arr)
    memory_histogram(mem_arr)


if __name__ == '__main__':
    make_histograms()
