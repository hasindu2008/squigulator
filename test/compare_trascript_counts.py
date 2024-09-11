#! /usr/bin/env python3

import os
import sys
import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['figure.dpi'] = 80
rcParams['figure.figsize'] = [20.0, 12.0]
# matplotlib.rcParams['savefig.dpi'] = 300
import numpy as np


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)



def main():
    '''
    main function
    '''

    parser = MyParser(
        description="compare transcript counts between two files")
    #group = parser.add_mutually_exclusive_group()
    parser.add_argument("--c1",
                        help="counts 1 file")
    parser.add_argument("--c2",
                        help="counts 2 file")
    parser.add_argument("-n", "--num_col", type=int, default="0",
                        help="column number (0-indexed) that has the counts")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="plot c1 vs c2")

    args = parser.parse_args()

    # print help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    c1_file = args.c1.split(".")
    # doing [-1:] rather than [-1] gives a list back
    name, ext = [".".join(c1_file[:-1])], c1_file[-1:]
    c1_name = ".".join(name + ["pass"] + ext)
    c2_file = args.c2.split(".")
    name, ext = [".".join(c2_file[:-1])], c2_file[-1:]
    c2_name = ".".join(name + ["pass"] + ext)


    c1 = {}
    with open(args.c1, 'r') as f:
        for line in f:
            l = line.split()
            if args.num_col == 0:
                count = l[0]
                name = l[1]
            else:
                count = l[1]
                name = l[0]
            c1[name] = count

    c2 = {}
    with open(args.c2, 'r') as f:
        for line in f:
            l = line.split()
            if args.num_col == 0:
                count = l[0]
                name = l[1]
            else:
                count = l[1]
                name = l[0]
            c2[name] = count

    both = {}
    for name in c1:
        if name in c2:
            both[name] = [c1[name], c2[name]]
        else:
            both[name] = [c1[name], 0]

    for name in c2:
        if name not in both:
            both[name] = [0, c2[name]]

    for name in both:
        print("{}\t{}\t{}".format(name, both[name][0], both[name][1]))

    if args.plot:
        # prep data
        x = []
        y = []
        for name in both:
            x.append(int(both[name][0]))
            y.append(int(both[name][1]))

        x_arr = np.array(x)
        y_arr = np.array(y)
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        fig.suptitle("{} vs {}".format(c1_name, c2_name), fontsize=16)
        #ax.axline((0, 0), slope=1)
        ax.set_xscale('log', base=2)
        ax.set_yscale('log', base=2)
        plt.scatter(x_arr, y_arr)
        plt.show()
        plt.savefig("plot.png")
        plt.clf()


if __name__ == '__main__':
    main()
