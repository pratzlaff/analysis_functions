import argparse
import matplotlib.pyplot as plt
import re

import response

def main():
    parser = argparse.ArgumentParser(
        description='Plot RMF file1 / file2.'
    )
    parser.add_argument('-o', '--outfile', help='Save plot to named file.')
    parser.add_argument('file1', help='RMF file 1.')
    parser.add_argument('file2', help='RMF file 2.')
    args = parser.parse_args()

    bin_lo1, bin_hi1, resp1 = response.read_rmf(args.file1)
    bin_lo2, bin_hi2, resp2 = response.read_rmf(args.file2)

    plt.plot(0.5*(bin_lo1+bin_hi1), resp1/resp2)

    plt.xlabel("λ (Å)")

    if args.outfile:
        plt.savefig(args.outfile)
    else:
        plt.show()
    
if __name__ == '__main__':
    main()
