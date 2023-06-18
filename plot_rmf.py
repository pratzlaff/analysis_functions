import argparse
import matplotlib.pyplot as plt
import re

import response

def main():
    parser = argparse.ArgumentParser(
        description='Plot RMFs.'
    )
    parser.add_argument('-o', '--outfile', help='Save plot to named file.')
    parser.add_argument('files', nargs='+', help='RMF file(s).')
    args = parser.parse_args()

    for i in range(len(args.files)):
        rmf = args.files[i]
        label = re.search('([^/]+$)', rmf).group(1)
        bin_lo, bin_hi, resp = response.read_rmf(rmf)
        plt.plot(0.5*(bin_lo+bin_hi), resp, label=label)

    plt.xlabel("λ (Å)")

    if len(args.files)>1:
        plt.legend()
    else:
        plt.title(label)

    if args.outfile:
        plt.savefig(args.outfile)
    else:
        plt.show()
    
if __name__ == '__main__':
    main()
