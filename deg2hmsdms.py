import argparse
import astropy.units as u
from astropy.coordinates import SkyCoord

def deg2hms(ra, dec):
    target = SkyCoord(ra=ra, dec=dec, unit=u.deg)
    print(target.to_string('hmsdms', sep=':'))

def main():
    parser = argparse.ArgumentParser(
        description='Convert ra/dec from degrees to [HD]MS'
    )
    parser.add_argument('ra', type=float)
    parser.add_argument('dec', type=float)
    args = parser.parse_args()
    deg2hms(args.ra, args.dec)

if __name__ == '__main__':
    main()
