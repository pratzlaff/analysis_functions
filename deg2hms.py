import argparse
import astropy.units as u
from astropy.coordinates import SkyCoord

def deg2hms(ra, dec):
    target = SkyCoord(ra=ra, dec=dec, unit=u.deg)
    print(f'{int(target.ra.hms.h):02d}:{int(target.ra.hms.m):02d}:{target.ra.hms.s:0.4f}',
          ' '
          f'{int(target.dec.dms.d):02d}:{int(target.dec.dms.m):02d}:{target.dec.dms.s:07.4f}'
          )

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
