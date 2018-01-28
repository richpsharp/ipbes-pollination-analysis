"""
This is "3b" from https://docs.google.com/document/d/1k5yyhisemNrjG7ZSFGfD7GbIGY3rcNWmNulL2Hcqr9Q/edit

Multiply by annual nutrient demand by age group and add up to get total demand
per grid cell for en, fe, fo, va at a degree scale raster.
"""
import logging
import glob
import os

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('ipbes-cv')
LOGGER.setLevel(logging.DEBUG)

POSSIBLE_DROPBOX_LOCATIONS = [
    r'D:\Dropbox',
    r'C:\Users\Rich\Dropbox',
    r'C:\Users\rpsharp\Dropbox',
    r'E:\Dropbox']

LOGGER.info("checking dropbox locations")
for path in POSSIBLE_DROPBOX_LOCATIONS:
    print path
    if os.path.exists(path):
        BASE_DROPBOX_DIR = path
        break
LOGGER.info("found %s", BASE_DROPBOX_DIR)


def main():
    """Entry point."""
    raster_path_list = glob.glob(os.path.join(
        BASE_DROPBOX_DIR, 'rps_bck_shared_stuff', 'ipbes stuff',
        'ipbes_pollination_1_27_2018', 'gwpop', '*.tif'))
    print raster_path_list
     #[cur | ssp[1 | 3 | 5]]_gpwpop_[014 | 1564 | 65p][f | m] * table[en | fe | fo | va ; 0-14 F | 0 -14 M |15-64 F  | 15-64 M | 65+ F | 65+ M]


if __name__ == '__main__':
    main()
