"""
This is "1" from https://docs.google.com/document/d/1k5yyhisemNrjG7ZSFGfD7GbIGY3rcNWmNulL2Hcqr9Q/edit
"""
import re
import pickle
import logging
import glob
import os

import numpy
from osgeo import gdal
from osgeo import osr
import taskgraph
import pandas
import pygeoprocessing

N_WORKERS = -1
REPORTING_INTERVAL = 30.0
DRY_RUN = False

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

WORKSPACE_DIR = 'cont_rasters_workspace'

POLL_SERV_DIR = glob.glob(os.path.join(
    BASE_DROPBOX_DIR, 'ipbes stuff', 'ipbes_pollination_1_27_2018',
    'raster_results', 'poll_serv*.tif'))

TD_DIR = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes stuff', 'ipbes_pollination_1_27_2018',
    'td_results')

#"D:\Dropbox\ipbes stuff\ipbes_pollination_1_27_2018\raster_results\poll_serv_cur_vita.tif"

long_to_short = {
    'energy': 'en',
    'vita': 'va',
    'fe': 'fe',
    'folate': 'fo'
}


def div_rasters(num, denom):
    result = numpy.empty_like(num)
    result[:] = 0
    valid_mask = denom != 0.0
    result[valid_mask] = num[valid_mask] / denom[valid_mask]
    return result


def change_rasters(num, denom):
    result = numpy.empty_like(num)
    result[:] = -1
    valid_mask = denom != 0.0
    result[valid_mask] = (
        num[valid_mask] - denom[valid_mask]) / denom[valid_mask]
    return result

def main():
    """Entry point."""
    try:
        os.makedirs(WORKSPACE_DIR)
    except OSError:
        pass

    task_graph = taskgraph.TaskGraph(
        os.path.join(WORKSPACE_DIR, 'taskgraph_cache'), N_WORKERS,
        reporting_interval=REPORTING_INTERVAL, dry_run=DRY_RUN)

    for poll_serv_path in POLL_SERV_DIR:
        scenario, nut = re.match('.*_(.*)_(.*)\.tif', poll_serv_path).groups()
        print poll_serv_path, nut, scenario
        td_path = os.path.join(TD_DIR, 'td_%s_%s.tif' % (
            long_to_short[nut], scenario))
        print td_path, os.path.exists(td_path)

        cont_d_path = os.path.join(
            WORKSPACE_DIR, 'cont_d_%s_%s.tif' % (
                long_to_short[nut], scenario))

        pygeoprocessing.raster_calculator(
            [(poll_serv_path, 1), (td_path, 1)], div_rasters, cont_d_path,
            gdal.GDT_Float32, 0)

    for scenario in ['ssp1', 'ssp3', 'ssp5']:
        for nut in ['en', 'fo', 'fe', 'va']:
            c_cont_d_path = os.path.join(
                WORKSPACE_DIR, 'c_cont_d_%s_%s.tif' % (nut, scenario))
            cont_d_ssp = os.path.join(
                WORKSPACE_DIR, 'cont_d_%s_%s.tif' % (nut, scenario))
            cont_d_cur = os.path.join(
                WORKSPACE_DIR, 'cont_d_%s_cur.tif' % (nut,))

            pygeoprocessing.raster_calculator(
                [(cont_d_ssp, 1), (cont_d_cur, 1)], change_rasters,
                c_cont_d_path, gdal.GDT_Float32, -1)

    task_graph.close()
    task_graph.join()



if __name__ == '__main__':
    main()
