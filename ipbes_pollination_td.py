"""
This is "3b" from https://docs.google.com/document/d/1k5yyhisemNrjG7ZSFGfD7GbIGY3rcNWmNulL2Hcqr9Q/edit

Multiply by annual nutrient demand by age group and add up to get total demand
per grid cell for en, fe, fo, va at a degree scale raster.
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

N_WORKERS = 4
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

DIETARY_REQUIREMENTS_TABLE_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes-data', 'Dietary requirements RNI EAR_annual.csv')

DEGREE_GRID_PATH = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes stuff',
    'summary table shapefile', 'degree_basedata', 'grid_1_degree.shp')

WORKSPACE_DIR = 'ipbes_dietary_workspace'


def zonal_stats_to_pickle(
        raster_path, shapefile_path, field_id, target_pickle_file):
    """Calc zonal stats over raster w/ shapefile & save to pickle file."""
    result = pygeoprocessing.zonal_statistics(
        (raster_path, 1), shapefile_path,
        field_id, polygons_might_overlap=False)

    with open(target_pickle_file, 'wb') as pickle_file:
        pickle.dump(result, pickle_file)


def create_nutrient_demand_raster(
        nutrient, dietary_table, pickle_path_list, target_raster_path):
    """Add pickled degree stats and mult by nut to raster.

    Parameters:
        nutrient (string): nutrient to calculate
        dietary_table (pandas.DataFrame): dietary table with 'age' and
            `nutrient` headers.
        pickle_path_list (list): list of pickled zonal stats dictionaries of
            the form 'scenario_ageid'
        target_raster_path(string): path to target raster"""
    nodata = 0
    base_array = numpy.empty((180, 360), dtype=numpy.float32)
    base_array[:] = nodata
    wgs84_sr = osr.SpatialReference()
    wgs84_sr.ImportFromEPSG(4326)
    driver = gdal.GetDriverByName('GTiff')
    print 'create raster'
    summary_raster = driver.Create(
        target_raster_path, 360, 180, 1, gdal.GDT_Float32)
    summary_raster.SetProjection(wgs84_sr.ExportToWkt())
    wgs84_gt = [-180., 1., 0, 90., 0., -1.]
    summary_raster.SetGeoTransform(wgs84_gt)
    summary_band = summary_raster.GetRasterBand(1)
    summary_band.SetNoDataValue(nodata)
    summary_band.Fill(nodata)
    for pickle_path in pickle_path_list:
        age_id = re.match('.*_(.*)', pickle_path).group(1)

        nutrient_demand_raw = dietary_table.loc[
            dietary_table['age'] == age_id][nutrient]
        print age_id, nutrient, nutrient_demand_raw
        nutrient_demand = float(nutrient_demand_raw)

        with open(pickle_path, 'r') as pickle_file:
            grid_stats = pickle.load(pickle_file)

        for grid_id in grid_stats:
            grid_x = (grid_id - 1) % 360
            grid_y = (grid_id - 1) // 360

            # value = sum(table[nut, age] * count_age, for all ages)
            value = grid_stats[grid_id]['sum'] * nutrient_demand
            base_array[grid_y, grid_x] += value

    summary_band.WriteArray(base_array)


def main():
    """Entry point."""
    try:
        os.makedirs(WORKSPACE_DIR)
    except OSError:
        pass

    task_graph = taskgraph.TaskGraph(
        os.path.join(WORKSPACE_DIR, 'taskgraph_cache'), N_WORKERS,
        reporting_interval=REPORTING_INTERVAL, dry_run=DRY_RUN)

    raster_path_list = glob.glob(os.path.join(
        BASE_DROPBOX_DIR, 'ipbes stuff',
        'ipbes_pollination_1_27_2018', 'gwpop', '*.tif'))

    dietary_table = pandas.read_csv(DIETARY_REQUIREMENTS_TABLE_PATH)

    pickle_task_list = []
    for raster_path in raster_path_list:
        pickle_path = os.path.join(
            WORKSPACE_DIR, os.path.splitext(os.path.basename(raster_path))[0])
        pickle_task = task_graph.add_task(
            func=zonal_stats_to_pickle,
            args=(raster_path, DEGREE_GRID_PATH, 'GRIDCODE', pickle_path),
            target_path_list=[pickle_path])
        pickle_task_list.append((pickle_path, pickle_task))

    for scenario in ['cur', 'ssp1', 'ssp3', 'ssp5']:
        for nutrient in  ['en', 'fe', 'fo', 'va']:
            print 'create 1 degree raster %s %s' % (scenario, nutrient)

            scenario_task_list = [
                x for x in pickle_task_list
                if scenario in x[0] and
                'rural' not in x[0] and
                'total' not in x[0]]
            td_raster_path = os.path.join(
                WORKSPACE_DIR, 'td_%s_%s.tif' % (nutrient, scenario))

            task_graph.add_task(
                func=create_nutrient_demand_raster,
                args=(
                    nutrient, dietary_table, [x[0] for x in scenario_task_list],
                    td_raster_path),
                target_path_list=[td_raster_path],
                dependent_task_list=[x[1] for x in scenario_task_list])

    #[cur | ssp[1 | 3 | 5]]_gpwpop_[014 | 1564 | 65p][f | m] * table[en | fe | fo | va ; 0-14 F | 0 -14 M |15-64 F  | 15-64 M | 65+ F | 65+ M]

    task_graph.close()
    task_graph.join()

if __name__ == '__main__':
    main()
