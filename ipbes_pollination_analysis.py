"""Pollination ipbes"""
import os
import logging

from osgeo import gdal
import taskgraph
import pygeoprocessing
import numpy


logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

WORKSPACE_DIR = 'pollination_workspace'

def add_arrays(*array_list):
    return numpy.sum(array_list, axis=0)

def main():
    """Entry point."""

    if not os.path.exists(WORKSPACE_DIR):
        os.makedirs(WORKSPACE_DIR)

    crop_path_list = [
        r"C:\Users\Rich\Dropbox\ipbes-pollination-analysis\LUH2_1KM\c3ann.flt",
        r"C:\Users\Rich\Dropbox\ipbes-pollination-analysis\LUH2_1KM\c3nfx.flt",
        r"C:\Users\Rich\Dropbox\ipbes-pollination-analysis\LUH2_1KM\c3per.flt",
        r"C:\Users\Rich\Dropbox\ipbes-pollination-analysis\LUH2_1KM\c4ann.flt",
        r"C:\Users\Rich\Dropbox\ipbes-pollination-analysis\LUH2_1KM\c4per.flt",
    ]

    forest_path_list = [
        r"C:\Users\Rich\Dropbox\ipbes-pollination-analysis\LUH2_1KM\primf.flt",
        r"C:\Users\Rich\Dropbox\ipbes-pollination-analysis\LUH2_1KM\primn.flt",
        r"C:\Users\Rich\Dropbox\ipbes-pollination-analysis\LUH2_1KM\secdf.flt",
        r"C:\Users\Rich\Dropbox\ipbes-pollination-analysis\LUH2_1KM\secdn.flt",
    ]

    crop_class_path = os.path.join(WORKSPACE_DIR, 'crop_classes.tif')
    forest_class_path = os.path.join(WORKSPACE_DIR, 'forest_classes.tif')

    task_graph = taskgraph.TaskGraph(
        os.path.join(WORKSPACE_DIR, 'taskgraph_cache'), 4)

    for raster_path_list, target_path in [
            (crop_path_list, crop_class_path),
            (forest_path_list, forest_class_path)]:
        task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(x, 1) for x in raster_path_list], add_arrays, target_path,
                gdal.GDT_Float32, -9999),
            target_path_list=[target_path]
            )
    task_graph.join()



if __name__ == '__main__':
    main()
