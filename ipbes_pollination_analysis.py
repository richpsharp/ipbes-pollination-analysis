"""This pollination analysis is guided by this doc https://docs.google.com/document/d/1k5yyhisemNrjG7ZSFGfD7GbIGY3rcNWmNulL2Hcqr9Q/edit"""
import os
import logging

from osgeo import osr
from osgeo import gdal
import taskgraph
import pygeoprocessing
import numpy


logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

BASE_DATA_DIR = r"C:\Users\Rich\Dropbox\ipbes-pollination-analysis"
WORKSPACE_DIR = 'pollination_workspace'
GTIFF_CREATION_OPTIONS = ('TILED=YES', 'BIGTIFF=IF_SAFER', 'COMPRESS=DEFLATE')


NODATA = -9999

# any proportion of ag above this is ag
CROP_THRESHOLD_VALUE = 0.5


def threshold_crop_raster(array):
    result = numpy.zeros(array.shape, dtype=numpy.int8)
    result[:] = 2
    valid_mask = array != NODATA
    result[valid_mask] = array[valid_mask] > CROP_THRESHOLD_VALUE
    return result


def mult_arrays(*array_list):
    stack = numpy.stack(array_list)
    valid_mask = (
        numpy.bitwise_and.reduce(stack != NODATA, axis=0))
    n_valid = numpy.count_nonzero(valid_mask)
    broadcast_valid_mask = numpy.broadcast_to(valid_mask, stack.shape)
    valid_stack = stack[broadcast_valid_mask].reshape(
        len(array_list), n_valid)
    result = numpy.empty(array_list[0].shape, dtype=numpy.float32)
    result[:] = NODATA
    result[valid_mask] = numpy.prod(valid_stack, axis=0)
    return result


def add_arrays(*array_list):
    stack = numpy.stack(array_list)
    valid_mask = (
        numpy.bitwise_and.reduce(stack != NODATA, axis=0))
    n_valid = numpy.count_nonzero(valid_mask)
    broadcast_valid_mask = numpy.broadcast_to(valid_mask, stack.shape)
    valid_stack = stack[broadcast_valid_mask].reshape(
        len(array_list), n_valid)
    result = numpy.empty(array_list[0].shape, dtype=numpy.float32)
    result[:] = NODATA
    result[valid_mask] = numpy.sum(valid_stack, axis=0)
    return result


def step_kernel(n_pixels, kernel_filepath):
    """Create a box kernel of 1+2*n_pixels."""

    kernel_size = 1+2*n_pixels
    driver = gdal.GetDriverByName('GTiff')
    kernel_dataset = driver.Create(
        kernel_filepath.encode('utf-8'), kernel_size, kernel_size, 1,
        gdal.GDT_Float32, options=[
            'BIGTIFF=IF_SAFER', 'TILED=YES', 'BLOCKXSIZE=256',
            'BLOCKYSIZE=256'])

    # Make some kind of geotransform, it doesn't matter what but
    # will make GIS libraries behave better if it's all defined
    kernel_dataset.SetGeoTransform([444720, 30, 0, 3751320, 0, -30])
    srs = osr.SpatialReference()
    srs.SetUTM(11, 1)
    srs.SetWellKnownGeogCS('NAD27')
    kernel_dataset.SetProjection(srs.ExportToWkt())

    kernel_band = kernel_dataset.GetRasterBand(1)
    kernel_band.SetNoDataValue(-9999)
    kernel_band.Fill(1)


def main():
    """Entry point."""

    if not os.path.exists(WORKSPACE_DIR):
        os.makedirs(WORKSPACE_DIR)

    crop_path_list = [
        os.path.join(BASE_DATA_DIR, "LUH2_1KM", "c3ann.flt"),
        os.path.join(BASE_DATA_DIR, "LUH2_1KM", "c3nfx.flt"),
        os.path.join(BASE_DATA_DIR, "LUH2_1KM", "c3per.flt"),
        os.path.join(BASE_DATA_DIR, "LUH2_1KM", "c4ann.flt"),
        os.path.join(BASE_DATA_DIR, "LUH2_1KM", "c4per.flt"),
    ]

    habitat_path_list = [
        os.path.join(BASE_DATA_DIR, "LUH2_1KM", "primf.flt"),
        os.path.join(BASE_DATA_DIR, "LUH2_1KM", "primn.flt"),
        os.path.join(BASE_DATA_DIR, "LUH2_1KM", "secdf.flt"),
        os.path.join(BASE_DATA_DIR, "LUH2_1KM", "secdn.flt"),
    ]

    grass_path_list = [
        os.path.join(BASE_DATA_DIR, "LUH2_1KM", "pastr.flt"),
        os.path.join(BASE_DATA_DIR, "LUH2_1KM", "range.flt")
    ]

    ag_proportion_path = os.path.join(WORKSPACE_DIR, 'ag_proportion.tif')

    nathabpath_id_map = {
        'nathab': os.path.join(WORKSPACE_DIR, 'nathab_proportion.tif'),
        'nathabgrass': os.path.join(
            WORKSPACE_DIR, 'nathabgrass_proportion.tif'),
    }

    task_graph = taskgraph.TaskGraph(
        os.path.join(WORKSPACE_DIR, 'taskgraph_cache'), 4)

    sumtask_id_map = {}
    for raster_path_list, target_path, task_id in [
            (crop_path_list, ag_proportion_path, 'ag'),
            (habitat_path_list, nathabpath_id_map['nathab'], 'nathab'),
            (habitat_path_list+grass_path_list,
             nathabpath_id_map['nathabgrass'], 'nathabgrass')]:
        sumtask_id_map[task_id] = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(x, 1) for x in raster_path_list], add_arrays, target_path,
                gdal.GDT_Float32, -9999),
            kwargs={'gtiff_creation_options': GTIFF_CREATION_OPTIONS},
            target_path_list=[target_path]
            )

    kernel_path = os.path.join(WORKSPACE_DIR, 'box_kernel.tif')
    kernel_task = task_graph.add_task(
        func=step_kernel,
        args=(2, kernel_path),
        target_path_list=[kernel_path])

    hab_in_2km_path_id_map = {
        'nathab': os.path.join(
            WORKSPACE_DIR, 'nathab_proportion_within_2km.tif'),
        'nathabgrass': os.path.join(
            WORKSPACE_DIR, 'nathabgrass_proportion_within_2km.tif'),
        }

    # classify ag_proportion_path into a binary mask
    classified_ag_path = os.path.join(
        WORKSPACE_DIR, 'ag_classified.tif')
    # task_list[0] is the class for the ag_proportion_path
    classify_cropland_task = task_graph.add_task(
        func=pygeoprocessing.raster_calculator,
        args=(
            [(ag_proportion_path, 1)],
            threshold_crop_raster, classified_ag_path,
            gdal.GDT_Byte, 2),
        kwargs={'gtiff_creation_options': GTIFF_CREATION_OPTIONS},
        target_path_list=[classified_ag_path],
        dependent_task_list=[sumtask_id_map['ag']]
        )

    masked_hab_path_id_map = {
        'nathab': os.path.join(
            WORKSPACE_DIR, 'nathab_agmasked_proportion_within_2km.tif'),
        'nathabgrass': os.path.join(
            WORKSPACE_DIR, 'nathabgrass_agmasked_proportion_within_2km.tif'),
    }

    for hab_id, hab_path in nathabpath_id_map.iteritems():
        # task_list[1] is the habitat classer task
        habitat_proportion_task = task_graph.add_task(
            func=pygeoprocessing.convolve_2d,
            args=(
                (hab_path, 1), (kernel_path, 1),
                hab_in_2km_path_id_map[hab_id]),
            kwargs={
                'ignore_nodata': True,
                'mask_nodata': True,
                'normalize_kernel': True,
                'target_datatype': gdal.GDT_Float32,
                'target_nodata': NODATA,
                'gtiff_creation_options': GTIFF_CREATION_OPTIONS,
                },
            target_path_list=[hab_in_2km_path_id_map[hab_id]],
            dependent_task_list=[kernel_task, sumtask_id_map[hab_id]])


        # task_list[0] is the task for ag_proportion_path
        mask_habitat_cropland_task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(classified_ag_path, 1),
                 (hab_in_2km_path_id_map[hab_id], 1)],
                mult_arrays, masked_hab_path_id_map[hab_id],
                gdal.GDT_Float32, -9999),
            kwargs={'gtiff_creation_options': GTIFF_CREATION_OPTIONS},
            target_path_list=[masked_hab_path_id_map[hab_id]],
            dependent_task_list=[
                classify_cropland_task, habitat_proportion_task]
            )

    task_graph.join()


if __name__ == '__main__':
    main()
