"""This pollination analysis is guided by this doc https://docs.google.com/document/d/1k5yyhisemNrjG7ZSFGfD7GbIGY3rcNWmNulL2Hcqr9Q/edit"""
import os
import logging
import subprocess

from osgeo import osr
from osgeo import gdal
import taskgraph
import pygeoprocessing
import numpy
import pandas


logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LUH2_BASE_DATA_DIR = r"C:\Users\Rich\Dropbox\ipbes-pollination-analysis\LUH2_1KM"
WORKSPACE_DIR = 'pollination_workspace'
BASE_CROP_DATA_DIR = r"C:\Users\Rich\Dropbox\Monfreda maps"

CROP_FILE_DIR = os.path.join(WORKSPACE_DIR, 'crop_geotiffs')
MICRONUTRIENT_DIR = os.path.join(WORKSPACE_DIR, 'micronutrient_working_files')
CROP_NUTRIENT_TABLE_PATH = os.path.join(BASE_CROP_DATA_DIR, "crop_info.csv")

GTIFF_CREATION_OPTIONS = ('TILED=YES', 'BIGTIFF=IF_SAFER', 'COMPRESS=DEFLATE')

NODATA = -9999

# any proportion of ag above this is classified as
CROP_THRESHOLD_VALUE = 0.05


def asc_to_tiff(base_path, target_path):
    """Convert base .asc.zip to target .tif w/ a WGS84 coordinate system."""
    tmp_asc_path = os.path.splitext(base_path)[0]
    if not os.path.exists(os.path.dirname(target_path)):
        os.makedirs(os.path.dirname(target_path))
    cmd = 'unzip -o "%s" -d "%s"' % (
        base_path, os.path.dirname(base_path))
    print cmd
    subprocess.call(cmd)
    with open(tmp_asc_path, 'rb') as base_file:
        for _ in xrange(6):
            base_file.readline()
        array = numpy.empty((2160, 4320), dtype=numpy.float32)
        running_array = []
        for line in base_file:
            running_array += [float(x) for x in ' '.join(line.split()).split()]
        print len(running_array), 2160*4320
        array[:] = numpy.array(
            [float(x) for x in running_array]).reshape(2160, 4320)
    os.remove(tmp_asc_path)
    driver = gdal.GetDriverByName('GTiff')
    target_raster = driver.Create(
        target_path.encode('utf-8'), 4320, 2160, 1, gdal.GDT_Float32,
        options=GTIFF_CREATION_OPTIONS)
    target_projection = osr.SpatialReference()
    target_projection.ImportFromEPSG(4326)
    target_raster.SetProjection(target_projection.ExportToWkt())
    target_raster.SetGeoTransform([-180.0000, 8.3333001E-02, 0, 90.00000, 0, -8.3333001E-02])
    target_band = target_raster.GetRasterBand(1)
    target_band.SetNoDataValue(-9999.0)
    target_band.WriteArray(array)
    target_band.FlushCache()
    target_band = None
    target_raster = None


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

    task_graph = taskgraph.TaskGraph(
        os.path.join(WORKSPACE_DIR, 'taskgraph_cache'), 0)

    crop_table = pandas.read_csv(CROP_NUTRIENT_TABLE_PATH)

    crop_table = crop_table[
        pandas.notnull(crop_table['crop'])]

    dep_pol_id_map = dict(zip(
        crop_table['crop'],
        crop_table['pol_dep']))

    percent_refuse_crop_map = dict(zip(
        crop_table['crop'],
        crop_table['Percentrefuse']/100.0))

    crop_yield_path_id_map = dict(
        [(crop_id, os.path.join(BASE_CROP_DATA_DIR, "%s_yield.asc.zip" % crop_id))
         for crop_id in dep_pol_id_map])

    crop_area_path_id_map = dict(
        [(crop_id, os.path.join(BASE_CROP_DATA_DIR, "%s_harea.asc.zip" % crop_id))
         for crop_id in dep_pol_id_map])

    target_crop_total_yield_path_id_map = dict(
        [(crop_id, os.path.join(CROP_FILE_DIR, "%s_tot_yield.tif" % crop_id))
         for crop_id in dep_pol_id_map])

    crop_yield_task_id_map = {}
    crop_area_task_id_map = {}
    for crop_id in dep_pol_id_map:
        base_crop_yield_path = crop_yield_path_id_map[crop_id]
        target_crop_yield_path = os.path.join(
            CROP_FILE_DIR, os.path.basename(
                base_crop_yield_path).split('.')[0] + '.tif')
        crop_yield_path_id_map[crop_id] = target_crop_yield_path
        crop_yield_task_id_map[crop_id] = task_graph.add_task(
            func=asc_to_tiff,
            args=(base_crop_yield_path, target_crop_yield_path),
            target_path_list=[target_crop_yield_path])

        base_crop_area_path = crop_area_path_id_map[crop_id]
        target_crop_area_path = os.path.join(
            CROP_FILE_DIR, os.path.basename(
                base_crop_area_path).split('.')[0] + '.tif')
        crop_area_path_id_map[crop_id] = target_crop_area_path
        crop_area_task_id_map[crop_id] = task_graph.add_task(
            func=asc_to_tiff,
            args=(base_crop_area_path, target_crop_area_path),
            target_path_list=[target_crop_area_path])

        total_yield_task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(target_crop_area_path, 1),
                 (target_crop_yield_path, 1)],
                mult_arrays, target_crop_total_yield_path_id_map[crop_id],
                gdal.GDT_Float32, -9999),
            kwargs={'gtiff_creation_options': GTIFF_CREATION_OPTIONS},
            target_path_list=[target_crop_total_yield_path_id_map[crop_id]],
            dependent_task_list=[
                crop_area_task_id_map[crop_id],
                crop_yield_task_id_map[crop_id]]
            )

    """Proportion of micronutrient production dependent on pollination
        For each crop (at 10 km):
            For Calories, Vitamin A, Fe, and Folate
            total micronutrient production = yield (EarthStat) x (100-percent refuse)/100 x area (EarthStat, convert proportion of gridcell to hectares) x micronutrient content per ton of crop
            pollinator-dependent micronutrient production = total micronutrient production x pollinator dependency ratio (0-0.95)
    """

    crop_path_list = [
        os.path.join(LUH2_BASE_DATA_DIR, "c3ann.flt"),
        os.path.join(LUH2_BASE_DATA_DIR, "c3nfx.flt"),
        os.path.join(LUH2_BASE_DATA_DIR, "c3per.flt"),
        os.path.join(LUH2_BASE_DATA_DIR, "c4ann.flt"),
        os.path.join(LUH2_BASE_DATA_DIR, "c4per.flt"),
    ]

    habitat_path_list = [
        os.path.join(LUH2_BASE_DATA_DIR, "primf.flt"),
        os.path.join(LUH2_BASE_DATA_DIR, "primn.flt"),
        os.path.join(LUH2_BASE_DATA_DIR, "secdf.flt"),
        os.path.join(LUH2_BASE_DATA_DIR, "secdn.flt"),
    ]

    grass_path_list = [
        os.path.join(LUH2_BASE_DATA_DIR, "pastr.flt"),
        os.path.join(LUH2_BASE_DATA_DIR, "range.flt")
    ]

    ag_proportion_path = os.path.join(WORKSPACE_DIR, 'ag_proportion.tif')

    nathabpath_id_map = {
        'nathab': os.path.join(WORKSPACE_DIR, 'nathab_proportion.tif'),
        'nathabgrass': os.path.join(
            WORKSPACE_DIR, 'nathabgrass_proportion.tif'),
    }

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

    task_graph.close()
    task_graph.join()


if __name__ == '__main__':
    main()
