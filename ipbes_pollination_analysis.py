"""This pollination analysis is guided by this doc https://docs.google.com/document/d/1k5yyhisemNrjG7ZSFGfD7GbIGY3rcNWmNulL2Hcqr9Q/edit"""
import os
import errno
import logging
import hashlib
import inspect
import datetime
import distutils.dir_util
import traceback

from osgeo import osr
from osgeo import gdal
import taskgraph
import pygeoprocessing
import numpy
import scipy.ndimage.morphology
import pandas

N_WORKERS = 4
POL_DEP_THRESHOLD = 0.3
GTIFF_CREATION_OPTIONS = ('TILED=YES', 'BIGTIFF=IF_SAFER', 'COMPRESS=LZW')
NODATA = -9999
MASK_NODATA = 2

# any proportion of ag above this is classified as a crop
CROP_THRESHOLD_VALUE = 0.05

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('ipbes_pollination_analysis')

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

GLOBIO_BASE_DATA_DIR = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes-data',
    'GLOBIO4_landuse_10sec_tifs_20171207_Idiv')
WORKSPACE_DIR = '%s_%s' % (
    os.path.basename(GLOBIO_BASE_DATA_DIR), '_pollination_dep_workspace')
BASE_CROP_DATA_DIR = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes-data', 'Monfreda maps')
BASE_CROP_RASTER_DIR = os.path.join(
    BASE_CROP_DATA_DIR, 'crop_rasters_as_geotiff')
CROP_NUTRIENT_TABLE_PATH = os.path.join(BASE_CROP_DATA_DIR, 'crop_info.csv')

GLOBIO_LU_MAPS = {
    '2015': os.path.join(
        GLOBIO_BASE_DATA_DIR, 'Current2015',
        'Globio4_landuse_10sec_2015_cropint.tif'),
    '2050_ssp1_rcp26': os.path.join(
        GLOBIO_BASE_DATA_DIR, 'SSP1_RCP26',
        'Globio4_landuse_10sec_2050_cropint.tif'),
    '2050_ssp3_rcp70': os.path.join(
        GLOBIO_BASE_DATA_DIR, 'SSP3_RCP70',
        'Globio4_landuse_10sec_2050_cropint.tif'),
    '2050_ssp5_rcp85': os.path.join(
        GLOBIO_BASE_DATA_DIR, 'SSP5_RCP85',
        'Globio4_landuse_10sec_2050_cropint.tif'),
}

GLOBIO_NATHAB_CODES = [6] + range(50, 181)
GLOBIO_SEMINAT_CODES = GLOBIO_NATHAB_CODES + [3, 4, 5]
GLOBIO_AG_CODES = [2, 230, 231, 232]

MICRONUTRIENT_LIST = ['Energy', 'VitA', 'Fe', 'Folate']

TARGET_GLOBIO_WORKING_DIR = os.path.join(
    WORKSPACE_DIR, 'globio_masks')
TARGET_CROP_FILE_DIR = os.path.join(WORKSPACE_DIR, 'crop_rasters')
TARGET_MICRONUTRIENT_YIELD_DIR = os.path.join(
    WORKSPACE_DIR, 'micronutrient_yield')
FINAL_TARGET_DIR = os.path.join(
    BASE_DROPBOX_DIR, 'rps_bck_shared_stuff', 'ipbes stuff', WORKSPACE_DIR)


def threshold_crop_raster(array):
    """Mask if non-nodata > CROP_THRESHOLD_VALUE."""
    result = numpy.zeros(array.shape, dtype=numpy.int8)
    result[:] = 2
    valid_mask = array != NODATA
    result[valid_mask] = array[valid_mask] > CROP_THRESHOLD_VALUE
    return result


def subtract_when_less_than_threshold_mask_op(
        base_array, delta_array, threshold_array):
    """Subtract delta from base when threshold < POL_DEP_THRESHOLD."""
    result = numpy.copy(base_array)
    threshold_mask = (
        (base_array != NODATA) &
        (delta_array != NODATA) &
        (threshold_array != NODATA) &
        (threshold_array < POL_DEP_THRESHOLD))
    result[threshold_mask] -= delta_array[threshold_mask]
    return result


def mask_geq_threshold_mask_op(base_array, threshold_array):
    """Mask base_array to nodata when threshold >= POL_DEP_THRESHOLD."""
    result = numpy.copy(base_array)
    result[threshold_array >= POL_DEP_THRESHOLD] = NODATA
    return result


def mult_arrays(*array_list):
    """Multiply arrays in list that don't have nodata in pixel stack."""
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
    """Add array stacks that don't have nodata in pixel stack."""
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


def add_arrays_passthrough_nodata(*array_list):
    """Add arrays in stack that have at least 1 valid value, nodata=0."""
    stack = numpy.stack(array_list)
    valid_mask = (
        numpy.bitwise_or.reduce(stack != NODATA, axis=0))
    n_valid = numpy.count_nonzero(valid_mask)
    broadcast_valid_mask = numpy.broadcast_to(valid_mask, stack.shape)
    valid_stack = stack[broadcast_valid_mask].reshape(
        len(array_list), n_valid)
    valid_stack[valid_stack == NODATA] = 0.0
    result = numpy.empty(array_list[0].shape, dtype=numpy.float32)
    result[:] = NODATA
    result[valid_mask] = numpy.sum(valid_stack, axis=0)
    return result


def step_kernel(n_pixels, kernel_filepath):
    """Create a box circle kernel of 1+2*n_pixels."""
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
    kernel_band.SetNoDataValue(NODATA)
    mask_array = numpy.ones((kernel_size, kernel_size))
    mask_array[n_pixels, n_pixels] = 0
    dist_array = scipy.ndimage.morphology.distance_transform_edt(mask_array)
    kernel_band.WriteArray(dist_array < n_pixels)


class MaskAtThreshold(object):
    """Pass base pixels if mask is < threshold."""
    def __init__(
            self, base_threshold_raster_path, base_raster_path, threshold,
            target_raster_path):
        """Pass through base raster where threshold is < threshold.

        Parameters:
            base_threshold_raster_path (str): path to a floating point raster
                used as the threshold cutoff.
            base_raster_path (str): path to base raster to pass through or
                mask if threshold raster is < threshold
            threshold (float): if a pixel in threshold_raster is below this
                value, pollinator dependent yield is subtracted from total.
                Otherwise total is passed through.
            target_raster_path (string): path to output raster.

        Returns:
            None.
        """
        try:
            self.__name__ = hashlib.sha1(inspect.getsource(
                MaskAtThreshold.__call__)).hexdigest()
        except IOError:
            # default to the classname if it doesn't work
            self.__name__ = MaskAtThreshold.__name__
        self.__name__ += str([
            base_threshold_raster_path, base_raster_path, threshold,
            target_raster_path])
        self.base_threshold_raster_path = base_threshold_raster_path
        self.base_raster_path = base_raster_path
        self.threshold = threshold
        self.target_raster_path = target_raster_path

    def __call__(self):
        try:
            os.makedirs(os.path.dirname(self.target_raster_path))
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        base_threshold_nodata = pygeoprocessing.get_raster_info(
            self.base_threshold_raster_path)['nodata'][0]

        def threshold_mask(
                threshold_array, base_array):
            """Pass through base if threshold < self.threshold."""
            result = numpy.empty_like(base_array)
            result[:] = NODATA
            threshold_mask = (
                (threshold_array != base_threshold_nodata) &
                (base_array != NODATA) &
                (threshold_array < self.threshold))
            result[threshold_mask] -= base_array[threshold_mask]
            return result

        pygeoprocessing.raster_calculator(
            [(self.base_threshold_raster_path, 1),
             (self.base_raster_path, 1)],
            threshold_mask, self.target_raster_path,
            gdal.GDT_Float32, NODATA,
            gtiff_creation_options=GTIFF_CREATION_OPTIONS)


class MicroNutrientRiskThreshold(object):
    """Subtract poll dep yield from total when LUH2 < threshold."""
    def __init__(
            self, base_raster_path_list, threshold, target_pol_dep_risk_path):
        """Mask subtract pollinator dependent yields from total.

        This happens when the fractional area of the functional type is less
        than the threshold.

        Parameters:
            base_raster_path_list (list): three raster paths:
                * LUH2 functional type proportion
                * total yield
                * pollinator dep yield.
                These rasters are all spatially aligned.
            threshold (float): if a piel in LUH2 is below this value,
                pollinator dependent yield is subtracted from total. Otherwise
                total is passed through.
            target_pol_dep_risk_path (string): path to output raster.

        Returns:
            None.
        """
        try:
            self.__name__ = hashlib.sha1(inspect.getsource(
                MicroNutrientRiskThreshold.__call__)).hexdigest()
        except IOError:
            # default to the classname if it doesn't work
            self.__name__ = MicroNutrientRiskThreshold.__name__
        self.__name__ += str([
            base_raster_path_list, threshold, target_pol_dep_risk_path])
        self.base_raster_path_list = base_raster_path_list
        self.threshold = threshold
        self.target_pol_dep_risk_path = target_pol_dep_risk_path

    def __call__(self):
        try:
            os.makedirs(os.path.dirname(self.target_pol_dep_risk_path))
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        luh2_nodata = pygeoprocessing.get_raster_info(
            self.base_raster_path_list[0])['nodata'][0]

        def calc_threshold_yield(
                luh2_array, total_yield_array, pol_dep_yield_array):
            """Pass through total yield unless luh2 < 0.3 then sub pol_dep."""
            result = numpy.empty_like(total_yield_array)
            result[:] = total_yield_array
            pol_dep_mask = (
                (luh2_array != luh2_nodata) &
                (pol_dep_yield_array != NODATA) &
                (luh2_array < self.threshold))
            result[pol_dep_mask] -= pol_dep_yield_array[pol_dep_mask]
            return result

        pygeoprocessing.raster_calculator(
            [(x, 1) for x in self.base_raster_path_list],
            calc_threshold_yield, self.target_pol_dep_risk_path,
            gdal.GDT_Float32, NODATA,
            gtiff_creation_options=GTIFF_CREATION_OPTIONS)


class MultRastersAndScalar(object):
    """Multiply list of rasters by themselves and a scalar."""
    def __init__(
            self, base_raster_path_list, scalar, target_path):
        try:
            self.__name__ = hashlib.sha1(inspect.getsource(
                MultRastersAndScalar.__call__)).hexdigest()
        except IOError:
            # default to the classname if it doesn't work
            self.__name__ = MultRastersAndScalar.__name__
        self.__name__ += str([
            base_raster_path_list, scalar, target_path])
        self.base_raster_path_list = base_raster_path_list
        self.target_path = target_path
        self.scalar = scalar

    def __call__(self):
        try:
            os.makedirs(os.path.dirname(self.target_path))
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        def mult_arrays_scalar(*array_list):
            stack = numpy.stack(array_list)
            valid_mask = (
                numpy.bitwise_and.reduce(stack != NODATA, axis=0))
            n_valid = numpy.count_nonzero(valid_mask)
            broadcast_valid_mask = numpy.broadcast_to(valid_mask, stack.shape)
            valid_stack = stack[broadcast_valid_mask].reshape(
                len(array_list), n_valid)
            result = numpy.empty(array_list[0].shape, dtype=numpy.float32)
            result[:] = NODATA
            result[valid_mask] = self.scalar * numpy.prod(
                valid_stack, axis=0)
            return result

        print '**********', self.target_path
        pygeoprocessing.raster_calculator(
            [(x, 1) for x in self.base_raster_path_list],
            mult_arrays_scalar, self.target_path,
            gdal.GDT_Float32, NODATA,
            gtiff_creation_options=GTIFF_CREATION_OPTIONS)


class MaskByRasterValue(object):
    """Mask a given raster by a given list of pixel values."""

    def __init__(
            self, base_raster_path_band, mask_list, target_path):
        try:
            self.__name__ = hashlib.sha1(inspect.getsource(
                MaskByRasterValue.__call__)).hexdigest()
        except IOError:
            # default to the classname if it doesn't work
            self.__name__ = MaskByRasterValue.__name__
        self.__name__ += str([
            base_raster_path_band, mask_list, target_path])
        self.base_raster_path_band = base_raster_path_band
        self.target_path = target_path
        self.mask_list = numpy.array(mask_list)

    def __call__(self):
        try:
            os.makedirs(os.path.dirname(self.target_path))
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        LOGGER.debug(
            "MaskByRasterValue\n\t%s\n\t%s\n\t%s", self.base_raster_path_band,
            self.mask_list, self.target_path)
        base_nodata = pygeoprocessing.get_raster_info(
            self.base_raster_path_band[0])['nodata'][
                self.base_raster_path_band[1]-1]

        def mask_values_op(base_array):
            valid_mask = base_array != base_nodata
            result = numpy.empty(base_array.shape, dtype=numpy.int8)
            result[:] = MASK_NODATA
            result[valid_mask] = numpy.in1d(
                base_array[valid_mask], self.mask_list)
            return result

        pygeoprocessing.raster_calculator(
            [self.base_raster_path_band],
            mask_values_op, self.target_path,
            gdal.GDT_Byte, MASK_NODATA,
            gtiff_creation_options=GTIFF_CREATION_OPTIONS)


def main():
    """Entry point."""
    LOGGER.info("starting %s" % __name__)
    if not os.path.exists(WORKSPACE_DIR):
        os.makedirs(WORKSPACE_DIR)
    if not os.path.exists(FINAL_TARGET_DIR):
        os.makedirs(FINAL_TARGET_DIR)
    readme_path = os.path.join(FINAL_TARGET_DIR, 'README.txt')
    readme_restarting = False
    if os.path.exists(readme_path):
        readme_restarting = True
    with open(os.path.join(
            FINAL_TARGET_DIR, 'README.txt'), 'a') as readme_file:
        if not readme_restarting:
            readme_file.write(
                "started `ipbes_pollination_analysis.py` on %s\n" %
                datetime.datetime.now())
        else:
            readme_file.write(
                "restarting `ipbes_pollination_analysis.py` on %s\n" %
                datetime.datetime.now())

    task_graph = taskgraph.TaskGraph(
        os.path.join(WORKSPACE_DIR, 'taskgraph_cache'), N_WORKERS)

    # for globio / ESA rasters
    globio_kernel_path = os.path.join(WORKSPACE_DIR, 'globio_box_kernel.tif')
    # 7 pixels since 2000 / 300 = 6.66
    globio_kernel_task = task_graph.add_task(
        func=step_kernel,
        args=(7, globio_kernel_path),
        target_path_list=[globio_kernel_path],
        task_name='globio_step_kernel')

    # Subtask I: Habitat Area Tasks
    habitat_area_path_task_map = {}
    for mask_hab_id, mask_list in [
            ('nathab', GLOBIO_NATHAB_CODES),
            ('seminat', GLOBIO_SEMINAT_CODES)]:
        for globio_raster_key, globio_raster_path in GLOBIO_LU_MAPS.iteritems():
            #Mask nathab/seminat/ag
            task_name = 'globio_%s_%s' % (mask_hab_id, globio_raster_key)
            globio_habmask_path = os.path.join(
                TARGET_GLOBIO_WORKING_DIR, '%s.tif' % task_name)
            globio_mask_task = task_graph.add_task(
                func=MaskByRasterValue(
                    (globio_raster_path, 1), mask_list, globio_habmask_path),
                target_path_list=[globio_habmask_path],
                task_name=task_name)

            globio_convolve_task_name = 'globio_prop_hab_in_2km_%s_%s' % (
                mask_hab_id, globio_raster_key)
            globio_hab_prop_path = os.path.join(
                TARGET_GLOBIO_WORKING_DIR,
                '%s.tif' % globio_convolve_task_name)
            habitat_proportion_task = task_graph.add_task(
                func=pygeoprocessing.convolve_2d,
                args=(
                    (globio_habmask_path, 1), (globio_kernel_path, 1),
                    globio_hab_prop_path),
                kwargs={
                    'ignore_nodata': False,
                    'mask_nodata': True,
                    'normalize_kernel': True,
                    'target_datatype': gdal.GDT_Float32,
                    'target_nodata': NODATA,
                    'gtiff_creation_options': GTIFF_CREATION_OPTIONS,
                    },
                target_path_list=[globio_hab_prop_path],
                dependent_task_list=[globio_kernel_task, globio_mask_task],
                task_name=globio_convolve_task_name)
            habitat_area_path_task_map[(mask_hab_id, globio_raster_key)] = (
                globio_hab_prop_path, habitat_proportion_task)

    # Subtask II: Nutrient Yield Tasks
    crop_table = pandas.read_csv(CROP_NUTRIENT_TABLE_PATH)
    crop_table = crop_table[pandas.notnull(crop_table['crop'])]

    dep_pol_id_map = dict(zip(
        crop_table['crop'],
        crop_table['pol_dep']))

    proportion_refuse_crop_map = dict(zip(
        crop_table['crop'],
        crop_table['Percentrefuse']/100.0))

    crop_yield_path_id_map = dict(
        [(crop_id, os.path.join(BASE_CROP_RASTER_DIR, '%s_yield.tif' % crop_id))
         for crop_id in dep_pol_id_map])

    harvested_proportion_path_id_map = dict(
        [(crop_id, os.path.join(BASE_CROP_RASTER_DIR, '%s_harea.tif' % crop_id))
         for crop_id in dep_pol_id_map])

    target_crop_grid_yield_path_id_map = dict(
        [(crop_id, os.path.join(
            TARGET_CROP_FILE_DIR, '%s_grid_yield.tif' % crop_id))
         for crop_id in dep_pol_id_map])

    pollinator_dependent_micronutrient_yield_path_id_map = {}
    total_micronutrient_yield_path_id_map = {}
    for crop_id in dep_pol_id_map:
        total_yield_task = task_graph.add_task(
            func=MultRastersAndScalar(
                [harvested_proportion_path_id_map[crop_id],
                 crop_yield_path_id_map[crop_id]],
                1.0-proportion_refuse_crop_map[crop_id],
                target_crop_grid_yield_path_id_map[crop_id]),
            target_path_list=[target_crop_grid_yield_path_id_map[crop_id]],
            task_name='MultRastersAndScalar_%s' % (
                target_crop_grid_yield_path_id_map[crop_id])
            )

        for micronutrient_id in MICRONUTRIENT_LIST:
            micronutrient_yield_path = os.path.splitext(
                target_crop_grid_yield_path_id_map[crop_id])[0] + (
                    '_%s.tif' % (micronutrient_id))

            # the 1e4 converts the Mg to g for nutrient units
            micronutrient_task = task_graph.add_task(
                func=MultRastersAndScalar(
                    [target_crop_grid_yield_path_id_map[crop_id]],
                    1e4 * float(crop_table[crop_table['crop'] == crop_id][
                        micronutrient_id]),
                    micronutrient_yield_path),
                target_path_list=[micronutrient_yield_path],
                dependent_task_list=[total_yield_task],
                task_name='MultRastersAndScalar_%s' % micronutrient_yield_path
                )

            # record the path and the task for later
            total_micronutrient_yield_path_id_map[
                (micronutrient_id, crop_id)] = (
                    micronutrient_yield_path, micronutrient_task)

            pollinator_dependent_micronutrient_yield_path = os.path.splitext(
                micronutrient_yield_path)[0] + '_%s.tif' % 'pol_dep'

            pollinator_dependent_micronutrient_task = task_graph.add_task(
                func=MultRastersAndScalar(
                    [micronutrient_yield_path],
                    float(crop_table[
                        crop_table['crop'] == crop_id]['pol_dep']),
                    pollinator_dependent_micronutrient_yield_path),
                target_path_list=[
                    pollinator_dependent_micronutrient_yield_path],
                dependent_task_list=[micronutrient_task],
                task_name='MultRastersAndScalar'
                )
            pollinator_dependent_micronutrient_yield_path_id_map[
                (micronutrient_id, crop_id)] = (
                    pollinator_dependent_micronutrient_yield_path,
                    pollinator_dependent_micronutrient_task)

    micro_yield_path_task_map = {}
    for micronutrient_id in MICRONUTRIENT_LIST:
        for pol_dep_id, micronutrient_yield_task_map in [(
                '', total_micronutrient_yield_path_id_map),
                ('_pol_dep', pollinator_dependent_micronutrient_yield_path_id_map)]:

            micronutrient_crop_path_list = [
                (micronutrient_yield_task_map[(
                    micronutrient_id, crop_id)][0], 1)
                for crop_id in dep_pol_id_map]
            micronutrient_crop_task_list = [
                micronutrient_yield_task_map[(micronutrient_id, crop_id)][1]
                for crop_id in dep_pol_id_map]

            total_micro_yield_path = os.path.join(
                TARGET_CROP_FILE_DIR, 'total%s_yield_%s.tif' % (
                    pol_dep_id, micronutrient_id))

            micronutrient_crop_sum_task = task_graph.add_task(
                func=add_arrays_passthrough_nodata,
                args=(
                    micronutrient_crop_path_list,
                    add_arrays_passthrough_nodata,
                    total_micro_yield_path, gdal.GDT_Float32, NODATA),
                target_path_list=[total_micro_yield_path],
                dependent_task_list=micronutrient_crop_task_list,
                task_name='add_total%s_crops_%s' % (
                    pol_dep_id, micronutrient_id))
            micro_yield_path_task_map[(micronutrient_id, pol_dep_id)] = (
                total_micro_yield_path, micronutrient_crop_sum_task)

    # Subtask III: micronutrient yield to GLOBIO scenarios
    # mask ag from all 4 globio lulc maps
    globio_ag_mask_path_task_map = {}
    for globio_raster_key, globio_raster_path in GLOBIO_LU_MAPS.iteritems():
        task_name = 'globio_ag_mask_%s' % (globio_raster_key)
        globio_agmask_path = os.path.join(
            TARGET_GLOBIO_WORKING_DIR, '%s.tif' % task_name)
        globio_ag_mask_task = task_graph.add_task(
            func=MaskByRasterValue(
                (globio_raster_path, 1), GLOBIO_AG_CODES,
                globio_agmask_path),
            target_path_list=[globio_agmask_path],
            task_name=task_name)
        globio_ag_mask_path_task_map[globio_raster_key] = (
            globio_agmask_path, globio_ag_mask_task)

    # project micronutrient yields onto ag maps
    base_globio_info = pygeoprocessing.get_raster_info(
        GLOBIO_LU_MAPS.itervalues().next())
    warp_micronutrient_path_task_map = {}
    for micronutrient_id in MICRONUTRIENT_LIST:
        for pol_dep_id in ['', '_pol_dep']:
            base_micronutrient_yield_path = micro_yield_path_task_map[
                (micronutrient_id, pol_dep_id)][0]
            target_aligned_micronutrient_yield_path = os.path.join(
                TARGET_MICRONUTRIENT_YIELD_DIR, 'globioaligned_%s.tif' % (
                    os.path.splitext(
                        os.path.basename(base_micronutrient_yield_path))[0]))

            warp_micronutrient_task = task_graph.add_task(
                func=pygeoprocessing.warp_raster,
                args=(
                    base_micronutrient_yield_path,
                    base_globio_info['pixel_size'],
                    target_aligned_micronutrient_yield_path,
                    'nearest'),
                kwargs={
                    'target_bb': base_globio_info['bounding_box'],
                    'target_sr_wkt': base_globio_info['projection']
                    },
                target_path_list=[target_aligned_micronutrient_yield_path],
                dependent_task_list=[
                    micro_yield_path_task_map[
                        (micronutrient_id, pol_dep_id)][1]],
                task_name='warp_micronutrient_%s%s' % (
                    micronutrient_id, pol_dep_id))
            warp_micronutrient_path_task_map[
                (micronutrient_id, pol_dep_id)] = (
                    target_aligned_micronutrient_yield_path,
                    warp_micronutrient_task)

    # multiply warp_micronutrient_path with globio_ag_mask_path
    total_micronutrient_path_task_map = {}
    for globio_raster_key in GLOBIO_LU_MAPS:
        globio_ag_mask_path, globio_ag_mask_task = (
            globio_ag_mask_path_task_map[globio_raster_key])
        for micronutrient_id in MICRONUTRIENT_LIST:
            for pol_dep_id in ['', '_pol_dep']:
                warp_micronutrient_path, warp_micronutrient_task = (
                    warp_micronutrient_path_task_map[
                        (micronutrient_id, pol_dep_id)])

                target_total_micronutrient_yield_path = os.path.join(
                    TARGET_MICRONUTRIENT_YIELD_DIR,
                    'total_%s%s_yield_%s.tif' % (
                        globio_raster_key, pol_dep_id, micronutrient_id))

                total_micronutrient_yield_task = task_graph.add_task(
                    func=pygeoprocessing.raster_calculator,
                    args=(
                        [(globio_ag_mask_path, 1),
                         (warp_micronutrient_path, 1)], mult_arrays,
                        target_total_micronutrient_yield_path, NODATA,),
                    target_path_list=[target_total_micronutrient_yield_path],
                    dependent_task_list=[
                        globio_ag_mask_task, warp_micronutrient_task],
                    task_name='total_%s%s_yield_%s' % (
                        globio_raster_key, pol_dep_id, micronutrient_id))
                total_micronutrient_path_task_map[
                    (globio_raster_key, micronutrient_id, pol_dep_id)] = (
                        target_total_micronutrient_yield_path,
                        total_micronutrient_yield_task)

    # Primary output of analysis 'total_%s_%s_stable_yield_pointthree_%s
    for mask_hab_id in ['nathab', 'seminat']:
        for globio_raster_key in GLOBIO_LU_MAPS:
            globio_ag_mask_path, globio_ag_mask_task = (
                globio_ag_mask_path_task_map[globio_raster_key])
            hab_area_path, hab_area_task = habitat_area_path_task_map[
                (mask_hab_id, globio_raster_key)]
            for micronutrient_id in MICRONUTRIENT_LIST:
                total_yield_path, total_yield_task = (
                    total_micronutrient_path_task_map[(
                        globio_raster_key, micronutrient_id, '')])
                pol_dep_total_yield_path, pol_dep_total_yield_task = (
                    total_micronutrient_path_task_map[(
                        globio_raster_key, micronutrient_id, '_pol_dep')])

                target_total_stable_yield_path = os.path.join(
                    WORKSPACE_DIR,
                    'total_%s_%s_stable_yield_pointthree_%s.tif' % (
                        globio_raster_key, mask_hab_id, micronutrient_id))

                task_graph.add_task(
                    func=pygeoprocessing.raster_calculator,
                    args=(
                        [(total_yield_path, 1), (pol_dep_total_yield_path, 1),
                         (globio_ag_mask_path, 1)],
                        subtract_when_less_than_threshold_mask_op,
                        target_total_stable_yield_path, gdal.GDT_Float32,
                        NODATA),
                    target_path_list=[target_total_stable_yield_path],
                    dependent_task_list=[
                        globio_ag_mask_task, pol_dep_total_yield_task,
                        total_yield_task],
                    task_name='total_%s_%s_stable_yield_pointthree_%s' % (
                        globio_raster_key, mask_hab_id, micronutrient_id))

                target_total_pol_dep_yield_loss_path = os.path.join(
                    WORKSPACE_DIR,
                    'total_%s_%s_pol_dep_yield_loss_pointthree_%s.tif' % (
                        globio_raster_key, mask_hab_id, micronutrient_id))

                task_graph.add_task(
                    func=pygeoprocessing.raster_calculator,
                    args=(
                        [(pol_dep_total_yield_path, 1),
                         (globio_ag_mask_path, 1)],
                        mask_geq_threshold_mask_op,
                        target_total_pol_dep_yield_loss_path,
                        gdal.GDT_Float32, NODATA),
                    target_path_list=[target_total_pol_dep_yield_loss_path],
                    dependent_task_list=[
                        globio_ag_mask_task, pol_dep_total_yield_task],
                    task_name=(
                        'total_%s_%s_pol_dep_yield_loss_pointthree_%s' % (
                            globio_raster_key, mask_hab_id,
                            micronutrient_id)))

    LOGGER.info("closing and joining taskgraph")
    task_graph.close()

    with open(os.path.join(
            FINAL_TARGET_DIR, 'README.txt'), 'a') as readme_file:
        readme_file.write(
            "taskgraph scheduled and closed on %s now we wait\n" %
            datetime.datetime.now())

    task_graph.join()

    with open(os.path.join(
            FINAL_TARGET_DIR, 'README.txt'), 'a') as readme_file:
        readme_file.write(
            "copying complete files to this directory on %s\n" %
            datetime.datetime.now())

    distutils.dir_util.copy_tree(WORKSPACE_DIR, FINAL_TARGET_DIR)

    with open(os.path.join(
            FINAL_TARGET_DIR, 'README.txt'), 'a') as readme_file:
        readme_file.write(
            "everything done on %s\n" % datetime.datetime.now())

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        with open(
                os.path.join(FINAL_TARGET_DIR, 'README.txt'),
                'a') as readme_file:
            readme_file.write(
                "something crashed! on %s\n" % datetime.datetime.now())
            readme_file.write(
                "here's the error: %s\n" % traceback.format_exec())
            readme_file.write(
                "program is ending now you won't see anything else until "
                "it restarts\n")
