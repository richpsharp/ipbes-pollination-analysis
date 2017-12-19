"""This pollination analysis is guided by this doc https://docs.google.com/document/d/1k5yyhisemNrjG7ZSFGfD7GbIGY3rcNWmNulL2Hcqr9Q/edit"""
import os
import errno
import logging
import hashlib
import inspect
import datetime

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
for path in POSSIBLE_DROPBOX_LOCATIONS:
    if os.path.exists(path):
        BASE_DROPBOX_DIR = path
        break
LUH2_BASE_DATA_DIR = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes-data', 'LUH2_1KM_ag_and_cover_as_geotiff')
GLOBIO_BASE_DATA_DIR = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes-data',
    'GLOBIO4_landuse_10sec_tifs_20171207_Idiv')
WORKSPACE_DIR = os.path.join(BASE_DROPBOX_DIR, 'ipbes_pollination_analysis')
BASE_CROP_DATA_DIR = os.path.join(
    BASE_DROPBOX_DIR, 'ipbes-data', 'Monfreda maps')
BASE_CROP_RASTER_DIR = os.path.join(
    BASE_CROP_DATA_DIR, 'crop_rasters_as_geotiff')
CROP_NUTRIENT_TABLE_PATH = os.path.join(BASE_CROP_DATA_DIR, 'crop_info.csv')
CROP_CATEGORIES_TABLE_PATH = os.path.join(
    BASE_CROP_DATA_DIR, 'earthstat_to_luh_categories.csv')

GLOBIO_LU_MAPS = {
    '2015': os.path.join(
        GLOBIO_BASE_DATA_DIR, 'Current2015',
        'Globio4_landuse_10sec_2015_cropint.tif'),
    '2050_ssp1_rcp26': os.path.join(
        GLOBIO_BASE_DATA_DIR, 'SSP1_RCP26',
        'Globio4_landuse_10sec_2050_cropint.tif')
    '2050_ssp3_rcp70': : os.path.join(
        GLOBIO_BASE_DATA_DIR, 'SSP3_RCP70',
        'Globio4_landuse_10sec_2050_cropint.tif')
    '2050_ssp5_rcp85': os.path.join(
        GLOBIO_BASE_DATA_DIR, 'SSP5_RCP85',
        'Globio4_landuse_10sec_2050_cropint.tif')
}

GLOBIO_NATHAB_CODES = [6] + range(50, 181)
GLOBIO_SEMINAT_CODES = GLOBIO_NATHAB_CODES + [3, 4, 5]
GLOBIO_AG_CODES = [2, 230, 231, 232]

MICRONUTRIENT_LIST = ['Energy', 'VitA', 'Fe', 'Folate']

RASTER_FUNCTIONAL_TYPE_MAP = {
    ("C3", "annual"): os.path.join(LUH2_BASE_DATA_DIR, 'c3ann.tif'),
    ("C3", "perennial"): os.path.join(LUH2_BASE_DATA_DIR, 'c3per.tif'),
    ("C4", "annual"): os.path.join(LUH2_BASE_DATA_DIR, 'c4ann.tif'),
    ("C4", "perennial"): os.path.join(LUH2_BASE_DATA_DIR, 'c4per.tif'),
    ("N-fixer", None): os.path.join(LUH2_BASE_DATA_DIR, 'c3nfx.tif')
}

TARGET_GLOBIO_WORKING_DIR = os.path.join(
    WORKSPACE_DIR, 'globio_masks')
TARGET_CROP_FILE_DIR = os.path.join(
    WORKSPACE_DIR, 'crop_geotiffs')
TARGET_MICRONUTRIENT_DIR = os.path.join(
    WORKSPACE_DIR, 'micronutrient_working_files')


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


def add_arrays_passthrough_nodata(*array_list):
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
    def __init__(
            self, base_raster_path_list, threshold, target_pol_dep_risk_path):
        """Mask subtract pollinator dependent yields from total.

        This happens when the fractional area of the functional type is less
        than the threshold.

        Parameters:
            base_raster_path_list (list): list of base paths in order of LUH2
                functional type proportion, total yield, pollinator dep yield.
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

    with open(os.path.join(WORKSPACE_DIR, 'README.txt'), 'w') as readme_file:
        readme_file.write(
            "output of `ipbes_pollination_analysis.py` on %s" %
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

            globio_convolve_task_name = 'globio_prop_hab_%s_%s' % (
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
    #Convolve nathab/seminat
    #Mask convolution of nathab/seminat w/ ag

    crop_table = pandas.read_csv(CROP_NUTRIENT_TABLE_PATH)

    # First we do this:
    """Proportion of micronutrient production dependent on pollination
        For each crop (at 10 km):
            For Calories, Vitamin A, Fe, and Folate
            total micronutrient production = yield (EarthStat) x (100-percent refuse)/100 x area (EarthStat, convert proportion of gridcell to hectares) x micronutrient content per ton of crop
            pollinator-dependent micronutrient production = total micronutrient production x pollinator dependency ratio (0-0.95)
    """

    # this filters out the first partially blank row that's used for unit
    # notation
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

    crop_area_path_id_map = dict(
        [(crop_id, os.path.join(BASE_CROP_RASTER_DIR, '%s_harea.tif' % crop_id))
         for crop_id in dep_pol_id_map])

    target_crop_total_yield_path_id_map = dict(
        [(crop_id, os.path.join(
            WORKSPACE_DIR, 'crop_yields', '%s_tot_yield.tif' % crop_id))
         for crop_id in dep_pol_id_map])

    pollinator_dependent_micronutrient_yield_path_id_map = {}
    total_micronutrient_yield_path_id_map = {}
    for crop_id in dep_pol_id_map:
        total_yield_task = task_graph.add_task(
            func=MultRastersAndScalar(
                [crop_area_path_id_map[crop_id],
                 crop_yield_path_id_map[crop_id]],
                1.0-proportion_refuse_crop_map[crop_id],
                target_crop_total_yield_path_id_map[crop_id]),
            target_path_list=[target_crop_total_yield_path_id_map[crop_id]],
            task_name='MultRastersAndScalar_%s' % (
                target_crop_total_yield_path_id_map[crop_id])
            )

        for micronutrient_id in MICRONUTRIENT_LIST:
            micronutrient_yield_path = os.path.splitext(
                target_crop_total_yield_path_id_map[crop_id])[0] + (
                    '_%s.tif' % (micronutrient_id))

            # the 1e4 converts the Mg to g for nutrient units
            micronutrient_task = task_graph.add_task(
                func=MultRastersAndScalar(
                    [target_crop_total_yield_path_id_map[crop_id]],
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

    # Now we do this:
    """Sum up all the crops in each functional group: c3ann, c3per, c4ann, c4per, c3nfx (Becky to classify in table) = c3ann vitamin A total and pollinator dependent production - for current, at 10 km resolution"""
    crop_categories_table = pandas.read_csv(CROP_CATEGORIES_TABLE_PATH)

    crop_id_functional_type_map = {}

    for c_type, period in RASTER_FUNCTIONAL_TYPE_MAP:
        crop_filter_series = crop_categories_table['c_type'] == c_type
        if period is not None:
            # only the N-Fix is missing a period
            crop_filter_series &= (crop_categories_table['period'] == period)
        crop_id_functional_type_map[(c_type, period)] = list(
            crop_categories_table[crop_filter_series][
                'earthstat_filename_prefix'])

    micronutrient_functional_yield_map = {}
    for (c_type, period), crop_id_list in (
            crop_id_functional_type_map.iteritems()):
        for micronutrient_id in MICRONUTRIENT_LIST:
            # calculate both the total and pollinator dependent micronutrient
            pol_dep_micronutrient_functional_yield_path = os.path.join(
                WORKSPACE_DIR, 'functional_group_yields',
                '%s_%s_%s_yield_pol_dep.tif' % (
                    micronutrient_id, c_type, period))

            # create a directory for the ouput file if it doesn't exist
            try:
                os.makedirs(os.path.dirname(
                    pol_dep_micronutrient_functional_yield_path))
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

            pollinator_functional_crop_path_list = [
                pollinator_dependent_micronutrient_yield_path_id_map[
                    (micronutrient_id, crop_id)][0]
                for crop_id in crop_id_list if crop_id in dep_pol_id_map]
            pollinator_functional_crop_task_list = [
                pollinator_dependent_micronutrient_yield_path_id_map[
                    (micronutrient_id, crop_id)][1]
                for crop_id in crop_id_list if crop_id in dep_pol_id_map]

            pol_dep_micro_task = task_graph.add_task(
                func=pygeoprocessing.raster_calculator,
                args=(
                    [(x, 1) for x in pollinator_functional_crop_path_list],
                    add_arrays_passthrough_nodata,
                    pol_dep_micronutrient_functional_yield_path,
                    gdal.GDT_Float32, NODATA),
                kwargs={'gtiff_creation_options': GTIFF_CREATION_OPTIONS},
                target_path_list=[pol_dep_micronutrient_functional_yield_path],
                dependent_task_list=pollinator_functional_crop_task_list,
                task_name='raster_calculator_sum_pol_dep_micronutrient'
            )

            total_micronutrient_functional_yield_path = os.path.join(
                WORKSPACE_DIR, 'functional_group_yields',
                '%s_%s_%s_yield_total.tif' % (
                    micronutrient_id, c_type, period))

            total_functional_crop_path_list = [
                total_micronutrient_yield_path_id_map[
                    (micronutrient_id, crop_id)][0]
                for crop_id in crop_id_list if crop_id in dep_pol_id_map]

            total_functional_crop_task_list = [
                total_micronutrient_yield_path_id_map[
                    (micronutrient_id, crop_id)][1]
                for crop_id in crop_id_list if crop_id in dep_pol_id_map]

            total_micro_task = task_graph.add_task(
                func=pygeoprocessing.raster_calculator,
                args=(
                    [(x, 1) for x in total_functional_crop_path_list],
                    add_arrays_passthrough_nodata,
                    total_micronutrient_functional_yield_path,
                    gdal.GDT_Float32, NODATA),
                kwargs={'gtiff_creation_options': GTIFF_CREATION_OPTIONS},
                target_path_list=[total_micronutrient_functional_yield_path],
                dependent_task_list=total_functional_crop_task_list,
                task_name='raster_calculator_sum_total_micronutrient'
            )

            # we'll use this later when combining with the proportion of
            # natural habitat and grassland habitat
            micronutrient_functional_yield_map[
                (micronutrient_id, c_type, period)] = (
                    pol_dep_micronutrient_functional_yield_path,
                    pol_dep_micro_task,
                    total_micronutrient_functional_yield_path,
                    total_micro_task)

    crop_path_list = [
        os.path.join(LUH2_BASE_DATA_DIR, "c3ann.tif"),
        os.path.join(LUH2_BASE_DATA_DIR, "c3nfx.tif"),
        os.path.join(LUH2_BASE_DATA_DIR, "c3per.tif"),
        os.path.join(LUH2_BASE_DATA_DIR, "c4ann.tif"),
        os.path.join(LUH2_BASE_DATA_DIR, "c4per.tif"),
    ]

    habitat_path_list = [
        os.path.join(LUH2_BASE_DATA_DIR, "primf.tif"),
        os.path.join(LUH2_BASE_DATA_DIR, "primn.tif"),
        os.path.join(LUH2_BASE_DATA_DIR, "secdf.tif"),
        os.path.join(LUH2_BASE_DATA_DIR, "secdn.tif"),
    ]

    grass_path_list = [
        os.path.join(LUH2_BASE_DATA_DIR, "pastr.tif"),
        os.path.join(LUH2_BASE_DATA_DIR, "range.tif")
    ]

    ag_proportion_path = os.path.join(WORKSPACE_DIR, 'ag_proportion.tif')

    nathabpath_id_map = {
        'nathab': os.path.join(WORKSPACE_DIR, 'nathab_proportion.tif'),
        'seminat': os.path.join(
            WORKSPACE_DIR, 'seminat_proportion.tif'),
    }

    sumtask_id_map = {}
    for raster_path_list, target_path, task_id in [
            (crop_path_list, ag_proportion_path, 'ag'),
            (habitat_path_list, nathabpath_id_map['nathab'], 'nathab'),
            (habitat_path_list+grass_path_list,
             nathabpath_id_map['seminat'], 'seminat')]:
        sumtask_id_map[task_id] = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(x, 1) for x in raster_path_list], add_arrays, target_path,
                gdal.GDT_Float32, NODATA),
            kwargs={'gtiff_creation_options': GTIFF_CREATION_OPTIONS},
            target_path_list=[target_path],
            task_name='raster_calculator_add_arrays'
            )

    kernel_path = os.path.join(WORKSPACE_DIR, 'box_kernel.tif')
    kernel_task = task_graph.add_task(
        func=step_kernel,
        args=(2, kernel_path),
        target_path_list=[kernel_path],
        task_name='step_kernel')

    hab_in_2km_path_id_map = {
        'nathab': os.path.join(
            WORKSPACE_DIR, 'nathab_proportion_within_2km.tif'),
        'seminat': os.path.join(
            WORKSPACE_DIR, 'seminat_proportion_within_2km.tif'),
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
        dependent_task_list=[sumtask_id_map['ag']],
        task_name='raster_calculator'
        )

    masked_hab_path_id_map = {
        'nathab': os.path.join(
            WORKSPACE_DIR, 'nathab_agmasked_proportion_within_2km.tif'),
        'seminat': os.path.join(
            WORKSPACE_DIR, 'seminat_agmasked_proportion_within_2km.tif'),
    }
    for hab_id, hab_path in nathabpath_id_map.iteritems():
        # task_list[1] is the habitat classer task
        habitat_proportion_task = task_graph.add_task(
            func=pygeoprocessing.convolve_2d,
            args=(
                (hab_path, 1), (kernel_path, 1),
                hab_in_2km_path_id_map[hab_id]),
            kwargs={
                'ignore_nodata': False,
                'mask_nodata': True,
                'normalize_kernel': True,
                'target_datatype': gdal.GDT_Float32,
                'target_nodata': NODATA,
                'gtiff_creation_options': GTIFF_CREATION_OPTIONS,
                },
            target_path_list=[hab_in_2km_path_id_map[hab_id]],
            dependent_task_list=[kernel_task, sumtask_id_map[hab_id]],
            task_name='convolve_2d')

        # task_list[0] is the task for ag_proportion_path
        mask_habitat_cropland_task = task_graph.add_task(
            func=pygeoprocessing.raster_calculator,
            args=(
                [(classified_ag_path, 1),
                 (hab_in_2km_path_id_map[hab_id], 1)],
                mult_arrays, masked_hab_path_id_map[hab_id],
                gdal.GDT_Float32, NODATA),
            kwargs={'gtiff_creation_options': GTIFF_CREATION_OPTIONS},
            target_path_list=[masked_hab_path_id_map[hab_id]],
            dependent_task_list=[
                classify_cropland_task, habitat_proportion_task],
            task_name='raster_calculator'
            )

    # At 1 km LUH2 data, subtract pollinator dependent from total micronutrient production for any grid cells < 0.30; keep it at total for micronutrient production for any grid cells >0.30
    for micronutrient_key in micronutrient_functional_yield_map:
        (pol_dep_micronutrient_functional_yield_path,
         pol_dep_micro_task,
         total_micronutrient_functional_yield_path,
         total_micro_task) = micronutrient_functional_yield_map[
            micronutrient_key]
        """micronutrient_functional_yield_map[
            (micronutrient_id, c_type, period)] = (
                pol_dep_micronutrient_functional_yield_path,
                pol_dep_micro_task,
                total_micronutrient_functional_yield_path,
                total_micro_task)"""

        # put aligned rasters in a subdirectory
        base_luh_raster_info = pygeoprocessing.get_raster_info(
            crop_path_list[0])
        target_aligned_rasters = [
            os.path.join(
                WORKSPACE_DIR, 'aligned_rasters',
                'aligned_%s_rasters' % str(micronutrient_key),
                os.path.basename(x)) for x in [
                    nathabpath_id_map['nathab'],
                    nathabpath_id_map['seminat'],
                    total_micronutrient_functional_yield_path,
                    pol_dep_micronutrient_functional_yield_path]]

        # clip and align the yield and LUH2 rasters to be the intersection
        # of bounding boxes and the pixel size of LUH2
        align_raster_task = task_graph.add_task(
            func=pygeoprocessing.align_and_resize_raster_stack,
            args=(
                [nathabpath_id_map['nathab'],
                 nathabpath_id_map['seminat'],
                 total_micronutrient_functional_yield_path,
                 pol_dep_micronutrient_functional_yield_path],
                target_aligned_rasters,
                ['nearest']*4, base_luh_raster_info['pixel_size'],
                'intersection'),
            kwargs={
                'raster_align_index': 0,
                'gtiff_creation_options': GTIFF_CREATION_OPTIONS},
            target_path_list=target_aligned_rasters,
            dependent_task_list=sumtask_id_map.values()+[
                pol_dep_micro_task, total_micro_task],
            task_name='align_resize_raster_%s' % str(
                micronutrient_key))

        for hab_id, raster_index in [('nathab', 0), ('seminat', 1)]:
            hab_risk_path = os.path.join(
                WORKSPACE_DIR, '%s_risk_%s_yield.tif' % (
                    hab_id, micronutrient_key))

            # the 0.3 comes from Becky saying that's where we should threshold
            nathab_risk_task = task_graph.add_task(
                func=MicroNutrientRiskThreshold(
                    [target_aligned_rasters[raster_index]] +
                    target_aligned_rasters[2:], POL_DEP_THRESHOLD,
                    hab_risk_path),
                target_path_list=[hab_risk_path],
                dependent_task_list=[align_raster_task],
                task_name='MicroNutrientRiskThreshold_%s' % str(
                    (hab_id, micronutrient_key)))

            # calculate the pure loss
            pol_dep_loss_path = os.path.join(
                WORKSPACE_DIR, '%s_pol_dep_loss_%s_yield.tif' % (
                    hab_id, micronutrient_key))
            pol_dep_mask = task_graph.add_task(
                func=MaskAtThreshold(
                    target_aligned_rasters[raster_index],
                    target_aligned_rasters[3],
                    POL_DEP_THRESHOLD, pol_dep_loss_path),
                target_path_list=[pol_dep_loss_path],
                dependent_task_list=[align_raster_task],
                task_name='MaskAtThreshold_%s' % str(
                    (hab_id, micronutrient_key)))

    LOGGER.info("closing and joining taskgraph")
    task_graph.close()
    task_graph.join()


if __name__ == '__main__':
    main()
