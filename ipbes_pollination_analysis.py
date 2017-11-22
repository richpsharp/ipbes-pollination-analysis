"""This pollination analysis is guided by this doc https://docs.google.com/document/d/1k5yyhisemNrjG7ZSFGfD7GbIGY3rcNWmNulL2Hcqr9Q/edit"""
import os
import logging
import subprocess
import hashlib
import inspect

from osgeo import osr
from osgeo import gdal
import taskgraph
import pygeoprocessing
import numpy
import pandas

N_WORKERS = 4

logging.basicConfig(
    format='%(asctime)s %(name)-10s %(levelname)-8s %(message)s',
    level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('ipbes_pollination_analysis')

LUH2_BASE_DATA_DIR = r"C:\Users\Rich\Dropbox\ipbes-pollination-analysis\LUH2_1KM"
WORKSPACE_DIR = 'pollination_workspace'
BASE_CROP_DATA_DIR = r"C:\Users\Rich\Dropbox\Monfreda maps"
CROP_NUTRIENT_TABLE_PATH = os.path.join(BASE_CROP_DATA_DIR, "crop_info.csv")
CROP_CATEGORIES_TABLE_PATH = os.path.join(
    BASE_CROP_DATA_DIR, "earthstat_to_luh_categories.csv")

MICRONUTRIENT_LIST = ['Energy', 'VitA', 'Fe', 'Folate']

RASTER_FUNCTIONAL_TYPE_MAP = {
    ("C3", "annual"): os.path.join(LUH2_BASE_DATA_DIR, 'c3ann.flt'),
    ("C3", "perennial"): os.path.join(LUH2_BASE_DATA_DIR, 'c3per.flt'),
    ("C4", "annual"): os.path.join(LUH2_BASE_DATA_DIR, 'c4ann.flt'),
    ("C4", "perennial"): os.path.join(LUH2_BASE_DATA_DIR, 'c4per.flt'),
    ("N-fixer", None): os.path.join(LUH2_BASE_DATA_DIR, 'c3nfx.flt')
}

TARGET_CROP_FILE_DIR = os.path.join(WORKSPACE_DIR, 'crop_geotiffs')
TARGET_MICRONUTRIENT_DIR = os.path.join(WORKSPACE_DIR, 'micronutrient_working_files')

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
        # this is hard-coded from the Montfrida maps so putting assert here
        # in case I blindly copy it somewhere else someday
        assert(len(running_array) == 2160*4320)
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
    target_raster.SetGeoTransform([
        -180.0000, 8.3333001E-02, 0, 90.00000, 0, -8.3333001E-02])
    target_band = target_raster.GetRasterBand(1)
    target_band.SetNoDataValue(NODATA)
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
    kernel_band.SetNoDataValue(NODATA)
    kernel_band.Fill(1)


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

        luh2_nodata = pygeoprocessing.get_raster_info(
            self.base_raster_path_list[0])['nodata'][0]

        def calc_threshold_yield(
                luh2_array, total_yield_array, pol_dep_yield_array):
            """Pass through total yield unless luh2 < 0.3 then sub pol_dep."""
            result = numpy.empty_like(total_yield_array)
            pol_dep_mask = (
                (luh2_array != luh2_nodata) &
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

def main():
    """Entry point."""
    if not os.path.exists(WORKSPACE_DIR):
        os.makedirs(WORKSPACE_DIR)

    task_graph = taskgraph.TaskGraph(
        os.path.join(WORKSPACE_DIR, 'taskgraph_cache'), N_WORKERS)

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
        [(crop_id, os.path.join(BASE_CROP_DATA_DIR, "%s_yield.asc.zip" % crop_id))
         for crop_id in dep_pol_id_map])

    crop_area_path_id_map = dict(
        [(crop_id, os.path.join(BASE_CROP_DATA_DIR, "%s_harea.asc.zip" % crop_id))
         for crop_id in dep_pol_id_map])

    target_crop_total_yield_path_id_map = dict(
        [(crop_id, os.path.join(TARGET_CROP_FILE_DIR, "%s_tot_yield.tif" % crop_id))
         for crop_id in dep_pol_id_map])

    crop_yield_task_id_map = {}
    crop_area_task_id_map = {}
    pollinator_dependent_micronutrient_yield_path_id_map = {}
    total_micronutrient_yield_path_id_map = {}
    for crop_id in dep_pol_id_map:
        base_crop_yield_path = crop_yield_path_id_map[crop_id]
        target_crop_yield_path = os.path.join(
            TARGET_CROP_FILE_DIR, os.path.basename(
                base_crop_yield_path).split('.')[0] + '.tif')
        crop_yield_path_id_map[crop_id] = target_crop_yield_path
        crop_yield_task_id_map[crop_id] = task_graph.add_task(
            func=asc_to_tiff,
            args=(base_crop_yield_path, target_crop_yield_path),
            target_path_list=[target_crop_yield_path],
            task_name='asc_to_tiff')

        base_crop_area_path = crop_area_path_id_map[crop_id]
        target_crop_area_path = os.path.join(
            TARGET_CROP_FILE_DIR, os.path.basename(
                base_crop_area_path).split('.')[0] + '.tif')
        crop_area_path_id_map[crop_id] = target_crop_area_path
        crop_area_task_id_map[crop_id] = task_graph.add_task(
            func=asc_to_tiff,
            args=(base_crop_area_path, target_crop_area_path),
            target_path_list=[target_crop_area_path],
            task_name='asc_to_tiff')

        total_yield_task = task_graph.add_task(
            func=MultRastersAndScalar(
                [target_crop_area_path, target_crop_yield_path],
                1.0-proportion_refuse_crop_map[crop_id],
                target_crop_total_yield_path_id_map[crop_id]),
            target_path_list=[target_crop_total_yield_path_id_map[crop_id]],
            dependent_task_list=[
                crop_area_task_id_map[crop_id],
                crop_yield_task_id_map[crop_id]],
            task_name='MultRastersAndScalar'
            )

        for micronutrient_id in MICRONUTRIENT_LIST:
            micronutrient_yield_path = os.path.splitext(
                target_crop_total_yield_path_id_map[crop_id])[0] + (
                    '_%s.tif' % (micronutrient_id))
            total_micronutrient_yield_path_id_map[
                (micronutrient_id, crop_id)] = micronutrient_yield_path

            # the 1e4 converts the Mg to g for nutrient units
            micronutrient_task = task_graph.add_task(
                func=MultRastersAndScalar(
                    [target_crop_total_yield_path_id_map[crop_id]],
                    1e4 * float(crop_table[crop_table['crop'] == crop_id][
                        micronutrient_id]),
                    micronutrient_yield_path),
                target_path_list=[micronutrient_yield_path],
                dependent_task_list=[total_yield_task],
                task_name='MultRastersAndScalar'
                )

            pollinator_dependent_micronutrient_yield_path = os.path.splitext(
                micronutrient_yield_path)[0] + '_%s.tif' % 'pol_dep'

            pollinator_dependent_micronutrient_yield_path_id_map[
                (micronutrient_id, crop_id)] = (
                pollinator_dependent_micronutrient_yield_path)
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

    for (c_type, period), crop_id_list in (
            crop_id_functional_type_map.iteritems()):
        for micronutrient_id in MICRONUTRIENT_LIST:
            # calculate both the total and pollinator dependent micronutrient
            pol_dep_micronutrient_functional_yield_path = os.path.join(
                WORKSPACE_DIR, 'pol_dep_%s_%s_%s_yield.tif' % (
                    micronutrient_id, c_type, period))

            pollinator_functional_crop_path_list = [
                pollinator_dependent_micronutrient_yield_path_id_map[
                    (micronutrient_id, crop_id)]
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
                task_name='raster_calculator_sum_pol_dep_micronutrient'
            )

            total_micronutrient_functional_yield_path = os.path.join(
                WORKSPACE_DIR, 'total_%s_%s_%s_yield.tif' % (
                    micronutrient_id, c_type, period))

            total_functional_crop_path_list = [
                total_micronutrient_yield_path_id_map[(micronutrient_id, crop_id)]
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
                task_name='raster_calculator_sum_total_micronutrient'
            )

            # At 1 km LUH2 data, subtract pollinator dependent from total micronutrient production for any grid cells < 0.30; keep it at total for micronutrient production for any grid cells >0.30

            # we need to align the lh2 and yield rasters so they're the
            # same overlap and pixel size
            luh2_raster_path = RASTER_FUNCTIONAL_TYPE_MAP[(c_type, period)]
            luh2_raster_info = pygeoprocessing.get_raster_info(
                luh2_raster_path)

            # put aligned rasters in a subdirectory
            target_aligned_rasters = [
                os.path.join(
                    WORKSPACE_DIR, 'aligned_%s_%s_%s_rasters' % (
                        c_type, period, micronutrient_id),
                    os.path.basename(x)) for x in [
                        luh2_raster_path,
                        total_micronutrient_functional_yield_path,
                        pol_dep_micronutrient_functional_yield_path]]

            # clip and align the yield and LUH2 rasters to be the intersection
            # of bounding boxes and the pixel size of LUH2
            align_luh2_task = task_graph.add_task(
                func=pygeoprocessing.align_and_resize_raster_stack,
                args=(
                    [luh2_raster_path,
                     total_micronutrient_functional_yield_path,
                     pol_dep_micronutrient_functional_yield_path],
                    target_aligned_rasters,
                    ['nearest']*3, luh2_raster_info['pixel_size'],
                    'intersection'),
                kwargs={
                    'raster_align_index': 0,
                    'gtiff_creation_options': GTIFF_CREATION_OPTIONS},
                target_path_list=target_aligned_rasters,
                dependent_task_list=[pol_dep_micro_task, total_micro_task],
                task_name='align_resize_raster_%s' % str(
                    (c_type, period, micronutrient_id)))

            pol_dep_risk_path = os.path.join(
                WORKSPACE_DIR, 'LUH2_risk_%s_%s_%s_yield.tif' % (
                    micronutrient_id, c_type, period))

            # the 0.3 comes from Becky saying that's where we should threshold
            pol_dep_risk_task = task_graph.add_task(
                func=MicroNutrientRiskThreshold(
                    target_aligned_rasters, 0.3, pol_dep_risk_path),
                target_path_list=[pol_dep_risk_path],
                dependent_task_list=[align_luh2_task],
                task_name='MicroNutrientRiskThreshold_%s' % str(
                    (micronutrient_id, c_type, period)))

    LOGGER.info("joining taskgraph")
    task_graph.join()
    return

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
        dependent_task_list=[sumtask_id_map['ag']],
        task_name='raster_calculator'
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

    task_graph.close()
    task_graph.join()


if __name__ == '__main__':
    main()
