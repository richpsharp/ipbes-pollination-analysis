"""
Pollination analysis for IPBES.

    From "IPBES Methods: Pollination Contribution to Human Nutrition."
    https://www.dropbox.com/s/gc4b1miw2zypuke/IPBES%20Methods_Pollination_RS.docx?dl=0
"""
import sys
import zipfile
import time
import os
import re
import logging

import google.cloud.client
import google.cloud.storage
from osgeo import gdal
from osgeo import osr
import pandas
import numpy
import scipy.ndimage.morphology
import reproduce
import taskgraph
import pygeoprocessing

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger('ipbes_pollination')

# the following are the globio landcover codes. A tuple (x, y) indicates the
# inclusive range from x to y.
GLOBIO_AG_CODES = [2, (230, 232)]
GLOBIO_NATURAL_CODES = [6, (50, 180)]

WORKING_DIR = 'workspace'
GOOGLE_BUCKET_KEY_PATH = "ecoshard-202992-key.json"
NODATA = -9999
N_WORKERS = 4


def main():
    """Entry point."""
    reproduce_env = reproduce.Reproduce(WORKING_DIR)
    LOGGER.debug(reproduce_env['CACHE_DIR'])
    task_graph = taskgraph.TaskGraph(
        reproduce_env['CACHE_DIR'], N_WORKERS, delayed_start=True,
        reporting_interval=5.0)

    # The following table is used for:
    # Pollinator habitat was defined as any natural land covers, as defined
    #  in Table X  (GLOBIO land-cover classes 6, secondary vegetation, and
    #  50-180, various types of primary vegetation). To test sensitivity to
    #  this definition we included "semi-natural" habitats (GLOBIO land-cover
    #  classes 3, 4, and 5; pasture, rangeland and forestry, respectively) in
    #  addition to "natural", and repeated all analyses with semi-natural
    #  plus natural habitats, but this did not substantially alter the results
    #  so we do not include it in our final analysis or code base.
    globio_class_table_url = (
        'https://storage.cloud.google.com/ecoshard-root/'
        'GLOBIOluclass_md5_4506b5c87fe70f7fba63eb4ee5b1e2d0.csv')
    globio_table_path = os.path.join(
        reproduce_env['DATA_DIR'], 'pollination_data',
        os.path.basename(globio_class_table_url))
    globio_class_table_fetch_task = task_graph.add_task(
        func=google_bucket_fetch_and_validate,
        args=(
            globio_class_table_url, GOOGLE_BUCKET_KEY_PATH,
            globio_table_path),
        target_path_list=[globio_table_path],
        task_name=f'fetch {os.path.basename(globio_class_table_url)}')

    # The proportional area of natural within 2 km was calculated for every
    #  pixel of agricultural land (GLOBIO land-cover classes 2, 230, 231, and
    #  232) at 10 arc seconds (~300 m) resolution. This 2 km scale represents
    #  the distance most commonly found to be predictive of pollination
    #  services (Kennedy et al. 2013).
    kernel_raster_path = os.path.join(
        reproduce_env['DATA_DIR'], 'radial_kernel.tif')
    kernel_task = task_graph.add_task(
        func=create_radial_convolution_mask,
        args=(0.00277778, 2000., kernel_raster_path),
        target_path_list=[kernel_raster_path],
        task_name='make convolution kernel')

    landcover_data = {
        'GLOBIO4_LU_10sec_2050_SSP5_RCP85': "https://storage.cloud.google.com/ecoshard-root/globio_landcover/GLOBIO4_LU_10sec_2050_SSP5_RCP85_md5_1b3cc1ce6d0ff14d66da676ef194f130.tif",
        'GLOBIO4_LU_10sec_2050_SSP1_RCP26': "https://storage.cloud.google.com/ecoshard-root/globio_landcover/GLOBIO4_LU_10sec_2050_SSP1_RCP26_md5_803166420f51e5ef7dcaa970faa98173.tif",
        'GLOBIO4_LU_10sec_2050_SSP3_RCP70': "https://storage.cloud.google.com/ecoshard-root/globio_landcover/GLOBIO4_LU_10sec_2050_SSP3_RCP70_md5_e77077a3220a36f7f0441bbd0f7f14ab.tif",
        'Globio4_landuse_10sec_1850': "https://storage.cloud.google.com/ecoshard-root/globio_landcover/Globio4_landuse_10sec_1850_md5_0b7fcb4b180d46b4fc2245beee76d6b9.tif",
        'Globio4_landuse_10sec_2015': "https://storage.cloud.google.com/ecoshard-root/globio_landcover/Globio4_landuse_10sec_2015_md5_939a57c2437cd09bd5a9eb472b9bd781.tif",
        'Globio4_landuse_10sec_1980': "https://storage.cloud.google.com/ecoshard-root/globio_landcover/Globio4_landuse_10sec_1980_md5_f6384eac7579318524439df9530ca1f4.tif",
        'Globio4_landuse_10sec_1945': "https://storage.cloud.google.com/ecoshard-root/globio_landcover/Globio4_landuse_10sec_1945_md5_52c7b4c38c26defefa61132fd25c5584.tif",
        'Globio4_landuse_10sec_1910': "https://storage.cloud.google.com/ecoshard-root/globio_landcover/Globio4_landuse_10sec_1910_md5_e7da8fa29db305ff63c99fed7ca8d5e2.tif",
        'Globio4_landuse_10sec_1900': "https://storage.cloud.google.com/ecoshard-root/globio_landcover/Globio4_landuse_10sec_1900_md5_f5db818a5b16799bf2cb627e574120a4.tif",
        }

    # mask landcover into agriculture and pollinator habitat
    for landcover_key, landcover_url in landcover_data.items():
        landcover_env_path = f'landcover/{os.path.basename(landcover_url)}'
        landcover_path = reproduce_env.predict_path(landcover_env_path)
        landcover_fetch_task = task_graph.add_task(
            func=google_bucket_fetch_and_validate,
            args=(landcover_url, GOOGLE_BUCKET_KEY_PATH, landcover_path),
            target_path_list=[landcover_path],
            task_name=f'fetch {landcover_key}')

        _ = task_graph.add_task(
            func=build_overviews,
            args=(landcover_path, 'mode'),
            target_path_list=[f'{landcover_path}.ovr'],
            dependent_task_list=[landcover_fetch_task],
            task_name=f'compress {os.path.basename(landcover_path)}')

        hab_task_path_list = []

        for mask_prefix, globio_codes in [
                ('ag', GLOBIO_AG_CODES), ('hab', GLOBIO_NATURAL_CODES)]:
            mask_key = f'{landcover_key}_{mask_prefix}_mask'
            mask_target_path = os.path.join(
                reproduce_env['DATA_DIR'],
                f'{mask_prefix}_mask/{mask_key}.tif')

            mask_task = task_graph.add_task(
                func=mask_raster,
                args=(
                    reproduce_env.predict_path(landcover_env_path),
                    globio_codes, mask_target_path),
                target_path_list=[mask_target_path],
                dependent_task_list=[landcover_fetch_task],
                task_name=f'mask {mask_key}',)

            _ = task_graph.add_task(
                func=build_overviews,
                args=(mask_target_path, 'mode'),
                target_path_list=[
                    f'{mask_target_path}.ovr'],
                dependent_task_list=[mask_task],
                task_name=f'compress {os.path.basename(mask_target_path)}')

            if mask_prefix == 'hab':
                hab_task_path_list.append(
                    (mask_task, mask_target_path))

        raster_tasks_to_threshold_list = []
        for mask_task, mask_path, in hab_task_path_list:
            proportional_hab_area_2km_path = os.path.join(
                reproduce_env['DATA_DIR'], 'proportional_area',
                f'{landcover_key}_hab_prop_area_2km.tif')
            convolve2d_task = task_graph.add_task(
                func=pygeoprocessing.convolve_2d,
                args=[
                    (mask_path, 1), (kernel_raster_path, 1),
                    proportional_hab_area_2km_path],
                kwargs={
                    'working_dir': reproduce_env['CACHE_DIR'],
                    'gtiff_creation_options': (
                        'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
                        'PREDICTOR=3', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256',
                        'NUM_THREADS=ALL_CPUS')},
                dependent_task_list=[mask_task, kernel_task],
                target_path_list=[proportional_hab_area_2km_path],
                task_name=(
                    'calculate proportional'
                    f' {os.path.basename(proportional_hab_area_2km_path)}'))

            _ = task_graph.add_task(
                func=build_overviews,
                args=(proportional_hab_area_2km_path, 'average'),
                target_path_list=[
                    f'{proportional_hab_area_2km_path}.ovr'],
                dependent_task_list=[convolve2d_task],
                task_name=f'compress {os.path.basename(proportional_hab_area_2km_path)}',
                )

            raster_tasks_to_threshold_list.append(
                (convolve2d_task, proportional_hab_area_2km_path))

        #  1.1.4.  Sufficiency threshold
        #  A threshold of 0.3 was set to evaluate whether there was sufficient
        #   pollinator habitat in the 2 km around farmland to provide pollination
        #   services, based on Kremen et al.'s (2005)  estimate of the area
        #   requirements for achieving full pollination. This produced a map of
        #   wild pollination sufficiency where every agricultural pixel was
        #   designated in a binary fashion: 0 if proportional area of habitat was
        #   less than 0.3; 1 if greater than 0.3. Maps of pollination sufficiency
        #   can be found at (permanent link to output), outputs "poll_suff_..."
        #   below.
        threshold_val = 0.3
        for convolve2d_task, proportional_hab_area_2km_path in (
                raster_tasks_to_threshold_list):
            thresholded_path = os.path.join(
                reproduce_env['DATA_DIR'],
                'thresholded_'
                f'{os.path.basename(proportional_hab_area_2km_path)}')
            threshold_task = task_graph.add_task(
                func=threshold_raster,
                args=(
                    proportional_hab_area_2km_path, threshold_val,
                    thresholded_path),
                target_path_list=[thresholded_path],
                dependent_task_list=[convolve2d_task],
                task_name=f'threshold {os.path.basename(thresholded_path)}')

            _ = task_graph.add_task(
                func=build_overviews,
                args=(thresholded_path, 'mode'),
                target_path_list=[
                    f'{thresholded_path}.ovr'],
                dependent_task_list=[threshold_task],
                task_name=f'compress {os.path.basename(thresholded_path)}',
                )

    # 1.2.    POLLINATION-DEPENDENT NUTRIENT PRODUCTION
    # Pollination-dependence of crops, crop yields, and crop micronutrient
    # content were combined in an analysis to calculate pollination-dependent
    # nutrient production, following Chaplin-Kramer et al. (2012).
    # 1.2.2.  Pollination dependency
    # Crop pollination dependency was determined for 115 crops (permanent link
    # to table). Dependency was defined as the percent by which yields are
    # reduced for each crop with inadequate pollination (ranging from 0-95%),
    # according to Klein et al. (2007).

    # Crop content of critical macro and micronutrients (KJ energy/100 g, IU
    #   Vitamin A/ 100 g and mcg Folate/100g) for the 115 crops were taken
    #   from USDA (2011) . The USDA (2011) data also provided estimated refuse
    #   of the food item (e.g., peels, seeds). The pollination-dependent yield
    #   was reduced by this refuse percent and then multiplied by nutrient
    #   content, and summed across all crops to derive pollination-dependent
    #   nutrient yields (KJ/ha, IU Vitamin A/ha, mcg Folate/ha) for each
    #   nutrient at 5 arc min. The full table used in this analysis can be
    # found at https://storage.googleapis.com/ecoshard-root/'
    # 'crop_nutrient_md5_d6e67fd79ef95ab2dd44ca3432e9bb4d.csv
    crop_nutrient_url = (
        'https://storage.googleapis.com/ecoshard-root/'
        'crop_nutrient_md5_2fbe7455357f8008a12827fd88816fc1.csv')
    crop_nutrient_table_path = os.path.join(
        reproduce_env['DATA_DIR'], 'pollination_data',
        os.path.basename(crop_nutrient_url))

    crop_nutrient_table_fetch_task = task_graph.add_task(
        func=google_bucket_fetch_and_validate,
        args=(
            crop_nutrient_url, GOOGLE_BUCKET_KEY_PATH,
            crop_nutrient_table_path),
        target_path_list=[crop_nutrient_table_path],
        task_name=f'fetch {os.path.basename(crop_nutrient_table_path)}')

    # 1.2.3.  Crop production
    # Spatially-explicit global crop yields (tons/ha) at 5 arc min (~10 km) were
    # taken from Monfreda et al. (2008) for 115 crops (permanent link to crop
    # yield folder). These yields were multiplied by crop pollination dependency
    # to calculate the pollination-dependent crop yield for each 5 min grid
    # cell.
    # Note the monfredia maps are in units of per-hectare yields
    yield_zip_url = (
        'https://storage.cloud.google.com/ecoshard-root/ipbes/'
        'monfreda_2008_observed_yield_md5_54c6b8e564973739ba75c8e54ac6f051.zip')
    yield_zip_path = os.path.join(
        reproduce_env['DATA_DIR'], 'observed_yields',
        os.path.basename(yield_zip_url))

    yield_zip_fetch_task = task_graph.add_task(
        func=google_bucket_fetch_and_validate,
        args=(
            yield_zip_url, GOOGLE_BUCKET_KEY_PATH,
            yield_zip_path),
        target_path_list=[yield_zip_path],
        priority=100,
        task_name=f'fetch {os.path.basename(yield_zip_path)}')

    zip_touch_file_path = os.path.join(
        os.path.dirname(yield_zip_path), 'monfreda_2008_observed_yield.txt')
    unzip_yield_task = task_graph.add_task(
        func=unzip_file,
        args=(
            yield_zip_path, os.path.dirname(yield_zip_path),
            zip_touch_file_path),
        target_path_list=[zip_touch_file_path],
        dependent_task_list=[yield_zip_fetch_task],
        task_name=f'unzip monfreda_2008_observed_yield',
        priority=100)

    yield_raster_dir = os.path.join(
        os.path.dirname(yield_zip_path), 'monfreda_2008_observed_yield')

    for nutrient_id, nutrient_name in [
            ('en', 'Energy'), ('va', 'VitA'), ('fo', 'Folate')]:
        tot_prod_nut_1d_path = os.path.join(
            WORKING_DIR, 'tot_prod_nut_prod_rasters',
            f'tot_prod_{nutrient_id}_1d.tif')

        tot_prod_nutrient_task = task_graph.add_task(
            func=create_prod_nutrient_1d_raster,
            args=(
                crop_nutrient_table_path, nutrient_name, yield_raster_dir,
                False, tot_prod_nut_1d_path),
            target_path_list=[tot_prod_nut_1d_path],
            dependent_task_list=[
                unzip_yield_task, crop_nutrient_table_fetch_task],
            task_name=f'total nut prod {nutrient_name}',
            priority=100)

        poll_serv_nutrient_1d_path = os.path.join(
            WORKING_DIR, 'poll_serv_nutrient_rasters',
            f'poll_serv_{nutrient_id}_1d.tif')
        poll_serv_nutrient_task = task_graph.add_task(
            func=create_prod_nutrient_1d_raster,
            args=(
                crop_nutrient_table_path, nutrient_name, yield_raster_dir,
                True, poll_serv_nutrient_1d_path),
            target_path_list=[poll_serv_nutrient_1d_path],
            dependent_task_list=[unzip_yield_task],
            task_name=f'total pol serv {nutrient_name}',
            priority=100)

    # 1.3.    NUTRITION PROVIDED BY WILD POLLINATORS
    # 1.3.1.  Overview
    #   Nutrition provided by wild pollinators on each pixel of agricultural
    # land was calculated according to pollination habitat sufficiency and the
    # pollination-dependent nutrient yields.

    # 1.3.2.  Nutrition production by wild pollinators
    #  Pollinator-dependent nutrient yields at 5 arc minutes were applied to
    # agricultural pixels in the GLOBIO land-use map at 10 arc seconds, and
    # multiplied by wild pollination sufficiency (1 if sufficient, 0 if not) to
    # report pollination-derived nutrient yields in each pixel. We call this
    # "pollinator-derived" instead of "pollination-dependent" because
    # "dependence" is the proportion that yields are reduced if not adequately
    # pollinated. Whether that dependent yield is actually produced is
    # determined by whether there is sufficient natural habitat around the
    # agricultural pixel in question to provide wild pollinators and hence
    # adequate pollination to crops.  These pollinator-derived nutrient yields
    # were then multiplied by the area of the pixel to convert yields to
    # pixel-level production for each nutrient.  Maps of pollination-derived
    # nutrient production for each nutrient in each scenario can be found at
    # (permanent link to output), outputs "poll_serv_..." below.

    task_graph.close()
    task_graph.join()


def create_radial_convolution_mask(
        pixel_size_degree, radius_meters, kernel_filepath):
    """Creates a radial mask to sample pixels in convolution filter.

    Parameters:
        pixel_size_degree (float): size of pixel in degrees.
        radius_meters (float): desired size of radial mask in meters.

    Returns:
        A 2D numpy array that can be used in a convolution to aggregate a
        raster while accounting for partial coverage of the circle on the
        edges of the pixel.
    """
    degree_len_0 = 110574  # length at 0 degrees
    degree_len_60 = 111412  # length at 60 degrees
    pixel_size_m = pixel_size_degree * (degree_len_0 + degree_len_60) / 2.0
    pixel_radius = numpy.ceil(radius_meters / pixel_size_m)
    n_pixels = (int(pixel_radius) * 2 + 1)
    sample_pixels = 200
    mask = numpy.ones((sample_pixels * n_pixels, sample_pixels * n_pixels))
    mask[mask.shape[0]//2, mask.shape[0]//2] = 0
    distance_transform = scipy.ndimage.morphology.distance_transform_edt(mask)
    mask = None
    stratified_distance = distance_transform * pixel_size_m / sample_pixels
    distance_transform = None
    in_circle = numpy.where(stratified_distance <= 2000.0, 1.0, 0.0)
    stratified_distance = None
    reshaped = in_circle.reshape(
        in_circle.shape[0] // sample_pixels, sample_pixels,
        in_circle.shape[1] // sample_pixels, sample_pixels)
    kernel_array = numpy.sum(reshaped, axis=(1, 3)) / sample_pixels**2
    normalized_kernel_array = kernel_array / numpy.sum(kernel_array)
    reshaped = None

    driver = gdal.GetDriverByName('GTiff')
    kernel_raster = driver.Create(
        kernel_filepath.encode('utf-8'), n_pixels, n_pixels, 1,
        gdal.GDT_Float32, options=[
            'BIGTIFF=IF_SAFER', 'TILED=YES', 'BLOCKXSIZE=256',
            'BLOCKYSIZE=256'])

    # Make some kind of geotransform, it doesn't matter what but
    # will make GIS libraries behave better if it's all defined
    kernel_raster.SetGeoTransform([-180, 1, 0, 90, 0, -1])
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    kernel_raster.SetProjection(srs.ExportToWkt())
    kernel_band = kernel_raster.GetRasterBand(1)
    kernel_band.SetNoDataValue(NODATA)
    kernel_band.WriteArray(normalized_kernel_array)


def google_bucket_fetch_and_validate(url, json_key_path, target_path):
    """Create a function to download a Google Blob to a given path.

    Parameters:
        url (string): url to blob, matches the form
            '^https://storage.cloud.google.com/([^/]*)/(.*)$'
        json_key_path (string): path to Google Cloud private key generated
            by https://cloud.google.com/iam/docs/creating-managing-service-account-keys
        target_path (string): path to target file.

    Returns:
        a function with a single `path` argument to the target file. Invoking
            this function will download the Blob to `path`.

    """
    url_matcher = re.match(
        r'^https://[^/]*\.com/([^/]*)/(.*)$', url)
    LOGGER.debug(url)
    client = google.cloud.storage.client.Client.from_service_account_json(
        json_key_path)
    bucket_id = url_matcher.group(1)
    LOGGER.debug(f'parsing bucket {bucket_id} from {url}')
    bucket = client.get_bucket(bucket_id)
    blob_id = url_matcher.group(2)
    LOGGER.debug(f'loading blob {blob_id} from {url}')
    blob = google.cloud.storage.Blob(blob_id, bucket)
    LOGGER.info(f'downloading blob {target_path} from {url}')
    try:
        os.makedirs(os.path.dirname(target_path))
    except os.error:
        pass
    blob.download_to_filename(target_path)
    if not reproduce.valid_hash(target_path, 'embedded'):
        raise ValueError(f"{target_path}' does not match its expected hash")


class MaskCodes(object):
    def __init__(self, code_array):
        self.code_array = code_array

    def __call__(self, base_array):
        return numpy.isin(base_array, self.code_array)


class Threshold(object):
    def __init__(self, threshold_val, nodata_val, target_nodata):
        self.threshold_val = threshold_val
        self.nodata_val = nodata_val
        self.target_nodata = target_nodata

    def __call__(self, base_array):
        result = numpy.empty(base_array.shape, dtype=numpy.uint8)
        result[:] = self.target_nodata
        valid_mask = (base_array != self.nodata_val)
        result[valid_mask] = 0.0
        result[valid_mask & (base_array >= self.threshold_val)] = 1.0
        return result


def threshold_raster(base_path, threshold_val, target_path):
    """Threshold base to be 1 if >= `theshold_val`.

    Parameters:
        base_path (string): path to single band raster.
        threshold_val (numeric): value to use as threshold cutoff
        target_path (string): path to desired output raster, raster is a
            byte type with same dimensions and projection as `base_path`.
            A pixel in this raster will be 1 if the corresponding pixel in
            `base_path` is >= `threshold_val`, 0 in all other cases.

    Returns:
        None.

    """
    base_nodata = pygeoprocessing.get_raster_info(base_path)['nodata'][0]
    target_nodata = 2
    pygeoprocessing.raster_calculator(
        [(base_path, 1)],
        Threshold(threshold_val, base_nodata, target_nodata), target_path,
        gdal.GDT_Byte, target_nodata, gtiff_creation_options=(
            'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
            'PREDICTOR=2', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256',
            'NUM_THREADS=ALL_CPUS'))


def mask_raster(base_path, codes, target_path):
    """Mask `base_path` to 1 where values are in codes. 0 otherwise.

    Parameters:
        base_path (string): path to single band integer raster.
        codes (list): list of integer or tuple integer pairs. Membership in
            `codes` or within the inclusive range of a tuple in `codes`
            is sufficient to mask the corresponding raster integer value
            in `base_path` to 1 for `target_path`.
        target_path (string): path to desired mask raster. Any corresponding
            pixels in `base_path` that also match a value or range in
            `codes` will be masked to 1 in `target_path`. All other values
            are 0.

    Returns:
        None.

    """
    code_list = [
        item for sublist in [
            range(x[0], x[1]+1) if isinstance(x, tuple) else [x]
            for x in codes] for item in sublist]
    code_array = numpy.array(code_list)
    LOGGER.debug(f'expanded code array {code_array}')

    pygeoprocessing.raster_calculator(
        [(base_path, 1)], MaskCodes(code_array), target_path, gdal.GDT_Byte,
        2, gtiff_creation_options=(
            'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
            'PREDICTOR=2', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256',
            'NUM_THREADS=ALL_CPUS'))


def unzip_file(zipfile_path, target_dir, touchfile_path):
    """Unzip contents of `zipfile_path`.

    Parameters:
        zipfile_path (string): path to a zipped file.
        target_dir (string): path to extract zip file to.
        touchfile_path (string): path to a file to create if unzipping is
            successful.

    Returns:
        None.

    """
    with zipfile.ZipFile(zipfile_path, 'r') as zip_ref:
        zip_ref.extractall(target_dir)

    with open(touchfile_path, 'w') as touchfile:
        touchfile.write(f'unzipped {zipfile_path}')


def build_overviews(local_path, resample_method):
    """Builds as many overviews as possible using 'average' interpolation.

    Parameters:
        local_path (string): path to GTiff raster to build overviews on. The
            overview will be built externally in a .ovr file in the same
            directory.
        resample_method (string): interpolation mode for overviews, must be one of
            'nearest', 'average', 'gauss', 'cubic', 'cubicspline', 'lanczos',
            'average_mp', 'average_magphase', 'mode'.

    Returns:
        None.

    """
    min_dimension = min(
        pygeoprocessing.get_raster_info(local_path)['raster_size'])
    LOGGER.info(f"min min_dimension {min_dimension}")
    raster_copy = gdal.Open(local_path)

    overview_levels = []
    current_level = 2
    while True:
        if min_dimension // current_level == 0:
            break
        overview_levels.append(current_level)
        current_level *= 2
    LOGGER.info(f'level list: {overview_levels}')
    gdal.SetConfigOption('COMPRESS_OVERVIEW', 'LZW')
    raster_copy.BuildOverviews(
        resample_method, overview_levels, callback=_make_logger_callback(
            f'build overview for {os.path.basename(local_path)} '
            '%.2f%% complete'))


def _make_logger_callback(message):
    """Build a timed logger callback that prints `message` replaced.

    Parameters:
        message (string): a string that expects a %f placement variable,
            for % complete.

    Returns:
        Function with signature:
            logger_callback(df_complete, psz_message, p_progress_arg)
    """
    def logger_callback(df_complete, psz_message, p_progress_arg):
        """The argument names come from the GDAL API for callbacks."""
        try:
            current_time = time.time()
            if ((current_time - logger_callback.last_time) > 5.0 or
                    (df_complete == 1.0 and
                     logger_callback.total_time >= 5.0)):
                LOGGER.info(message, df_complete * 100)
                logger_callback.last_time = current_time
                logger_callback.total_time += current_time
        except AttributeError:
            logger_callback.last_time = time.time()
            logger_callback.total_time = 0.0

    return logger_callback


def create_prod_nutrient_1d_raster(
        crop_nutrient_df_path, nutrient_name, yield_raster_dir,
        consider_pollination, target_production_path):
    """Create total production yield for a nutrient for all crops.

    Parameters:
        crop_nutrient_df_path (str): path to CSV with at least the
            column `filenm`, `nutrient_name`, `Percent refuse crop`, and
            `Pollination dependence crop`.
        nutrient_name (str): nutrient name to use to index into the crop
            data frame.
        yield_raster_dir (str): path to a directory that has files of the
            format `[crop_name]_yield_map.tif` where `crop_name` is a value
            in the `filenm` column of `crop_nutrient_df`.
        consider_pollination (bool): if True, multiply yields by pollinator
            dependence ratio.
        target_production_path (str): path to target raster that will contain
            a per-pixel amount of pollinator produced `nutrient_name`
            calculated as the sum(
                crop_yield_map * (100-Percent refuse crop) *
                (Pollination dependence crop) * nutrient) * (ha / pixel map))

    Returns:
        None.

    """
    crop_nutrient_df = pandas.read_csv(crop_nutrient_df_path)
    yield_raster_path_band_list = []
    pollination_yield_factor_list = []
    for _, row in crop_nutrient_df.iterrows():
        yield_raster_path = os.path.join(
            yield_raster_dir, f"{row['filenm']}_yield_map.tif")
        if os.path.exists(yield_raster_path):
            yield_raster_path_band_list.append(
                (yield_raster_path, 1))
            pollination_yield_factor_list.append(
                (1. - row['Percent refuse'] / 100.) * row[nutrient_name])
            if consider_pollination:
                pollination_yield_factor_list[-1] *= (
                    row['Pollination dependence'])
        else:
            raise ValueError(f"not found {yield_raster_path}")

    yield_raster_info = pygeoprocessing.get_raster_info(
        yield_raster_path_band_list[0][0])
    y_lat_array = numpy.linspace(
        yield_raster_info['geotransform'][3],
        yield_raster_info['geotransform'][3] +
        yield_raster_info['geotransform'][5] *
        yield_raster_info['raster_size'][1],
        yield_raster_info['raster_size'][1])

    y_ha_array = area_of_pixel(
        abs(yield_raster_info['geotransform'][1]),
        y_lat_array) / 10000.0
    y_ha_column = y_ha_array.reshape((y_ha_array.size, 1))

    yield_nodata = yield_raster_info['nodata'][0]

    def production_op(pixel_ha, *crop_yield_array_list):
        """Calculate total nutrient production."""
        result = numpy.empty_like(crop_yield_array_list[0])
        result[:] = 0.0
        all_valid = numpy.zeros(result.shape, dtype=numpy.bool)

        for crop_index, crop_array in enumerate(crop_yield_array_list):
            valid_mask = crop_array != yield_nodata
            all_valid |= valid_mask
            result[valid_mask] += (
                crop_array[valid_mask] * pixel_ha[valid_mask] *
                pollination_yield_factor_list[crop_index])
        result[~all_valid] = yield_nodata
        return result

    pygeoprocessing.raster_calculator(
        [y_ha_column] + yield_raster_path_band_list, production_op,
        target_production_path, gdal.GDT_Float32, yield_nodata)


def area_of_pixel(pixel_size, center_lat):
    """Calculate m^2 area of a wgs84 square pixel.

    Adapted from: https://gis.stackexchange.com/a/127327/2397

    Parameters:
        pixel_size (float): length of side of pixel in degrees.
        center_lat (float): latitude of the center of the pixel. Note this
            value +/- half the `pixel-size` must not exceed 90/-90 degrees
            latitude or an invalid area will be calculated.

    Returns:
        Area of square pixel of side length `pixel_size` centered at
        `center_lat` in m^2.

    """
    a = 6378137  # meters
    b = 6356752.3142  # meters
    e = numpy.sqrt(1-(b/a)**2)
    area_list = []
    for f in [center_lat+pixel_size/2, center_lat-pixel_size/2]:
        zm = 1 - e*numpy.sin(numpy.radians(f))
        zp = 1 + e*numpy.sin(numpy.radians(f))
        area_list.append(
            numpy.pi * b**2 * (
                numpy.log(zp/zm) / (2*e) +
                numpy.sin(numpy.radians(f)) / (zp*zm)))
    return pixel_size / 360. * (area_list[0]-area_list[1])


if __name__ == '__main__':
    main()
