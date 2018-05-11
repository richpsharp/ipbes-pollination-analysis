"""
Pollination analysis for IPBES.

    From "IPBES Methods: Pollination Contribution to Human Nutrition."
"""
import re
import logging
import functools

import pygeoprocessing
import google.cloud.client
import google.cloud.storage
from osgeo import gdal
from osgeo import osr
import pandas
import numpy
import scipy.ndimage.morphology
import reproduce
import taskgraph

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'))
LOGGER = logging.getLogger('ipbes_pollination')

# the following are the globio landcover codes. A tuple (x, y) indicates the
# inclusive range from x to y.
GLOBIO_AG_CODES = [2, (230, 232)]
GLOBIO_NATURAL_CODES = [6, (50, 180)]

WORKING_DIR = '.'
GOOGLE_BUCKET_KEY_PATH = "ecoshard-202992-key.json"
NODATA = -9999
N_WORKERS = 2


def main():
    """Entry point."""
    reproduce_env = reproduce.Reproduce(WORKING_DIR)
    LOGGER.debug(reproduce_env['CACHE_DIR'])
    task_graph = taskgraph.TaskGraph(
        reproduce_env['CACHE_DIR'], N_WORKERS)

    # The following table is used for:
    #  Crop pollination dependency was determined for 115 crops. Dependency
    #   was defined as the percent by which yields are reduced for each crop
    #   with inadequate pollination (ranging from 0-95%), according to
    #   Klein et al. (2007).
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
    reproduce_env.register_data(
        'crop_nutrient_table',
        'pollination_data/crop_nutrient.csv',
        expected_hash=('md5', '2fbe7455357f8008a12827fd88816fc1'),
        constructor_func=google_bucket_fetcher(
            'https://storage.googleapis.com/ecoshard-root/'
            'crop_nutrient_md5_2fbe7455357f8008a12827fd88816fc1.csv',
            GOOGLE_BUCKET_KEY_PATH))

    # The following table is used for:
    # Pollinator habitat was defined as any natural land covers, as defined
    #  in Table X  (GLOBIO land-cover classes 6, secondary vegetation, and
    #  50-180, various types of primary vegetation). To test sensitivity to
    #  this definition we included "semi-natural" habitats (GLOBIO land-cover
    #  classes 3, 4, and 5; pasture, rangeland and forestry, respectively) in
    #  addition to "natural", and repeated all analyses with semi-natural
    #  plus natural habitats, but this did not substantially alter the results
    #  so we do not include it in our final analysis or code base.
    reproduce_env.register_data(
        'globio_class_table',
        'pollination_data/GLOBIOluclass.csv',
        expected_hash=('md5', '4506b5c87fe70f7fba63eb4ee5b1e2d0'),
        constructor_func=google_bucket_fetcher(
            'https://storage.cloud.google.com/ecoshard-root/'
            'GLOBIOluclass_md5_4506b5c87fe70f7fba63eb4ee5b1e2d0.csv',
            GOOGLE_BUCKET_KEY_PATH))

    # The proportional area of natural within 2 km was calculated for every
    #  pixel of agricultural land (GLOBIO land-cover classes 2, 230, 231, and
    #  232) at 10 arc seconds (~300 m) resolution. This 2 km scale represents
    #  the distance most commonly found to be predictive of pollination
    #  services (Kennedy et al. 2013).
    reproduce_env.register_data(
        'kernel_raster',
        'kernel_raster.tif',
        constructor_func=(
            lambda kernel_filepath: create_radial_convolution_mask(
                0.00277778, 2000., kernel_filepath)))

    crop_nutrient_df = pandas.read_csv(reproduce_env['crop_nutrient_table'])
    globio_df = pandas.read_csv(reproduce_env['globio_class_table'])

    landcover_data = {
        'GLOBIO4_LU_10sec_2050_SSP5_RCP85': (
            "https://storage.cloud.google.com/ecoshard-root/globio_landcover/GLOBIO4_LU_10sec_2050_SSP5_RCP85_md5_1b3cc1ce6d0ff14d66da676ef194f130.tif",
            ('md5', '1b3cc1ce6d0ff14d66da676ef194f130')),
        'GLOBIO4_LU_10sec_2050_SSP1_RCP26': (
            "https://storage.cloud.google.com/ecoshard-root/globio_landcover/GLOBIO4_LU_10sec_2050_SSP1_RCP26_md5_803166420f51e5ef7dcaa970faa98173.tif",
            ('md5', '803166420f51e5ef7dcaa970faa98173')),
        'GLOBIO4_LU_10sec_2050_SSP3_RCP70': (
            "https://storage.cloud.google.com/ecoshard-root/globio_landcover/GLOBIO4_LU_10sec_2050_SSP3_RCP70_md5_e77077a3220a36f7f0441bbd0f7f14ab.tif",
            ('md5', 'e77077a3220a36f7f0441bbd0f7f14ab')),
        'Globio4_landuse_10sec_1850': (
            "https://storage.cloud.google.com/ecoshard-root/globio_landcover/Globio4_landuse_10sec_1850_md5_0b7fcb4b180d46b4fc2245beee76d6b9.tif",
            ('md5', '0b7fcb4b180d46b4fc2245beee76d6b9')),
        'Globio4_landuse_10sec_2015': (
            "https://storage.cloud.google.com/ecoshard-root/globio_landcover/Globio4_landuse_10sec_2015_md5_939a57c2437cd09bd5a9eb472b9bd781.tif",
            ('md5', '939a57c2437cd09bd5a9eb472b9bd781')),
        'Globio4_landuse_10sec_1980': (
            "https://storage.cloud.google.com/ecoshard-root/globio_landcover/Globio4_landuse_10sec_1980_md5_f6384eac7579318524439df9530ca1f4.tif",
            ('md5', 'f6384eac7579318524439df9530ca1f4')),
        'Globio4_landuse_10sec_1945': (
            "https://storage.cloud.google.com/ecoshard-root/globio_landcover/Globio4_landuse_10sec_1945_md5_52c7b4c38c26defefa61132fd25c5584.tif",
            ('md5', '52c7b4c38c26defefa61132fd25c5584')),
        'Globio4_landuse_10sec_1910': (
            "https://storage.cloud.google.com/ecoshard-root/globio_landcover/Globio4_landuse_10sec_1910_md5_e7da8fa29db305ff63c99fed7ca8d5e2.tif",
            ('md5', 'e7da8fa29db305ff63c99fed7ca8d5e2')),
        'Globio4_landuse_10sec_1900': (
            "https://storage.cloud.google.com/ecoshard-root/globio_landcover/Globio4_landuse_10sec_1900_md5_f5db818a5b16799bf2cb627e574120a4.tif",
            ('md5', 'f5db818a5b16799bf2cb627e574120a4')),
        }

    # mask landcover into agriculture and pollinator habitat
    for landcover_key, (landcover_url, expected_hash) in landcover_data.iteritems():
        landcover_local_path = 'landcover/%s.tif' % landcover_key
        landcover_fetch_task = task_graph.add_task(
            func=register_data,
            args=(
                reproduce_env,
                landcover_key,
                landcover_local_path,
                expected_hash,
                functools.partial(
                    google_bucket_fetcher,
                    landcover_url, GOOGLE_BUCKET_KEY_PATH)))

        hab_task_path_list = []

        for mask_prefix, globio_codes in [
                ('ag', GLOBIO_AG_CODES), ('hab', GLOBIO_NATURAL_CODES)]:
            mask_key = '%s_%s_mask' % (landcover_key, mask_prefix)
            local_mask_path = '%s_mask/%s.tif' % (mask_prefix, mask_key)
            mask_target_path = reproduce_env.predict_path(local_mask_path)

            mask_task = task_graph.add_task(
                func=mask_raster,
                args=(
                    reproduce_env.predict_path(landcover_local_path),
                    globio_codes, mask_target_path),
                target_path_list=[mask_target_path],
                dependent_task_list=[landcover_fetch_task])

            mask_register_task = task_graph.add_task(
                func=register_data,
                args=(
                    reproduce_env,
                    mask_key,
                    local_mask_path,
                    None,
                    None),
                dependent_task_list=[mask_task])

            if mask_prefix == 'hab':
                hab_task_path_list.append(
                    (mask_register_task, local_mask_path, landcover_key))

        for mask_task, local_mask_path, landcover_key in hab_task_path_list:
            proportional_hab_area_2km_key = (
                '%s_hab_prop_area_2km' % landcover_key)
            local_proportional_hab_area_2km_path = (
                'proportional_area/%s.tif' % proportional_hab_area_2km_key)
            task_graph.add_task(
                func=pygeoprocessing.convolve_2d,
                args=(
                    (reproduce_env.predict_path(local_mask_path), 1),
                    (reproduce_env['kernel_raster'], 1),
                    reproduce_env.predict_path(
                        local_proportional_hab_area_2km_path)),
                kwargs={'working_dir': reproduce_env['CACHE_DIR']},
                dependent_task_list=[mask_task])

    task_graph.close()
    task_graph.join()

    """
    pygeoprocessing.convolve_2d(
        signal_path_band, kernel_path_band, target_path,
        ignore_nodata=False, mask_nodata=True, normalize_kernel=False,
        target_datatype=gdal.GDT_Float64,
        target_nodata=None,
        gtiff_creation_options=_DEFAULT_GTIFF_CREATION_OPTIONS,
        working_dir=None)
    """

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
    # TODO: threshold


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
    degree_len_0 = 110574 # length at 0 degrees
    degree_len_60 = 111412 # length at 60 degrees
    pixel_size_m = pixel_size_degree * (degree_len_0 + degree_len_60) / 2.0
    pixel_radius = numpy.ceil(radius_meters / pixel_size_m)
    n_pixels = (int(pixel_radius) * 2 + 1)
    sample_pixels = 200
    mask = numpy.ones((sample_pixels * n_pixels, sample_pixels * n_pixels))
    mask[mask.shape[0]/2, mask.shape[0]/2] = 0
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
    kernel_band.WriteArray(kernel_array)


def google_bucket_fetcher(url, json_key_path):
    """Create a function to download a Google Blob to a given path.

    Parameters:
        url (string): url to blob, matches the form
            '^https://storage.cloud.google.com/([^/]*)/(.*)$'
        json_key_path (string): path to Google Cloud private key generated
            by https://cloud.google.com/iam/docs/creating-managing-service-account-keys

    Returns:
        a function with a single `path` argument to the target file. Invoking
            this function will download the Blob to `path`.

    """
    def _google_bucket_fetcher(path):
        """Fetch blob `url` to `path`."""
        url_matcher = re.match(
            '^https://[^/]*\.com/([^/]*)/(.*)$', url)
        LOGGER.debug(url)
        client = google.cloud.storage.client.Client.from_service_account_json(
            json_key_path)
        bucket_id = url_matcher.group(1)
        LOGGER.debug('parsing bucket %s from %s', bucket_id, url)
        bucket = client.get_bucket(bucket_id)
        blob_id = url_matcher.group(2)
        LOGGER.debug('loading blob %s from %s', blob_id, url)
        blob = google.cloud.storage.Blob(blob_id, bucket)
        LOGGER.info('downloading blob %s from %s' % (path, url))
        blob.download_to_filename(path)
    return _google_bucket_fetcher


class MaskCodes(object):
    def __init__(self, code_array):
        self.code_array = code_array

    def __call__(self, base_array):
        return numpy.isin(base_array, self.code_array)


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
    LOGGER.debug('expanded code array %s', code_array)

    pygeoprocessing.raster_calculator(
        [(base_path, 1)], MaskCodes(code_array), target_path, gdal.GDT_Byte, 2)


def register_data(
        reproduce_env, data_key, local_path, expected_hash, constructor_func):
    """Wrap a call to `reproduce_env.register_data` for pickling."""
    reproduce_env.register_data(
        data_key,
        local_path,
        expected_hash=expected_hash,
        constructor_func=constructor_func)


if __name__ == '__main__':
    main()
