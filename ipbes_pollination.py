"""
Pollination analysis for IPBES.

    From "IPBES Methods: Pollination Contribution to Human Nutrition."
"""
from osgeo import gdal
from osgeo import osr
import pandas
import numpy
import scipy.ndimage.morphology
import reproduce
import taskgraph


WORKING_DIR = '.'
NODATA = -9999
N_WORKERS = -1


def main():
    """Entry point."""
    reproduce_env = reproduce.Reproduce(WORKING_DIR)
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
    task_graph.add_task(
        func=reproduce_env.register_data, args=(
            'crop_nutrient_table',
            'pollination_data/crop_nutrient.csv'))

    # The following table is used for:
    # Pollinator habitat was defined as any natural land covers, as defined
    #  in Table X  (GLOBIO land-cover classes 6, secondary vegetation, and
    #  50-180, various types of primary vegetation). To test sensitivity to
    #  this definition we included "semi-natural" habitats (GLOBIO land-cover
    #  classes 3, 4, and 5; pasture, rangeland and forestry, respectively) in
    #  addition to "natural", and repeated all analyses with semi-natural
    #  plus natural habitats, but this did not substantially alter the results
    #  so we do not include it in our final analysis or code base.
    task_graph.add_task(
        func=reproduce_env.register_data, args=(
            'globio_class_table',
            'pollination_data/GLOBIOluclass.csv'))


    # The proportional area of natural within 2 km was calculated for every
    #  pixel of agricultural land (GLOBIO land-cover classes 2, 230, 231, and
    #  232) at 10 arc seconds (~300 m) resolution. This 2 km scale represents
    #  the distance most commonly found to be predictive of pollination
    #  services (Kennedy et al. 2013).
    # TODO: calculate proportional area of natural within 2km
    task_graph.add_task(
        func=reproduce_env.register_data,
        args=(
            'convolution_kernel_raster',
            'convolution_kernel.tif',
            lambda kernel_filepath: create_radial_convolution_mask(
                0.00277778, 2000., kernel_filepath)))

    task_graph.join()
    crop_nutrient_df = pandas.read_csv(reproduce_env['crop_nutrient_table'])
    globio_df = pandas.read_csv(reproduce_env['globio_class_table'])
    print crop_nutrient_df
    print globio_df
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


if __name__ == '__main__':
    main()
