"""
Pollination analysis for IPBES.

    From "IPBES Methods: Pollination Contribution to Human Nutrition."
    https://www.dropbox.com/s/gc4b1miw2zypuke/IPBES%20Methods_Pollination_RS.docx?dl=0
"""
import subprocess
import sys
import zipfile
import time
import os
import re
import logging
import itertools
import multiprocessing

from osgeo import gdal
import google.cloud.client
import google.cloud.storage
from osgeo import osr
import pandas
import numpy
import scipy.ndimage.morphology
import reproduce
import taskgraph
import pygeoprocessing

# set a 1GB limit for the cache
gdal.SetCacheMax(2**30)

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger('ipbes_pollination')

# the following are the globio landcover codes. A tuple (x, y) indicates the
# inclusive range from x to y. Pollinator habitat was defined as any natural
# land covers, as defined (GLOBIO land-cover classes 6, secondary vegetation,
# and  50-180, various types of primary vegetation). To test sensitivity to
# this definition we included "semi-natural" habitats (GLOBIO land-cover
# classes 3, 4, and 5; pasture, rangeland and forestry, respectively) in
# addition to "natural", and repeated all analyses with semi-natural  plus
# natural habitats, but this did not substantially alter the results  so we do
# not include it in our final analysis or code base.

GLOBIO_AG_CODES = [2, (230, 232)]
GLOBIO_NATURAL_CODES = [6, (50, 180)]

WORKING_DIR = 'workspace'
OUTPUT_DIR = os.path.join(WORKING_DIR, 'outputs')
ECOSHARD_DIR = os.path.join(WORKING_DIR, 'ecoshard_dir')
CHURN_DIR = os.path.join(WORKING_DIR, 'churn')
GOOGLE_BUCKET_KEY_PATH = "ecoshard-202992-key.json"
NODATA = -9999
N_WORKERS = max(1, multiprocessing.cpu_count())
DELAYED_START = N_WORKERS >= 0

GOOGLE_BUCKET_ID = 'ipbes-pollination-result'
MERCURIAL_HASH_ID = subprocess.check_output(
    ['hg', 'id', '--id']).strip().decode('utf-8')
if MERCURIAL_HASH_ID.endswith('+'):
    raise RuntimeError(
        f'Repository must be clean before running: {MERCURIAL_HASH_ID}')
BLOB_ROOT = f'''ipbes_pollination_result_{MERCURIAL_HASH_ID}'''


def main():
    """Entry point."""
    task_graph = taskgraph.TaskGraph(
        CHURN_DIR, N_WORKERS, delayed_start=DELAYED_START,
        reporting_interval=5.0)

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
        ECOSHARD_DIR, os.path.basename(crop_nutrient_url))

    crop_nutrient_table_fetch_task = task_graph.add_task(
        func=google_bucket_fetch_and_validate,
        args=(
            crop_nutrient_url, GOOGLE_BUCKET_KEY_PATH,
            crop_nutrient_table_path),
        target_path_list=[crop_nutrient_table_path],
        task_name=f'fetch {os.path.basename(crop_nutrient_table_path)}')

    landcover_data = {
        'GLOBIO4_LU_10sec_2050_SSP5_RCP85': (
            'https://storage.cloud.google.com/ecoshard-root/globio_landcover'
            '/GLOBIO4_LU_10sec_2050_SSP5_RCP85_'
            'md5_1b3cc1ce6d0ff14d66da676ef194f130.tif', 'ssp5'),
        'GLOBIO4_LU_10sec_2050_SSP1_RCP26': (
            'https://storage.cloud.google.com/ecoshard-root/globio_landcover'
            '/GLOBIO4_LU_10sec_2050_SSP1_RCP26_md5_'
            '803166420f51e5ef7dcaa970faa98173.tif', 'ssp1'),
        'GLOBIO4_LU_10sec_2050_SSP3_RCP70': (
            'https://storage.cloud.google.com/ecoshard-root/globio_landcover'
            '/GLOBIO4_LU_10sec_2050_SSP3_RCP70_md5_'
            'e77077a3220a36f7f0441bbd0f7f14ab.tif', 'ssp3'),
        'Globio4_landuse_10sec_1850': (
            'https://storage.cloud.google.com/ecoshard-root/globio_landcover'
            '/Globio4_landuse_10sec_1850_md5_'
            '0b7fcb4b180d46b4fc2245beee76d6b9.tif', '1850'),
        'Globio4_landuse_10sec_2015': (
            'https://storage.cloud.google.com/ecoshard-root/globio_landcover'
            '/Globio4_landuse_10sec_2015_md5_'
            '939a57c2437cd09bd5a9eb472b9bd781.tif', 'cur'),
        'Globio4_landuse_10sec_1980': (
            'https://storage.cloud.google.com/ecoshard-root/globio_landcover'
            '/Globio4_landuse_10sec_1980_md5_'
            'f6384eac7579318524439df9530ca1f4.tif', '1980'),
        'Globio4_landuse_10sec_1945': (
            'https://storage.cloud.google.com/ecoshard-root/globio_landcover'
            '/Globio4_landuse_10sec_1945_md5_'
            '52c7b4c38c26defefa61132fd25c5584.tif', '1945'),
        'Globio4_landuse_10sec_1910': (
            'https://storage.cloud.google.com/ecoshard-root/globio_landcover'
            '/Globio4_landuse_10sec_1910_md5_'
            'e7da8fa29db305ff63c99fed7ca8d5e2.tif', '1910'),
        'Globio4_landuse_10sec_1900': (
            'https://storage.cloud.google.com/ecoshard-root/globio_landcover'
            '/Globio4_landuse_10sec_1900_md5_'
            'f5db818a5b16799bf2cb627e574120a4.tif', '1900'),
        }

    # 1.2.3.  Crop production

    # Spatially-explicit global crop yields (tons/ha) at 5 arc min (~10 km)
    # were taken from Monfreda et al. (2008) for 115 crops (permanent link to
    # crop yield folder). These yields were multiplied by crop pollination
    # dependency to calculate the pollination-dependent crop yield for each 5
    # min grid cell. Note the monfredia maps are in units of per-hectare
    # yields

    yield_zip_url = (
        'https://storage.cloud.google.com/ecoshard-root/ipbes/'
        'monfreda_2008_observed_yield_md5_54c6b8e564973739ba75c8e54ac6f051.'
        'zip')

    yield_zip_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(yield_zip_url))

    yield_zip_fetch_task = task_graph.add_task(
        func=google_bucket_fetch_and_validate,
        args=(
            yield_zip_url, GOOGLE_BUCKET_KEY_PATH,
            yield_zip_path),
        target_path_list=[yield_zip_path],
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
        task_name=f'unzip monfreda_2008_observed_yield')

    # fetch a landcover map to use as a base for the dimensions of the
    # production raster, later on we'll fetch them all. TODO MOVE UP FETCH TO REPLACE THIS
    landcover_key, (landcover_url, _) = next(iter(landcover_data.items()))
    landcover_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(landcover_url))

    landcover_fetch_task = task_graph.add_task(
        func=google_bucket_fetch_and_validate,
        args=(landcover_url, GOOGLE_BUCKET_KEY_PATH, landcover_path),
        target_path_list=[landcover_path],
        task_name=f'fetch {landcover_key}')

    yield_raster_dir = os.path.join(
        os.path.dirname(yield_zip_path), 'monfreda_2008_observed_yield')

    prod_total_nut_10s_task_path_map = {}
    poll_dep_prod_nut_10s_task_path_map = {}
    for nutrient_id, nutrient_name in [
            ('en', 'Energy'), ('va', 'VitA'), ('fo', 'Folate')]:
        # total annual production of nutrient
        yield_total_nut_10km_path = os.path.join(
            CHURN_DIR, 'yield_nutrient_rasters',
            f'yield_total_{nutrient_id}_10km.tif')
        yield_total_nut_10s_path = os.path.join(
            CHURN_DIR, 'yield_nutrient_rasters',
            f'yield_total_{nutrient_id}_10s.tif')
        prod_total_nut_10s_path = os.path.join(
            CHURN_DIR, 'prod_nutrient_rasters',
            f'prod_total_{nutrient_id}_10s.tif')

        prod_total_task = task_graph.add_task(
            func=create_prod_nutrient_raster,
            args=(
                 crop_nutrient_table_path, nutrient_name,
                 yield_raster_dir, False, landcover_path,
                 yield_total_nut_10km_path, yield_total_nut_10s_path,
                 prod_total_nut_10s_path),
            target_path_list=[
                yield_total_nut_10km_path, yield_total_nut_10s_path,
                prod_total_nut_10s_path],
            dependent_task_list=[
                landcover_fetch_task, unzip_yield_task,
                crop_nutrient_table_fetch_task],
            task_name=f"""create prod raster {
                os.path.basename(prod_total_nut_10s_path)}""")
        for upload_path in [
                yield_total_nut_10km_path, yield_total_nut_10s_path,
                prod_total_nut_10s_path]:
            upload_blob(task_graph, upload_path, prod_total_task)
        prod_total_nut_10s_task_path_map[nutrient_id] = (
            prod_total_task, prod_total_nut_10s_path)

        # pollination-dependent annual production of nutrient
        poll_dep_yield_nut_10km_path = os.path.join(
            CHURN_DIR, 'yield_poll_dep_rasters',
            f'yield_poll_dep_{nutrient_id}_10km.tif')
        poll_dep_yield_nut_10s_path = os.path.join(
            CHURN_DIR, 'yield_poll_dep_rasters',
            f'yield_poll_dep_{nutrient_id}_10s.tif')
        poll_dep_prod_nut_10s_path = os.path.join(
            CHURN_DIR, 'prod_poll_dep_rasters',
            f'prod_poll_dep_total_{nutrient_id}_10s.tif')
        pol_dep_prod_task = task_graph.add_task(
            func=create_prod_nutrient_raster,
            args=(
                crop_nutrient_table_path, nutrient_name,
                yield_raster_dir, True, landcover_path,
                poll_dep_yield_nut_10km_path, poll_dep_yield_nut_10s_path,
                poll_dep_prod_nut_10s_path),
            target_path_list=[
                poll_dep_yield_nut_10km_path, poll_dep_yield_nut_10s_path,
                poll_dep_prod_nut_10s_path],
            dependent_task_list=[landcover_fetch_task, unzip_yield_task],
            task_name=f"""create poll dep production raster {
                os.path.basename(poll_dep_prod_nut_10s_path)}""")
        for upload_path in [
                poll_dep_yield_nut_10km_path, poll_dep_yield_nut_10s_path,
                poll_dep_prod_nut_10s_path]:
            upload_blob(task_graph, upload_path, pol_dep_prod_task)
        poll_dep_prod_nut_10s_task_path_map[nutrient_id] = (
            pol_dep_prod_task, poll_dep_prod_nut_10s_path)

    # The proportional area of natural within 2 km was calculated for every
    #  pixel of agricultural land (GLOBIO land-cover classes 2, 230, 231, and
    #  232) at 10 arc seconds (~300 m) resolution. This 2 km scale represents
    #  the distance most commonly found to be predictive of pollination
    #  services (Kennedy et al. 2013).
    kernel_raster_path = os.path.join(CHURN_DIR, 'radial_kernel.tif')
    kernel_task = task_graph.add_task(
        func=create_radial_convolution_mask,
        args=(0.00277778, 2000., kernel_raster_path),
        target_path_list=[kernel_raster_path],
        task_name='make convolution kernel')

    # mask landcover into agriculture and pollinator habitat
    for landcover_key, (landcover_url, landcover_short_suffix) in (
            landcover_data.items()):
        landcover_path = os.path.join(
            ECOSHARD_DIR, os.path.basename(landcover_url))
        landcover_fetch_task = task_graph.add_task(
            func=google_bucket_fetch_and_validate,
            args=(landcover_url, GOOGLE_BUCKET_KEY_PATH, landcover_path),
            target_path_list=[landcover_path],
            task_name=f'fetch {landcover_key}')
        schedule_build_overviews(
            task_graph, landcover_path, landcover_fetch_task)

        # This loop is so we don't duplicate code for 'ag' and 'hab' with the
        # only difference being
        for mask_prefix, globio_codes in [
                ('ag', GLOBIO_AG_CODES), ('hab', GLOBIO_NATURAL_CODES)]:
            mask_key = f'{landcover_key}_{mask_prefix}_mask'
            mask_target_path = os.path.join(
                CHURN_DIR, f'{mask_prefix}_mask', f'{mask_key}.tif')

            mask_task = task_graph.add_task(
                func=mask_raster,
                args=(landcover_path, globio_codes, mask_target_path),
                target_path_list=[mask_target_path],
                dependent_task_list=[landcover_fetch_task],
                task_name=f'mask {mask_key}',)
            upload_blob(task_graph, mask_target_path, mask_task)
            schedule_build_overviews(task_graph, mask_target_path, mask_task)

            if mask_prefix == 'hab':
                hab_task_path_tuple = (mask_task, mask_target_path)
            elif mask_prefix == 'ag':
                ag_task_path_tuple = (mask_task, mask_target_path)

        proportional_hab_area_2km_path = os.path.join(
            CHURN_DIR, 'proportional_area',
            f'{landcover_key}_hab_prop_area_2km.tif')
        prop_hab_area_2km_task = task_graph.add_task(
            func=pygeoprocessing.convolve_2d,
            args=[
                (hab_task_path_tuple[1], 1), (kernel_raster_path, 1),
                proportional_hab_area_2km_path],
            kwargs={
                'working_dir': CHURN_DIR,
                'gtiff_creation_options': (
                    'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
                    'PREDICTOR=3', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256',
                    'NUM_THREADS=2'),
                'n_threads': 2},
            dependent_task_list=[hab_task_path_tuple[0], kernel_task],
            target_path_list=[proportional_hab_area_2km_path],
            task_name=(
                'calculate proportional'
                f' {os.path.basename(proportional_hab_area_2km_path)}'))
        upload_blob(
            task_graph, proportional_hab_area_2km_path,
            prop_hab_area_2km_task)
        schedule_build_overviews(
            task_graph, proportional_hab_area_2km_path,
            prop_hab_area_2km_task)

        #  1.1.4.  Sufficiency threshold A threshold of 0.3 was set to
        #  evaluate whether there was sufficient pollinator habitat in the 2
        #  km around farmland to provide pollination services, based on
        #  Kremen et al.'s (2005)  estimate of the area requirements for
        #  achieving full pollination. This produced a map of wild
        #  pollination sufficiency where every agricultural pixel was
        #  designated in a binary fashion: 0 if proportional area of habitat
        #  was less than 0.3; 1 if greater than 0.3. Maps of pollination
        #  sufficiency can be found at (permanent link to output), outputs
        #  "poll_suff_..." below.

        threshold_val = 0.3
        pollinator_suff_hab_path = os.path.join(
            CHURN_DIR, 'poll_suff_hab_ag_coverage_rasters',
            f'poll_suff_ag_coverage_mask_10s_{landcover_short_suffix}.tif')
        poll_suff_task = task_graph.add_task(
            func=threshold_select_raster,
            args=(
                proportional_hab_area_2km_path,
                ag_task_path_tuple[1], threshold_val,
                pollinator_suff_hab_path),
            target_path_list=[pollinator_suff_hab_path],
            dependent_task_list=[
                prop_hab_area_2km_task, ag_task_path_tuple[0]],
            task_name=f"""poll_suff_ag_coverage_mask {
                os.path.basename(pollinator_suff_hab_path)}""")
        upload_blob(task_graph, pollinator_suff_hab_path, poll_suff_task)
        schedule_build_overviews(
            task_graph, pollinator_suff_hab_path, poll_suff_task)

        # tot_prod_en|va|fo_10s|1d_cur|ssp1|ssp3|ssp5
        # total annual production of energy (KJ/yr), vitamin A (IU/yr),
        # and folate (mg/yr)
        for nutrient_id in ['en', 'va', 'fo']:
            tot_prod_task, tot_prod_1d_path = (
                prod_total_nut_10s_task_path_map[nutrient_id])

            prod_total_path = os.path.join(
                OUTPUT_DIR, f'''prod_total_{
                    nutrient_id}_10s_{landcover_short_suffix}.tif''')
            tot_prod_nut_scenario_task = task_graph.add_task(
                func=mult_rasters,
                args=(
                    ag_task_path_tuple[1], tot_prod_1d_path,
                    prod_total_path),
                target_path_list=[prod_total_path],
                dependent_task_list=[tot_prod_task, ag_task_path_tuple[0]],
                task_name=(
                    f'tot_prod_{nutrient_id}_10s_{landcover_short_suffix}'))
            upload_blob(
                task_graph, prod_total_path, tot_prod_nut_scenario_task)
            schedule_build_overviews(
                task_graph, prod_total_path, tot_prod_nut_scenario_task)

            poll_dep_prod_task, poll_dep_prod_path = (
                poll_dep_prod_nut_10s_task_path_map[nutrient_id])

            prod_poll_dep_total_nut_scenario_path = os.path.join(
                WORKING_DIR, OUTPUT_DIR,
                f'prod_poll_dep_total_{nutrient_id}_10s_'
                f'{landcover_short_suffix}.tif')
            poll_dep_prod_nut_scenario_task = task_graph.add_task(
                func=mult_rasters,
                args=(
                    ag_task_path_tuple[1], poll_dep_prod_path,
                    prod_poll_dep_total_nut_scenario_path),
                target_path_list=[prod_poll_dep_total_nut_scenario_path],
                dependent_task_list=[
                    poll_dep_prod_task, ag_task_path_tuple[0]],
                task_name=(
                    f'poll_dep_prod_{nutrient_id}_'
                    f'10s_{landcover_short_suffix}'))
            upload_blob(
                task_graph, prod_poll_dep_total_nut_scenario_path,
                poll_dep_prod_nut_scenario_task)
            schedule_build_overviews(
                task_graph, prod_poll_dep_total_nut_scenario_path,
                poll_dep_prod_nut_scenario_task)

            # prod_poll_dep_realized_en|va|fo_10s|1d_cur|ssp1|ssp3|ssp5:
            # pollination-dependent annual production of energy (KJ/yr),
            # vitamin A (IU/yr), and folate (mg/yr) that can be met by wild
            # pollinators due to the proximity of sufficient habitat.
            prod_poll_dep_realized_nut_scenario_path = os.path.join(
                WORKING_DIR, OUTPUT_DIR,
                f'prod_poll_dep_realized_{nutrient_id}_10s_'
                f'{landcover_short_suffix}.tif')
            prod_poll_dep_realized_nut_scenario_task = task_graph.add_task(
                func=mult_rasters,
                args=(
                    prod_poll_dep_total_nut_scenario_path,
                    pollinator_suff_hab_path,
                    prod_poll_dep_realized_nut_scenario_path),
                target_path_list=[prod_poll_dep_realized_nut_scenario_path],
                dependent_task_list=[
                    poll_suff_task, poll_dep_prod_nut_scenario_task],
                task_name=(
                    f'prod_poll_dep_realized_{nutrient_id}_'
                    f'10s_{landcover_short_suffix}'))
            upload_blob(
                task_graph, prod_poll_dep_realized_nut_scenario_path,
                prod_poll_dep_realized_nut_scenario_task)
            schedule_build_overviews(
                task_graph, prod_poll_dep_realized_nut_scenario_path,
                prod_poll_dep_realized_nut_scenario_task)

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

    spatial_population_scenarios_url = (
        'https://storage.cloud.google.com/ecoshard-root/ipbes/'
        'Spatial_population_scenarios_GeoTIFF_'
        'md5_1c4b6d87cb9a167585e1fc49914248fd.zip')

    spatial_population_scenarios_path = os.path.join(
        ECOSHARD_DIR, 'spatial_population_scenarios',
        os.path.basename(spatial_population_scenarios_url))

    spatial_population_scenarios_fetch_task = task_graph.add_task(
        func=google_bucket_fetch_and_validate,
        args=(
            spatial_population_scenarios_url, GOOGLE_BUCKET_KEY_PATH,
            spatial_population_scenarios_path),
        target_path_list=[spatial_population_scenarios_path],
        task_name=f"""fetch {os.path.basename(
            spatial_population_scenarios_path)}""",
        priority=100)

    spatial_scenarios_pop_zip_touch_file_path = os.path.join(
        os.path.dirname(spatial_population_scenarios_path),
        'Spatial_population_scenarios_GeoTIFF.txt')
    unzip_spatial_population_scenarios_task = task_graph.add_task(
        func=unzip_file,
        args=(
            spatial_population_scenarios_path,
            os.path.dirname(spatial_population_scenarios_path),
            spatial_scenarios_pop_zip_touch_file_path),
        target_path_list=[spatial_scenarios_pop_zip_touch_file_path],
        dependent_task_list=[spatial_population_scenarios_fetch_task],
        task_name=f'unzip Spatial_population_scenarios_GeoTIFF',
        priority=100)

    gpw_urls = {
        'gpw_v4_e_a000_014ft_2010_dens': (
            'https://storage.cloud.google.com/ecoshard-root/ipbes/'
            'gpw_v4_e_a000_014ft_2010_dens_30_sec_md5_'
            '8b2871967b71534d56d2df83e413bf33.tif'),
        'gpw_v4_e_a000_014mt_2010_dens': (
            'https://storage.cloud.google.com/ecoshard-root/ipbes/'
            'gpw_v4_e_a000_014mt_2010_dens_30_sec_md5_'
            '8ddf234eabc0025efd5678776e2ae792.tif'),
        'gpw_v4_e_a065plusft_2010_dens': (
            'https://storage.cloud.google.com/ecoshard-root/ipbes/'
            'gpw_v4_e_a065plusft_2010_dens_30_sec_md5_'
            'c0cfbc8685dbf16cbe520aac963cd6b6.tif'),
        'gpw_v4_e_a065plusmt_2010_dens': (
            'https://storage.cloud.google.com/ecoshard-root/ipbes/'
            'gpw_v4_e_a065plusmt_2010_dens_30_sec_md5_'
            '1d36e79aa083ee25de7295a4a7d10a4b.tif'),
        'gpw_v4_e_atotpopft_2010_dens': (
            'https://storage.cloud.google.com/ecoshard-root/ipbes/'
            'gpw_v4_e_atotpopft_2010_dens_30_sec_md5_'
            '5bbb72a050c76264e1b6a3c7530fedce.tif'),
        'gpw_v4_e_atotpopmt_2010_dens': (
            'https://storage.cloud.google.com/ecoshard-root/ipbes/'
            'gpw_v4_e_atotpopmt_2010_dens_30_sec_md5_'
            '31637ca784b8b917222661d4a915ead6.tif'),
    }

    gpw_task_path_id_map = {}
    for gpw_id, gpw_url in gpw_urls.items():
        gpw_dens_path = os.path.join(
            ECOSHARD_DIR, 'gpw_pop_densities', os.path.basename(gpw_url))
        gpw_fetch_task = task_graph.add_task(
            func=google_bucket_fetch_and_validate,
            args=(gpw_url, GOOGLE_BUCKET_KEY_PATH, gpw_dens_path),
            target_path_list=[gpw_dens_path],
            task_name=f"""fetch {os.path.basename(gpw_dens_path)}""",
            priority=100)
        schedule_build_overviews(task_graph, gpw_dens_path, gpw_fetch_task)

        gpw_count_path = os.path.join(
            CHURN_DIR, 'gpw_count', f"""{gpw_id[:-4]}count.tif""")
        gpw_count_task = task_graph.add_task(
            func=calc_pop_count,
            args=(gpw_dens_path, gpw_count_path),
            target_path_list=[gpw_count_path],
            dependent_task_list=[gpw_fetch_task],
            task_name=f"""pop count {os.path.basename(gpw_count_path)}""")
        gpw_task_path_id_map[f'{gpw_id[:-4]}count'] = (
            gpw_count_task, gpw_count_path)
        upload_blob(task_graph, gpw_count_path, gpw_count_task)
        schedule_build_overviews(task_graph, gpw_count_path, gpw_count_task)

    # calculate 15-65 population gpw count by subtracting total from
    # 0-14 and 65plus
    for gender_id in ['f', 'm']:
        gpw_v4_e_a15_65t_2010_count_path = os.path.join(
            WORKING_DIR, 'gpw_count',
            f'gpw_v4_e_a015_065{gender_id}_2010_count.tif')
        gpw_15_65f_count_task = task_graph.add_task(
            func=subtract_rasters,
            args=(
                gpw_task_path_id_map[
                    f'gpw_v4_e_atotpop{gender_id}t_2010_count'][1],
                gpw_task_path_id_map[
                    f'gpw_v4_e_a000_014{gender_id}t_2010_count'][1],
                gpw_task_path_id_map[
                    f'gpw_v4_e_a065plus{gender_id}t_2010_count'][1],
                gpw_v4_e_a15_65t_2010_count_path),
            target_path_list=[gpw_v4_e_a15_65t_2010_count_path],
            dependent_task_list=[
                gpw_task_path_id_map[
                    f'gpw_v4_e_atotpop{gender_id}t_2010_count'][0],
                gpw_task_path_id_map[
                    f'gpw_v4_e_a000_014{gender_id}t_2010_count'][0],
                gpw_task_path_id_map[
                    f'gpw_v4_e_a065plus{gender_id}t_2010_count'][0]],
            task_name=f'calc gpw 15-65 {gender_id}')
        gpw_task_path_id_map[f'gpw_v4_e_a015_065{gender_id}t_2010_count'] = (
            gpw_15_65f_count_task, gpw_v4_e_a15_65t_2010_count_path)

        schedule_build_overviews(
            task_graph, gpw_v4_e_a15_65t_2010_count_path,
            gpw_15_65f_count_task)

    # we need to warp SSP rasters to match the GPW rasters
    # we need to calcualte 15-65 pop by subtracting 0-14 and 65 plus from tot
    # then we can calculate SSP future for 0-14, 15-65, and 65 plus
    # then we can calculate nutritional needs for cur, and ssp scenarios
    ssp_task_pop_map = {}
    for ssp_id in (1, 3, 5):
        spatial_pop_dir = os.path.join(
            os.path.dirname(spatial_population_scenarios_path),
            'Spatial_population_scenarios_GeoTIFF',
            f'SSP{ssp_id}_GeoTIFF', 'total', 'GeoTIFF')
        cur_ssp_path = os.path.join(spatial_pop_dir, f'ssp{ssp_id}_2010.tif')
        fut_ssp_path = os.path.join(spatial_pop_dir, f'ssp{ssp_id}_2050.tif')

        cur_ssp_warp_path = f'{os.path.splitext(cur_ssp_path)[0]}_gpwwarp.tif'
        fut_ssp_warp_path = f'{os.path.splitext(fut_ssp_path)[0]}_gpwwarp.tif'
        warp_task_list = []
        # we need a canonical base gpw raster to warp to, just grab the
        # "first" one off the dict.
        (gpw_base_tot_count_task, gpw_base_tot_count_path) = (
            gpw_task_path_id_map.values().__iter__().__next__())
        for base_path, warp_path in [
                (cur_ssp_path, cur_ssp_warp_path),
                (fut_ssp_path, fut_ssp_warp_path)]:
            warp_task = task_graph.add_task(
                func=warp_to_raster,
                args=(
                    base_path, gpw_base_tot_count_path, 'bilinear',
                    warp_path),
                dependent_task_list=[
                    gpw_base_tot_count_task,
                    unzip_spatial_population_scenarios_task],
                target_path_list=[warp_path],
                task_name=f'warp to raster {os.path.basename(warp_path)}')
            warp_task_list.append(warp_task)
            schedule_build_overviews(task_graph, warp_path, warp_task)

        for gpw_id in [
                'gpw_v4_e_a000_014ft_2010_count',
                'gpw_v4_e_a000_014mt_2010_count',
                'gpw_v4_e_a065plusft_2010_count',
                'gpw_v4_e_a065plusmt_2010_count',
                'gpw_v4_e_a015_065ft_2010_count',
                'gpw_v4_e_a015_065mt_2010_count',
                ]:

            gpw_task, gpw_tot_count_path = (
                gpw_task_path_id_map[gpw_id])
            ssp_pop_path = os.path.join(
                WORKING_DIR, 'gpw_ssp_rasters', f'ssp{ssp_id}_{gpw_id}.tif')

            ssp_pop_task = task_graph.add_task(
                func=calculate_future_pop,
                args=(
                    cur_ssp_warp_path, fut_ssp_warp_path, gpw_tot_count_path,
                    ssp_pop_path),
                dependent_task_list=warp_task_list + [gpw_task],
                target_path_list=[ssp_pop_path],
                task_name=f'ssp pop {os.path.basename(ssp_pop_path)}')
            ssp_task_pop_map[(ssp_id, gpw_id)] = (ssp_pop_task, ssp_pop_path)
            schedule_build_overviews(task_graph, ssp_pop_path, ssp_pop_task)

    # 2)
    # calculate the total nutritional needs per pixel for cur ssp1..5 scenario
    # tot_req_en|va|fo_10s|1d_cur|ssp1|ssp3|ssp5

    # this table comes from section "1.4.3 Dietary requirements"
    # units are 'va': Vitamin A (mcg RE)
    # 'fo': folate (mcg DFE)
    # 'en': Energy (kcal)
    nutritional_needs_map = {
        'gpw_v4_e_a000_014ft_2010_count': {'va': 450, 'fo': 250, 'en': 1531},
        'gpw_v4_e_a000_014mt_2010_count': {'va': 483, 'fo': 250, 'en': 1648},
        'gpw_v4_e_a015_065ft_2010_count': {'va': 516, 'fo': 408, 'en': 2153},
        'gpw_v4_e_a015_065mt_2010_count': {'va': 600, 'fo': 400, 'en': 2675},
        'gpw_v4_e_a065plusft_2010_count': {'va': 500, 'fo': 400, 'en': 1876},
        'gpw_v4_e_a065plusmt_2010_count': {'va': 600, 'fo': 400, 'en': 2318},
    }
    for nut_id in ('en', 'va', 'fo'):
        # calculate 'cur' needs
        pop_task_path_list, nut_need_list = zip(*[
            (gpw_task_path_id_map[gpw_id],
             nutritional_needs_map[gpw_id][nut_id])
            for gpw_id in nutritional_needs_map])
        pop_task_list, pop_path_list = zip(*pop_task_path_list)

        tot_nut_requirements_path = os.path.join(
            OUTPUT_DIR,
            f'nut_req_local_{nut_id}_10s_cur.tif')
        total_requirements_task = task_graph.add_task(
            func=calculate_total_requirements,
            args=(
                pop_path_list, nut_need_list, tot_nut_requirements_path),
            target_path_list=[tot_nut_requirements_path],
            dependent_task_list=pop_task_list,
            task_name=f"""tot nut requirements {
                os.path.basename(tot_nut_requirements_path)}""",
            priority=100,)
        upload_blob(
            task_graph, tot_nut_requirements_path, total_requirements_task)
        schedule_build_overviews(
                task_graph, tot_nut_requirements_path, total_requirements_task)

        # calculate ssp needs
        for ssp_id in (1, 3, 5):
            pop_task_path_list, nut_need_list = zip(*[
                (ssp_task_pop_map[(ssp_id, gpw_id)],
                 nutritional_needs_map[gpw_id][nut_id])
                for gpw_id in nutritional_needs_map])
            pop_task_list, pop_path_list = zip(*pop_task_path_list)

            nut_req_path = os.path.join(
                WORKING_DIR, OUTPUT_DIR,
                f'nut_req_{nut_id}_10s_ssp{ssp_id}.tif')
            nut_req_task = task_graph.add_task(
                func=calculate_total_requirements,
                args=(
                    pop_path_list, nut_need_list, nut_req_path),
                target_path_list=[nut_req_path],
                dependent_task_list=pop_task_list,
                task_name=f"""tot nut requirements {
                    os.path.basename(nut_req_path)}""",
                priority=100,)
            upload_blob(task_graph, nut_req_path, nut_req_task)
            schedule_build_overviews(
                task_graph, nut_req_path, nut_req_task)

    task_graph.close()
    task_graph.join()
    # END MAIN


def schedule_build_overviews(task_graph, base_raster_path, base_raster_task):
    """Build overviews of raster path using taskgraph to schedule.

    Creates an external overview of `base_raster_path` raster in the same
    directory and filename as the original with a '.ovr' extension.

    Parameters:
        task_graph (TaskGraph): TaskGraph object to schedule through.
        base_raster_path (str): path to raster to build overviews on.
        base_raster_task (Task): task to wait for `base_raster_path` to be
            created on.

    Returns:
        None.

    """
    target_overview_path = f'{base_raster_path}.ovr'
    overview_task = task_graph.add_task(
        func=build_overviews,
        priority=-100,
        args=(base_raster_path, 'nearest'),
        target_path_list=[target_overview_path],
        dependent_task_list=[base_raster_task],
        task_name=f'overviews {os.path.basename(base_raster_path)}')
    upload_blob(task_graph, target_overview_path, overview_task)


def upload_blob(task_graph, base_path, dependent_task):
    """Upload file to blob root.

    Parameters:
        task_graph (taskgraph.TaskGraph): TaskGraph object to schedule.
        base_path (str): path to local file to upload. Blob id will end with
            the base filename of this path.
        dependent_task (taskgraph.Task): task that creates ``base_path``.

    Returns:
        None.

    """
    # get the relative blob path from the workspace and use Google blob
    # notation for directories
    blob_path = (
        os.path.join(BLOB_ROOT, os.path.relpath(
            base_path, WORKING_DIR)).replace(os.sep, '/'))
    _ = task_graph.add_task(
        func=google_bucket_upload,
        args=(
            base_path, GOOGLE_BUCKET_ID, blob_path,
            GOOGLE_BUCKET_KEY_PATH),
        dependent_task_list=[dependent_task],
        task_name=f'google bucket upload {blob_path}')


def calculate_total_requirements(
        pop_path_list, nut_need_list, target_path):
    """Calculate total nutrient requirements.

    Create a new raster by summing all rasters in `pop_path_list` multiplied
    by their corresponding scalar in `nut_need_list`.

    Parameters:
        pop_path_list (list of str): list of paths to population counts.
        nut_need_list (list): list of scalars that correspond in order to
            the per-count nutrient needs of `pop_path_list`.
        target_path (str): path to target file.

    Return:
        None.

    """
    nodata = -1
    pop_nodata = pygeoprocessing.get_raster_info(
        pop_path_list[0])['nodata'][0]

    def mult_and_sum(*arg_list):
        """Arg list is an (array0, scalar0, array1, scalar1,...) list.

        Returns:
            array0*scalar0 + array1*scalar1 + .... but ignore nodata.

        """
        result = numpy.empty(arg_list[0].shape, dtype=numpy.float32)
        result[:] = nodata
        array_stack = numpy.array(arg_list[0::2])
        scalar_list = numpy.array(arg_list[1::2])
        # make a valid mask as big as a single array
        valid_mask = numpy.logical_and.reduce(
            array_stack == pop_nodata, axis=0)

        # mask out all invalid elements but reshape so there's still the same
        # number of arrays
        valid_array_elements = (
            array_stack[numpy.broadcast_to(valid_mask, array_stack.shape)])
        array_stack = None

        # sometimes this array is empty, check first before reshaping
        if valid_array_elements:
            valid_array_elements = valid_array_elements.reshape(
                -1, numpy.count_nonzero(valid_mask))
            # multiply each element of the scalar with each row of the valid
            # array stack, then sum along the 0 axis to get the result
            result[valid_mask] = numpy.sum(
                (valid_array_elements.T * scalar_list).T, axis=0)
        scalar_list = None
        valid_mask = None
        valid_array_elements = None
        return result

    pygeoprocessing.raster_calculator(list(itertools.chain(*[
        ((path, 1), (scalar, 'raw')) for path, scalar in zip(
            pop_path_list, nut_need_list)])), mult_and_sum, target_path,
        gdal.GDT_Float32, nodata)


def subtract_rasters(
        raster_path_a, raster_path_b, raster_path_c, target_path):
    """Calculate target = a-b-c and ignore nodata."""
    a_nodata = pygeoprocessing.get_raster_info(raster_path_a)['nodata'][0]
    b_nodata = pygeoprocessing.get_raster_info(raster_path_b)['nodata'][0]
    c_nodata = pygeoprocessing.get_raster_info(raster_path_c)['nodata'][0]
    target_nodata = -9999

    def sub_op(a_array, b_array, c_array):
        """Sub a-b-c as arrays."""
        result = numpy.empty(a_array.shape, dtype=numpy.float32)
        result[:] = target_nodata
        valid_mask = (
            (a_array != a_nodata) &
            (b_array != b_nodata) &
            (c_array != c_nodata))
        result[valid_mask] = (
            a_array[valid_mask] - b_array[valid_mask] - c_array[valid_mask])
        return result

    pygeoprocessing.raster_calculator(
        [(raster_path_a, 1), (raster_path_b, 1), (raster_path_c, 1)],
        sub_op, target_path, gdal.GDT_Float32, target_nodata)


def warp_to_raster(
        base_raster_path, canonical_raster_path, resample_method,
        target_raster_path):
    """Warp base raster to canonical example.

    Parameters:
        base_raster_path (str): path to base raster
        canonical_raster_path (str),
        resample_method (str): one of nearest, bilinear, cubic, cubic_spline,
            lanczos, average, mode, max, min, med, q1, q3.
        target_raster_path (str): path to target warped raster.

    Returns:
        None.

    """
    canonical_raster_info = pygeoprocessing.get_raster_info(
        canonical_raster_path)
    pygeoprocessing.warp_raster(
        base_raster_path, canonical_raster_info['pixel_size'],
        target_raster_path,
        resample_method, target_bb=canonical_raster_info['bounding_box'],
        n_threads=2)


def calculate_future_pop(
        cur_ssp_path, fut_ssp_path, gpw_tot_count_path, target_ssp_pop_path):
    """Calculate future population raster.

    Multiply `gpw_tot_count` by the fut_ssp/cur_ssp ratio.
    """
    target_nodata = -1
    ssp_cur_nodata = pygeoprocessing.get_raster_info(
        cur_ssp_path)['nodata'][0]
    ssp_fut_nodata = pygeoprocessing.get_raster_info(
        fut_ssp_path)['nodata'][0]
    count_nodata = pygeoprocessing.get_raster_info(
        gpw_tot_count_path)['nodata'][0]

    def _future_pop_op(cur_array, fut_array, count_array):
        """Calculate future pop by dividing fut/cur*cur_count."""
        result = numpy.empty(cur_array.shape, dtype=numpy.float32)
        result[:] = count_nodata
        zero_mask = cur_array == 0
        valid_mask = (
            (cur_array != ssp_cur_nodata) &
            (fut_array != ssp_fut_nodata) &
            (count_array != count_nodata) & ~zero_mask)
        result[valid_mask] = (
            (fut_array[valid_mask] / cur_array[valid_mask]).astype(
                numpy.float32) * count_array[valid_mask])
        # assume if denominator is 0 we don't want to mask anything out, just
        # get it working
        result[zero_mask] = 1.0
        return result

    pygeoprocessing.raster_calculator(
        [(cur_ssp_path, 1), (fut_ssp_path, 1), (gpw_tot_count_path, 1)],
        _future_pop_op, target_ssp_pop_path, gdal.GDT_Float32, target_nodata)


def calc_pop_count(gpw_dens_path, gpw_count_path):
    """Calculate population count from density."""
    gpw_dens_info = pygeoprocessing.get_raster_info(gpw_dens_path)
    y_lat_array = numpy.linspace(
        gpw_dens_info['geotransform'][3],
        gpw_dens_info['geotransform'][3] +
        gpw_dens_info['geotransform'][5] *
        gpw_dens_info['raster_size'][1],
        gpw_dens_info['raster_size'][1])

    y_ha_array = area_of_pixel(
        abs(gpw_dens_info['geotransform'][1]),
        y_lat_array) / 10000.0
    y_ha_column = y_ha_array.reshape((y_ha_array.size, 1))

    nodata = gpw_dens_info['nodata'][0]
    pygeoprocessing.raster_calculator(
        [(gpw_dens_path, 1), y_ha_column, (nodata, 'raw')],
        density_to_value_op, gpw_count_path, gdal.GDT_Float32,
        nodata)


def create_radial_convolution_mask(
        pixel_size_degree, radius_meters, kernel_filepath):
    """Create a radial mask to sample pixels in convolution filter.

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
        json_key_path (string): path to Google Cloud private key generated by
        https://cloud.google.com/iam/docs/creating-managing-service-account-keys
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
    blob = google.cloud.storage.Blob(
        blob_id, bucket, chunk_size=2**24)
    LOGGER.info(f'downloading blob {target_path} from {url}')
    try:
        os.makedirs(os.path.dirname(target_path))
    except os.error:
        pass
    blob.download_to_filename(target_path)
    if not reproduce.valid_hash(target_path, 'embedded'):
        raise ValueError(f"{target_path}' does not match its expected hash")


def google_bucket_upload(
        local_file_path, bucket_id, blob_path, json_key_path):
    """Upload a local file to a given bucket.

    Uploads ``local_file_path`` to a blob registered on the ``json_key_path``
    account. The blob will be named as the filename of ``local_file_path``
    prefixed with the directory of ``__file__[mercurial_hash]/``

    Parameters:
        local_file_path (string): path to file to upload.
        bucket_id (string): name of Google Bucket.
        blob_id (string): name of the blob to upload to the bucket.
        json_key_path (string): path to Google Cloud private key generated by
            https://cloud.google.com/iam/docs/creating-managing-service-account-keys

    Returns:
        None.

    """
    client = google.cloud.storage.client.Client.from_service_account_json(
        json_key_path)
    bucket = client.get_bucket(bucket_id)
    # limit chunk size so we don't try to load the entire thing into memory
    blob = bucket.blob(blob_path, chunk_size=2**24)
    if not blob.exists():
        LOGGER.info(f'uploading blob {local_file_path} to {blob_path}')
        blob.upload_from_filename(local_file_path)
    else:
        LOGGER.info(f'not uploading {blob_path} exists')


def mask_codes_op(base_array, codes_array):
    """Return a bool raster if value in base_array is in codes_array."""
    return numpy.isin(base_array, codes_array)


def threshold_select_raster(
        base_raster_path, select_raster_path, threshold_val, target_path):
    """Select `select` if `base` >= `threshold_val`.

    Parameters:
        base_raster_path (string): path to single band raster that will be
            used to determine the threshold mask to select from
            `select_raster_path`.
        select_raster_path (string): path to single band raster to pass
            through to target if aligned `base` pixel is >= `threshold_val`
            0 otherwise, or nodata if base == nodata. Must be the same
            shape as `base_raster_path`.
        threshold_val (numeric): value to use as threshold cutoff
        target_path (string): path to desired output raster, raster is a
            byte type with same dimensions and projection as
            `base_raster_path`. A pixel in this raster will be `select` if
            the corresponding pixel in `base_raster_path` is >=
            `threshold_val`, 0 otherwise or nodata if `base` == nodata.

    Returns:
        None.

    """
    base_nodata = pygeoprocessing.get_raster_info(
        base_raster_path)['nodata'][0]
    target_nodata = 2

    def threshold_select_op(
            base_array, select_array, threshold_val, base_nodata,
            target_nodata):
        result = numpy.empty(select_array.shape, dtype=numpy.float32)
        result[:] = target_nodata
        valid_mask = base_array != base_nodata
        result[valid_mask] = numpy.where(
            base_array[valid_mask] >= threshold_val,
            select_array[valid_mask], 0)
        return result

    pygeoprocessing.raster_calculator(
        [(base_raster_path, 1), (select_raster_path, 1),
         (threshold_val, 'raw'), (base_nodata, 'raw'),
         (target_nodata, 'raw')], threshold_select_op,
        target_path, gdal.GDT_Byte, target_nodata, gtiff_creation_options=(
            'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
            'PREDICTOR=2', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256',
            'NUM_THREADS=2'))


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
    code_list = numpy.array([
        item for sublist in [
            range(x[0], x[1]+1) if isinstance(x, tuple) else [x]
            for x in codes] for item in sublist])
    LOGGER.debug(f'expanded code array {code_list}')

    pygeoprocessing.raster_calculator(
        [(base_path, 1), (code_list, 'raw')], mask_codes_op, target_path,
        gdal.GDT_Byte, 2, gtiff_creation_options=(
            'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
            'PREDICTOR=2', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256',
            'NUM_THREADS=2'))


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
    """Build as many overviews as possible.

    Parameters:
        local_path (string): path to GTiff raster to build overviews on. The
            overview will be built externally in a .ovr file in the same
            directory.
        resample_method (string): interpolation mode for overviews, must be
        one of 'nearest', 'average', 'gauss', 'cubic', 'cubicspline',
        'lanczos', 'average_mp', 'average_magphase', 'mode'.

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
        """Log updates using GDAL API for callbacks."""
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


def total_yield_op(
        yield_nodata, pollination_yield_factor_list, *crop_yield_array_list):
    """Calculate total yield."""
    result = numpy.empty(crop_yield_array_list[0].shape, dtype=numpy.float32)
    result[:] = 0.0
    all_valid = numpy.zeros(result.shape, dtype=numpy.bool)

    for crop_index, crop_array in enumerate(crop_yield_array_list):
        valid_mask = crop_array != yield_nodata
        all_valid |= valid_mask
        result[valid_mask] += (
            crop_array[valid_mask] *
            pollination_yield_factor_list[crop_index])
    result[~all_valid] = yield_nodata
    return result


def density_to_value_op(density_array, y_ha_array, density_nodata):
    """Calculate production."""
    result = numpy.empty(density_array.shape, dtype=numpy.float32)
    result[:] = density_nodata
    valid_mask = density_array != density_nodata
    result[valid_mask] = density_array[valid_mask] * y_ha_array[valid_mask]
    return result


def create_prod_nutrient_raster(
        crop_nutrient_df_path, nutrient_name, yield_raster_dir,
        consider_pollination, sample_target_raster_path,
        target_10km_yield_path, target_10s_yield_path,
        target_10s_production_path):
    """Create total production & yield for a nutrient for all crops.

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
        sample_target_raster_path (path): path to a file that has the raster
            pixel size and dimensions of the desired
            `target_10s_production_path`.
        sample_target_fetch_task (Task): must be complete before
            `sample_target_raster_path` is available.
        target_10km_yield_path (str): path to target raster that will
            contain total yield
        target_10s_yield_path (str): path to a resampled
            `target_10km_yield_path` at 10s resolution.
        target_10s_production_path (str): path to target raster that will
            contain a per-pixel amount of pollinator produced `nutrient_name`
            calculated as the sum(
                crop_yield_map * (100-Percent refuse crop) *
                (Pollination dependence crop) * nutrient) * (ha / pixel map))

    Returns:
        None.

    """
    crop_nutrient_df = pandas.read_csv(crop_nutrient_df_path)
    yield_raster_path_list = []
    pollination_yield_factor_list = []
    for _, row in crop_nutrient_df.iterrows():
        yield_raster_path = os.path.join(
            yield_raster_dir, f"{row['filenm']}_yield_map.tif")
        if os.path.exists(yield_raster_path):
            yield_raster_path_list.append(yield_raster_path)
            pollination_yield_factor_list.append(
                (1. - row['Percent refuse'] / 100.) * row[nutrient_name])
            if consider_pollination:
                pollination_yield_factor_list[-1] *= (
                    row['Pollination dependence'])
        else:
            raise ValueError(f"not found {yield_raster_path}")

    sample_target_raster_info = pygeoprocessing.get_raster_info(
        sample_target_raster_path)

    yield_raster_info = pygeoprocessing.get_raster_info(
        yield_raster_path_list[0])
    yield_nodata = yield_raster_info['nodata'][0]

    pygeoprocessing.raster_calculator(
        [(yield_nodata, 'raw'), (pollination_yield_factor_list, 'raw')] +
        [(x, 1) for x in yield_raster_path_list], total_yield_op,
        target_10km_yield_path, gdal.GDT_Float32, yield_nodata),

    y_lat_array = numpy.linspace(
        sample_target_raster_info['geotransform'][3],
        sample_target_raster_info['geotransform'][3] +
        sample_target_raster_info['geotransform'][5] *
        sample_target_raster_info['raster_size'][1],
        sample_target_raster_info['raster_size'][1])

    y_ha_array = area_of_pixel(
        abs(sample_target_raster_info['geotransform'][1]),
        y_lat_array) / 10000.0
    y_ha_column = y_ha_array.reshape((y_ha_array.size, 1))

    pygeoprocessing.warp_raster(
        target_10km_yield_path,
        sample_target_raster_info['pixel_size'], target_10s_yield_path,
        'cubicspline', target_bb=sample_target_raster_info['bounding_box'],
        n_threads=2)

    pygeoprocessing.raster_calculator(
        [(target_10s_yield_path, 1), y_ha_column,
         (yield_nodata, 'raw')], density_to_value_op,
        target_10s_production_path, gdal.GDT_Float32, yield_nodata)


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


def _mult_raster_op(array_a, array_b, nodata_a, nodata_b, target_nodata):
    """Multiply a by b and skip nodata."""
    result = numpy.empty(array_a.shape, dtype=numpy.float32)
    result[:] = target_nodata
    valid_mask = (array_a != nodata_a) & (array_b != nodata_b)
    result[valid_mask] = array_a[valid_mask] * array_b[valid_mask]
    return result


def mult_rasters(raster_a_path, raster_b_path, target_path):
    """Multiply a by b and skip nodata."""
    nodata_a = pygeoprocessing.get_raster_info(raster_a_path)['nodata'][0]
    nodata_b = pygeoprocessing.get_raster_info(raster_b_path)['nodata'][0]

    target_nodata = -1.0
    pygeoprocessing.raster_calculator(
        [(raster_a_path, 1), (raster_b_path, 1), (nodata_a, 'raw'),
         (nodata_b, 'raw'), (target_nodata, 'raw')], _mult_raster_op,
        target_path, gdal.GDT_Float32, target_nodata)


if __name__ == '__main__':
    main()
