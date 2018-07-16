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
import tempfile

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
GOOGLE_BUCKET_KEY_PATH = "ecoshard-202992-key.json"
NODATA = -9999
N_WORKERS = 4
DELAYED_START = N_WORKERS >= 0


def main():
    """Entry point."""
    reproduce_env = reproduce.Reproduce(WORKING_DIR)
    LOGGER.debug(reproduce_env['CACHE_DIR'])
    task_graph = taskgraph.TaskGraph(
        reproduce_env['CACHE_DIR'], N_WORKERS, delayed_start=DELAYED_START,
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
        reproduce_env['DATA_DIR'], 'pollination_data',
        os.path.basename(crop_nutrient_url))

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
        reproduce_env['DATA_DIR'], 'observed_yields',
        os.path.basename(yield_zip_url))

    yield_zip_fetch_task = task_graph.add_task(
        func=google_bucket_fetch_and_validate,
        args=(
            yield_zip_url, GOOGLE_BUCKET_KEY_PATH,
            yield_zip_path),
        target_path_list=[yield_zip_path],
        task_name=f'fetch {os.path.basename(yield_zip_path)}',
        skip_if_target_exists=True)

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
        skip_if_target_exists=True)

    # fetch a landcover map to use as a base for the dimensions of the
    # production raster
    landcover_key, (landcover_url, _) = next(iter(landcover_data.items()))
    landcover_env_path = f'landcover/{os.path.basename(landcover_url)}'
    landcover_path = reproduce_env.predict_path(landcover_env_path)
    landcover_fetch_task = task_graph.add_task(
        func=google_bucket_fetch_and_validate,
        args=(landcover_url, GOOGLE_BUCKET_KEY_PATH, landcover_path),
        target_path_list=[landcover_path],
        task_name=f'fetch {landcover_key}',
        skip_if_target_exists=True)

    yield_raster_dir = os.path.join(
        os.path.dirname(yield_zip_path), 'monfreda_2008_observed_yield')

    nut_task_path_tot_prod_1d_map = {}
    nut_task_path_poll_dep_prod_map = {}
    for nutrient_id, nutrient_name in [
            ('en', 'Energy'), ('va', 'VitA'), ('fo', 'Folate')]:
        # total annual production of nutrient
        tot_yield_nut_10km_path = os.path.join(
            WORKING_DIR, 'total_nutrient_rasters',
            f'tot_yield_{nutrient_id}_10km.tif')
        tot_yield_nut_10s_path = os.path.join(
            WORKING_DIR, 'total_nutrient_rasters',
            f'tot_yield_{nutrient_id}_10s.tif')
        tot_prod_nut_10s_path = os.path.join(
            WORKING_DIR, 'total_nutrient_rasters',
            f'tot_prod_{nutrient_id}_10s.tif')

        nut_prod_task = task_graph.add_task(
            func=create_prod_nutrient_raster,
            args=(
                 crop_nutrient_table_path, nutrient_name,
                 yield_raster_dir, False, landcover_path,
                 tot_yield_nut_10km_path, tot_yield_nut_10s_path,
                 tot_prod_nut_10s_path),
            target_path_list=[
                tot_yield_nut_10km_path, tot_yield_nut_10s_path,
                tot_prod_nut_10s_path],
            dependent_task_list=[landcover_fetch_task],
            task_name=f"""create prod raster {
                os.path.basename(tot_prod_nut_10s_path)}""")

        nut_task_path_tot_prod_1d_map[nutrient_id] = (
            nut_prod_task, tot_prod_nut_10s_path)

        # pollination-dependent annual production of nutrient
        poll_dep_yield_nut_10km_path = os.path.join(
            WORKING_DIR, 'total_poll_dependent_production',
            f'poll_dep_yield_{nutrient_id}_10km.tif')
        poll_dep_yield_nut_10s_path = os.path.join(
            WORKING_DIR, 'total_poll_dependent_production',
            f'poll_dep_yield_{nutrient_id}_10s.tif')
        poll_dep_prod_nut_10s_path = os.path.join(
            WORKING_DIR, 'total_poll_dependent_production',
            f'poll_dep_prod_{nutrient_id}_10s.tif')
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
            dependent_task_list=[landcover_fetch_task],
            task_name=f"""create poll dep production raster {
                os.path.basename(poll_dep_prod_nut_10s_path)}""")
        nut_task_path_poll_dep_prod_map[nutrient_id] = (
            pol_dep_prod_task, poll_dep_prod_nut_10s_path)

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

    # mask landcover into agriculture and pollinator habitat
    for landcover_key, (landcover_url, landcover_short_suffix) in (
            landcover_data.items()):
        landcover_env_path = f'landcover/{os.path.basename(landcover_url)}'
        landcover_path = reproduce_env.predict_path(landcover_env_path)
        landcover_fetch_task = task_graph.add_task(
            func=google_bucket_fetch_and_validate,
            args=(landcover_url, GOOGLE_BUCKET_KEY_PATH, landcover_path),
            target_path_list=[landcover_path],
            task_name=f'fetch {landcover_key}',
            skip_if_target_exists=True)

        _ = task_graph.add_task(
            func=build_overviews,
            priority=-100,
            args=(landcover_path, 'mode'),
            target_path_list=[f'{landcover_path}.ovr'],
            dependent_task_list=[landcover_fetch_task],
            task_name=f'compress {os.path.basename(landcover_path)}',
            skip_if_target_exists=True)

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
                priority=-100,
                args=(mask_target_path, 'mode'),
                target_path_list=[
                    f'{mask_target_path}.ovr'],
                dependent_task_list=[mask_task],
                task_name=f'compress {os.path.basename(mask_target_path)}')

            if mask_prefix == 'hab':
                hab_task_path_tuple = (mask_task, mask_target_path)
            elif mask_prefix == 'ag':
                ag_task_path_tuple = (mask_task, mask_target_path)

        proportional_hab_area_2km_path = os.path.join(
            reproduce_env['DATA_DIR'], 'proportional_area',
            f'{landcover_key}_hab_prop_area_2km.tif')
        prop_hab_area_2km_task = task_graph.add_task(
            func=pygeoprocessing.convolve_2d,
            args=[
                (hab_task_path_tuple[1], 1), (kernel_raster_path, 1),
                proportional_hab_area_2km_path],
            kwargs={
                'working_dir': reproduce_env['CACHE_DIR'],
                'gtiff_creation_options': (
                    'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=DEFLATE',
                    'PREDICTOR=3', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256',
                    'NUM_THREADS=ALL_CPUS')},
            dependent_task_list=[hab_task_path_tuple[0], kernel_task],
            target_path_list=[proportional_hab_area_2km_path],
            task_name=(
                'calculate proportional'
                f' {os.path.basename(proportional_hab_area_2km_path)}'),
            skip_if_target_exists=True)

        _ = task_graph.add_task(
            func=build_overviews,
            priority=-100,
            args=(proportional_hab_area_2km_path, 'average'),
            target_path_list=[
                f'{proportional_hab_area_2km_path}.ovr'],
            dependent_task_list=[prop_hab_area_2km_task],
            task_name=(
                'compress '
                f'{os.path.basename(proportional_hab_area_2km_path)}'),
            skip_if_target_exists=True)

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
        poll_suff_path = os.path.join(
            WORKING_DIR, 'main_outputs',
            f'poll_suff_10s_{landcover_short_suffix}.tif')
        poll_suff_task = task_graph.add_task(
            func=threshold_select_raster,
            args=(
                proportional_hab_area_2km_path,
                ag_task_path_tuple[1], threshold_val,
                poll_suff_path),
            target_path_list=[poll_suff_path],
            dependent_task_list=[
                prop_hab_area_2km_task, ag_task_path_tuple[0]],
            task_name=f'threshold {os.path.basename(poll_suff_path)}')

        _ = task_graph.add_task(
            func=build_overviews,
            priority=-100,
            args=(poll_suff_path, 'mode'),
            target_path_list=[
                f'{poll_suff_path}.ovr'],
            dependent_task_list=[poll_suff_task],
            task_name=f'compress {os.path.basename(poll_suff_path)}',
            )

        # tot_prod_en|va|fo_10s|1d_cur|ssp1|ssp3|ssp5
        # total annual production of energy (KJ/yr), vitamin A (IU/yr),
        # and folate (mg/yr)
        cont_prod_nutrient_task_path_list = []
        for nutrient_id in ['en', 'va', 'fo']:
            tot_prod_task, tot_prod_1d_path = (
                nut_task_path_tot_prod_1d_map[nutrient_id])

            tot_prod_nut_scenario_path = os.path.join(
                WORKING_DIR, 'main_outputs',
                f'tot_prod_{nutrient_id}_10s_{landcover_short_suffix}.tif')

            tot_prod_nut_scenario_task = task_graph.add_task(
                func=mult_rasters,
                args=(
                    ag_task_path_tuple[1], tot_prod_1d_path,
                    tot_prod_nut_scenario_path),
                target_path_list=[tot_prod_nut_scenario_path],
                dependent_task_list=[tot_prod_task, ag_task_path_tuple[0]],
                task_name=(
                    f'tot_prod_{nutrient_id}_10s_{landcover_short_suffix}'))

            _ = task_graph.add_task(
                func=build_overviews,
                priority=-100,
                args=(tot_prod_nut_scenario_path, 'average'),
                target_path_list=[f'{tot_prod_nut_scenario_path}.ovr'],
                dependent_task_list=[tot_prod_nut_scenario_task],
                task_name=('compress' +
                           os.path.basename(tot_prod_nut_scenario_path)))

            poll_dep_prod_task, poll_dep_prod_path = (
                nut_task_path_poll_dep_prod_map[nutrient_id])

            poll_dep_prod_nut_scenario_path = os.path.join(
                WORKING_DIR, 'main_outputs',
                f'poll_dep_prod_{nutrient_id}_10s_'
                f'{landcover_short_suffix}.tif')

            poll_dep_prod_nut_scenario_task = task_graph.add_task(
                func=mult_rasters,
                args=(
                    ag_task_path_tuple[1], poll_dep_prod_path,
                    poll_dep_prod_nut_scenario_path),
                target_path_list=[poll_dep_prod_nut_scenario_path],
                dependent_task_list=[
                    poll_dep_prod_task, ag_task_path_tuple[0]],
                task_name=(
                    f'poll_dep_prod_{nutrient_id}_'
                    f'10s_{landcover_short_suffix}'))

            _ = task_graph.add_task(
                func=build_overviews,
                priority=-100,
                args=(poll_dep_prod_nut_scenario_path, 'average'),
                target_path_list=[f'{poll_dep_prod_nut_scenario_path}.ovr'],
                dependent_task_list=[poll_dep_prod_nut_scenario_task],
                task_name=('compress' +
                           os.path.basename(poll_dep_prod_nut_scenario_path)))

            # poll_serv_en|va|fo_10s|1d_cur|ssp1|ssp3|ssp5: wild pollination-
            # derived annual production of energy (KJ/yr), vitamin A (IU/yr),
            # and folate (mg/yr) (given the habitat and whether it is
            # sufficient to meet pollination needs)
            poll_serv_prod_nut_scenario_path = os.path.join(
                WORKING_DIR, 'main_outputs',
                f'poll_serv_prod_{nutrient_id}_10s_'
                f'{landcover_short_suffix}.tif')

            poll_serv_prod_nut_scenario_task = task_graph.add_task(
                func=mult_rasters,
                args=(
                    poll_dep_prod_nut_scenario_path, poll_suff_path,
                    poll_serv_prod_nut_scenario_path),
                target_path_list=[poll_serv_prod_nut_scenario_path],
                dependent_task_list=[
                    poll_suff_task, poll_dep_prod_nut_scenario_task],
                task_name=(
                    f'poll_serv_prod_{nutrient_id}_'
                    f'10s_{landcover_short_suffix}'))

            _ = task_graph.add_task(
                func=build_overviews,
                priority=-100,
                args=(poll_serv_prod_nut_scenario_path, 'average'),
                target_path_list=[f'{poll_serv_prod_nut_scenario_path}.ovr'],
                dependent_task_list=[poll_serv_prod_nut_scenario_task],
                task_name=('compress' + os.path.basename(
                    poll_serv_prod_nut_scenario_path)))

            # cont_prod_en|va|fo_10s|1d_cur|ssp1|ssp3|ssp5
            # contribution of wild pollination to total annual micronutrient
            # production, as a proportion of total energy, vitamin or folate
            target_cont_prod_nutrient_path = os.path.join(
                WORKING_DIR, 'main_outputs',
                f'cont_poll_serv_prod_{nutrient_id}_10s_'
                f'{landcover_short_suffix}.tif')
            cont_prod_task = task_graph.add_task(
                func=create_cont_prod_nutrient_raster,
                args=(
                    poll_serv_prod_nut_scenario_path,
                    tot_prod_nut_scenario_path,
                    target_cont_prod_nutrient_path),
                target_path_list=[target_cont_prod_nutrient_path],
                dependent_task_list=[
                    poll_serv_prod_nut_scenario_task,
                    tot_prod_nut_scenario_task],
                task_name=(
                    f'cont_poll_serv_prod_{nutrient_id}_10s_'
                    f'{landcover_short_suffix}.tif'))
            cont_prod_nutrient_task_path_list.append(
                (cont_prod_task, target_cont_prod_nutrient_path))

            _ = task_graph.add_task(
                func=build_overviews,
                priority=-100,
                args=(target_cont_prod_nutrient_path, 'average'),
                target_path_list=[f'{target_cont_prod_nutrient_path}.ovr'],
                dependent_task_list=[cont_prod_task],
                task_name=('compress' + os.path.basename(
                    target_cont_prod_nutrient_path)))

        # cont_poll_serv_prod_avg_10s|1d_cur|ssp1|ssp3|ssp5
        # average contribution of wild pollination to total annual
        # micronutrient production, across all three nutrients
        cont_poll_serv_prod_avg_10s_path = os.path.join(
            WORKING_DIR, 'main_outputs',
            f'cont_poll_serv_prod_avg_10s_{landcover_short_suffix}.tif')
        cont_poll_serv_prod_avg_task = task_graph.add_task(
            func=create_avg_raster,
            args=(
                [task_path_tuple[1] for task_path_tuple in
                 cont_prod_nutrient_task_path_list], cont_poll_serv_prod_avg_10s_path),
            target_path_list=[cont_poll_serv_prod_avg_10s_path],
            dependent_task_list=(
                [task_path_tuple[0] for task_path_tuple in
                 cont_prod_nutrient_task_path_list]),
            task_name=f'cont_poll_serv_prod_avg_10s_{landcover_short_suffix}')

        _ = task_graph.add_task(
            func=build_overviews,
            priority=-100,
            args=(cont_poll_serv_prod_avg_10s_path, 'average'),
            target_path_list=[f'{cont_poll_serv_prod_avg_10s_path}.ovr'],
            dependent_task_list=[cont_poll_serv_prod_avg_task],
            task_name=('compress' + os.path.basename(
                cont_poll_serv_prod_avg_10s_path)))

        # c_cont_poll_serv_prod_avg_10s|1d_ssp1|ssp3|ssp5

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
        reproduce_env['DATA_DIR'], 'spatial_population_scenarios',
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
        skip_if_target_exists=True,
        priority=100)

    gpw_urls = {
        'gpw_v4_e_a000_014bt_2010_dens': 'https://storage.cloud.google.com/ecoshard-root/ipbes/gpw_v4_e_a000_014bt_2010_dens_30_sec_md5_f4129fcbe7a2f65fddfb378a93cc4d67.tif',
        'gpw_v4_e_a000_014ft_2010_dens': 'https://storage.cloud.google.com/ecoshard-root/ipbes/gpw_v4_e_a000_014ft_2010_dens_30_sec_md5_8b2871967b71534d56d2df83e413bf33.tif',
        'gpw_v4_e_a000_014mt_2010_dens': 'https://storage.cloud.google.com/ecoshard-root/ipbes/gpw_v4_e_a000_014mt_2010_dens_30_sec_md5_8ddf234eabc0025efd5678776e2ae792.tif',
        'gpw_v4_e_a065plusbt_2010_dens': 'https://storage.cloud.google.com/ecoshard-root/ipbes/gpw_v4_e_a065plusbt_2010_dens_30_sec_md5_87a734bf80b87aa773734122471e2cf1.tif',
        'gpw_v4_e_a065plusmt_2010_dens': 'https://storage.cloud.google.com/ecoshard-root/ipbes/gpw_v4_e_a065plusmt_2010_dens_30_sec_md5_1d36e79aa083ee25de7295a4a7d10a4b.tif',
        'gpw_v4_e_atotpopbt_2010_dens': 'https://storage.cloud.google.com/ecoshard-root/ipbes/gpw_v4_e_atotpopbt_2010_dens_30_sec_md5_3202da812b3a98bb6df1440eaa28fead.tif',
        'gpw_v4_e_atotpopft_2010_dens': 'https://storage.cloud.google.com/ecoshard-root/ipbes/gpw_v4_e_atotpopft_2010_dens_30_sec_md5_5bbb72a050c76264e1b6a3c7530fedce.tif',
        'gpw_v4_e_atotpopmt_2010_dens': 'https://storage.cloud.google.com/ecoshard-root/ipbes/gpw_v4_e_atotpopmt_2010_dens_30_sec_md5_31637ca784b8b917222661d4a915ead6.tif',
    }

    for gpw_id, gpw_url in gpw_urls.items():
        gpw_dens_path = os.path.join(
            reproduce_env['DATA_DIR'], 'gpw_pop_densities',
            os.path.basename(gpw_url))
        gpw_fetch_task = task_graph.add_task(
            func=google_bucket_fetch_and_validate,
            args=(gpw_url, GOOGLE_BUCKET_KEY_PATH, gpw_dens_path),
            target_path_list=[gpw_dens_path],
            task_name=f"""fetch {os.path.basename(gpw_dens_path)}""",
            priority=100)

        gpw_count_path = os.path.join(
            WORKING_DIR, 'gpw_count', f"""{gpw_id[:-4]}count.tif""")
        gpw_count_task = task_graph.add_task(
            func=calc_pop_count,
            args=(gpw_dens_path, gpw_count_path),
            target_path_list=[gpw_count_path],
            dependent_task_list=[gpw_fetch_task],
            task_name=f"""pop count {os.path.basename(gpw_count_path)}""")

    # calculate the ratios of the spatial scenarios from cur to ssp1..5
    # then multiply by the population counts to get ssp1..5 counts

    # calculate gpw count

    task_graph.close()
    task_graph.join()
    # END MAIN


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
        [(gpw_dens_path, 1), y_ha_column, nodata],
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
        result = numpy.empty_like(select_array)
        result[:] = target_nodata
        valid_mask = base_array != base_nodata
        result[valid_mask] = numpy.where(
            base_array[valid_mask] >= threshold_val,
            select_array[valid_mask], 0)
        return result

    pygeoprocessing.raster_calculator(
        [(base_raster_path, 1), (select_raster_path, 1),
         threshold_val, base_nodata, target_nodata], threshold_select_op,
        target_path, gdal.GDT_Byte, target_nodata, gtiff_creation_options=(
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
    result = numpy.empty_like(crop_yield_array_list[0])
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
    result = numpy.empty_like(density_array)
    result[:] = 0.0
    valid_mask = density_array != density_nodata
    result[valid_mask] = density_array[valid_mask] * y_ha_array[valid_mask]
    return result


def create_prod_nutrient_raster(
        crop_nutrient_df_path, nutrient_name, yield_raster_dir,
        consider_pollination, sample_target_raster_path,
        target_10km_yield_path, target_10s_yield_path,
        target_10s_production_path):
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
        [yield_nodata, (pollination_yield_factor_list, 'raw')] +
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
        'bilinear', target_bb=sample_target_raster_info['bounding_box'])

    pygeoprocessing.raster_calculator(
        [(target_10s_yield_path, 1), y_ha_column, yield_nodata],
        density_to_value_op, target_10s_production_path, gdal.GDT_Float32,
        yield_nodata)


def create_cont_prod_nutrient_raster(
        poll_serv_nutrient_path, tot_prod_nut_path,
        target_cont_prod_nutrient_path):
    """Calculate proportion of poll_serv to tot_prod.

    Contribution of wild pollination to total annual micronutrient production,
    as a proportion of nutrient.

    Parameters:
        poll_serv_nutrient_path (string): path to pollination service
            raster (portion of nutrient production that is pollinator
            dependent).
        tot_prod_nut_path (string): path to total nutrient production.
        target_cont_prod_nutrient_path (string): path to target file that
            will calculate the proportion of pollinator production to total
            production. If any production is 0, the value will be 0.

    Returns:
        None.

    """
    poll_serv_nodata = pygeoprocessing.get_raster_info(
        poll_serv_nutrient_path)['nodata'][0]
    tot_prod_nodata = pygeoprocessing.get_raster_info(
        tot_prod_nut_path)['nodata'][0]
    target_nodata = -1.

    def calc_proportion(poll_serv_array, tot_prod_array):
        result = numpy.empty_like(poll_serv_array)
        result[:] = target_nodata
        zero_mask = tot_prod_array == 0
        valid_mask = (
            tot_prod_array != tot_prod_nodata) & (
            poll_serv_array != poll_serv_nodata) & ~zero_mask
        result[zero_mask] = 0.0
        result[valid_mask] = (
            poll_serv_array[valid_mask] / tot_prod_array[valid_mask])
        return result

    pygeoprocessing.raster_calculator(
        [(poll_serv_nutrient_path, 1), (tot_prod_nut_path, 1)],
        calc_proportion, target_cont_prod_nutrient_path, gdal.GDT_Float32,
        target_nodata)


def create_avg_raster(base_path_list, target_avg_path):
    """Calculate the per-pixel average of the base path list.

    Parameters:
        base_path_list (list): list of raster paths to average.
        target_avg_path (str): path to desired target file that will have a
            per-pixel average of all the input rasters.

    Returns:
        None

    """
    nodata_set = set([
        pygeoprocessing.get_raster_info(path)['nodata'][0]
        for path in base_path_list])
    assert(len(nodata_set) == 1)
    nodata = nodata_set.pop()

    def average_op(*array_list):
        array_stack = numpy.array(array_list)
        nodata_mask = array_stack == nodata
        array_stack[nodata_mask] = 0
        result = numpy.sum(array_stack, axis=0)
        valid_count = numpy.sum(~nodata_mask, axis=0)
        valid_mask = valid_count > 0
        result[valid_mask] /= valid_count[valid_mask]
        result[~valid_mask] = nodata
        return result

    pygeoprocessing.raster_calculator(
        [(path, 1) for path in base_path_list], average_op,
        target_avg_path, gdal.GDT_Float32, nodata)


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
    result = numpy.empty_like(array_a)
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
        [(raster_a_path, 1), (raster_b_path, 1), nodata_a, nodata_b,
         target_nodata], _mult_raster_op, target_path, gdal.GDT_Float32,
        target_nodata)


if __name__ == '__main__':
    main()
