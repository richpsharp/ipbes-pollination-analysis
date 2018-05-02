"""
Pollination analysis for IPBES.

    From "IPBES Methods: Pollination Contribution to Human Nutrition."
"""
import pandas

import reproduce

WORKING_DIR = '.'

def main():
    """Entry point."""
    reproduce_env = reproduce.Reproduce(WORKING_DIR)

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
        reproduce.url_fetcher(
            'https://storage.googleapis.com/ecoshard-root/'
            'crop_nutrient_md5_d6e67fd79ef95ab2dd44ca3432e9bb4d.csv'))

    crop_nutrient_df = pandas.read_csv(reproduce_env['crop_nutrient_table'])
    print crop_nutrient_df

    # TODO: put in a link to the data crop files here (section 1.2):

if __name__ == '__main__':
    main()
