"""
Pollination analysis for IPBES.

    From "IPBES Methods: Pollination Contribution to Human Nutrition."
"""
import pandas


def main():
    """Entry point."""

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
    # found at (permanent link to table).
    crop_nutrient_df = pandas.read_csv('./pollination_data/crop_nutrient.csv')
    print crop_nutrient_df

    # TODO: put in a link to the data crop files here (section 1.2):

if __name__ == '__main__':
    main()
