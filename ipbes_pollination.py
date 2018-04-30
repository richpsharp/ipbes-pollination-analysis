"""
Pollination analysis for IPBES.

    From "IPBES Methods: Pollination Contribution to Human Nutrition."
"""
import pandas


def main():
    crop_nutrient_df = pandas.read_csv('./pollination_data/crop_nutrient.csv')
    print crop_nutrient_df

if __name__ == '__main__':
    main()
