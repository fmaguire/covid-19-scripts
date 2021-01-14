#!/usr/bin/env python

import argparse
import pgeocode
from pathlib import Path
import pandas as pd

def check_file(path: str) -> Path:
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")


def add_canada_lat_long(metadata, lat_long):
    """
    Use pgeocode to add lat and longs for canadian locations
    """
    # get all locations in lat long file to check which need added
    loc_lat_longs = set(lat_longs.query('geo_scale=="location"')['geo_loc'].unique())

    # get canada locs that are postcodes and not in lat long already
    canada = metadata.query('country=="Canada"')
    postcodes = canada.loc[canada['location'].str.match("^[a-zA-Z][0-9][a-zA-Z]$").fillna(False), 'location'].unique()
    postcodes = set(postcodes) - loc_lat_longs

    # get geocode database
    geocode = pgeocode.Nominatim('CA')
    postcode_lat_long = {'geo_scale': [],
                         'geo_loc': [],
                         'lat': [],
                         'long': []}
    for postcode in postcodes:
        geoloc_for_postcode = geocode.query_postal_code(postcode)
        postcode_lat_long['geo_scale'].append('location')
        postcode_lat_long['geo_loc'].append(postcode)
        postcode_lat_long['lat'].append(geoloc_for_postcode['latitude'])
        postcode_lat_long['long'].append(geoloc_for_postcode['longitude'])

    # update lat_longs with postcode locations
    updated_lat_longs = lat_longs.append(pd.DataFrame(postcode_lat_long))
    return updated_lat_longs


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Update ncov lat/long file with "
                                     "geoloc data for Canadian postcodes in "
                                     "nextmeta location field")
    parser.add_argument("-m", "--metadata", required=True, type=check_file,
                         help="nextmeta formatted metadata file")
    parser.add_argument("-l", "--lat_longs", required=True, type=check_file,
                         help="ncov lat_long file e.g., ncov/defaults/lat_longs.tsv")
    parser.add_argument("-o", "--output", default="custom_lat_longs.tsv",
                        help="Path to write updated lat_longs file")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep='\t')
    lat_longs = pd.read_csv(args.lat_longs, sep='\t', names=["geo_scale",
                                                             "geo_loc",
                                                             "lat",
                                                             "long"])

    postcode_added_lat_longs = add_canada_lat_long(metadata, lat_longs)

    postcode_added_lat_longs.to_csv(args.output, sep='\t', header=False,
                                    index=False)
