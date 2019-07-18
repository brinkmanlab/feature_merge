#!/usr/bin/env python
import sys

from . import get_args, load_data, merge_all

def main():
    paths, merge_strategy, merge_order, *args = get_args(sys.argv[1:])

    try:
        db = load_data(paths, merge_strategy)
    except ValueError as e:
        # Catch empty data, exit normally
        print(e, file=sys.stderr)
        exit(0)

    merge_all(db, merge_order, *args)

    # Output header
    print("##gff-version 3")

    for feature in db.all_features(order_by=merge_order):
        # Set source
        if ',' in feature.source:
            feature.attributes["sources"] = feature.source.split(',')
            feature.source = "feature_merge"

        print(feature)

if __name__ == "__main__":
    main()

