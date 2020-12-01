#!/usr/bin/env python
import sys

from . import get_args, load_data, merge_all, update

def main():
    paths, merge_strategy, merge_order, *args = get_args(sys.argv[1:])

    try:
        db = load_data(paths, merge_strategy)
    except ValueError as e:
        # Catch empty data, exit normally
        print(e, file=sys.stderr)
        exit(0)

    merged_features = merge_all(db, merge_order, *args)

    for feature in merged_features:
        feature.attributes["sources"] = feature.source.split(',')
        feature.source = "feature_merge"

    update(db, merged_features, make_backup=False, merge_strategy='replace')

    # Output header
    print("##gff-version 3")

    for feature in db.all_features(order_by=merge_order):
        print(feature)

if __name__ == "__main__":
    main()

