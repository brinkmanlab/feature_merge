#!/usr/bin/env python
import gffutils
import sys
import getopt

usage = "Usage: feature_merge.py [-i] [-e] <input1> <input2> [...<input_n>]. Accepts GFF or GTF format.\n" \
        "-i Ignore strand, merge feature regardless of strand.\n" \
        "-e Exclude component features from output."

def merge(target_db, features, ignore_strand=False):
    """
    Merge overlapping features together.

    Parameters
    ----------

    features : iterator of Feature instances

    ignore_strand : bool
        If True, features on multiple strands will be merged, and the final
        strand will be set to '.'.  Otherwise, ValueError will be raised if
        trying to merge features on differnt strands.

    Returns
    -------
    A generator object that yields :class:`Feature` objects representing
    the newly merged features.
    """

    # Consume iterator up front...
    features = list(features)

    if len(features) == 0:
        raise StopIteration

    # To start, we create a merged feature of just the first feature.
    current_merged_start = features[0].start
    current_merged_stop = features[0].stop
    seqid = features[0].seqid
    strand = features[0].strand
    featuretype = features[0].featuretype
    feature_components = []

    # We don't need to check the first one, so start at feature #2.
    for feature in features[1:]:
        # Does this feature start within the currently merged feature?...
        if feature.start <= current_merged_stop + 1 and feature.seqid == seqid and (ignore_strand or feature.strand == strand):
            feature_components += feature
            # ...It starts within, so leave current_merged_start where it
            # is.  Does it extend any farther?
            if feature.stop >= current_merged_stop:
                # Extends further, so set a new stop position
                current_merged_stop = feature.stop
            else:
                # If feature.stop < current_merged_stop, it's completely
                # within the previous feature.  Nothing more to do.
                continue
        else:
            # The start position is outside the merged feature, so we're
            # done with the current merged feature.  Prepare for output...
            merged_feature = dict(
                seqid=seqid,
                source='',
                featuretype=featuretype,
                start=current_merged_start,
                end=current_merged_stop,
                score='.',
                strand='.' if ignore_strand else strand,
                frame='.',
                attributes='')
            yield target_db._feature_returner(**merged_feature), feature_components

            # and we start a new one, initializing with this feature's
            # start and stop.
            current_merged_start = feature.start
            current_merged_stop = feature.stop
            seqid = feature.seqid
            strand = feature.strand
            featuretype = feature.featuretype
            feature_components = []

    # need to yield the last one.
    if len(features) == 1:
        feature = features[0]
    merged_feature = dict(
        seqid=feature.seqid,
        source='',
        featuretype=feature.featuretype,
        start=current_merged_start,
        end=current_merged_stop,
        score='.',
        strand=feature.strand,
        frame='.',
        attributes='')
    yield target_db._feature_returner(**merged_feature), [feature]

if __name__ == '__main__':
    ignore_strand = False
    exclude_components = False

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'ie')
        for opt, val in opts:
            if opt == '-i':
                ignore_strand = True
            elif opt == '-e':
                exclude_components = True
    except getopt.GetoptError as err:
        args = []

    if len(args) < 2:
        print(usage, file=sys.stderr)
        exit(1)

    db = gffutils.create_db(args[1], ":memory:", merge_strategy="merge")
    for input in args[2:]:
        db = db.update(input, make_backup=False)

    for merged, components in merge(db, db.all_features(order_by=('seqid', 'strand', 'start')), ignore_strand):
        if len(components) == 1 or not exclude_components:
            for component in components:
                component.attributes['Parent'] = merged.id
                print(component)
        if len(components) > 1: # Don't output merged record of single record
            for component in components:
                merged.source += component.source
            print(merged)
