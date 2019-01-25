#!/usr/bin/env python
import gffutils
import sys
import getopt

usage = "Usage: feature_merge.py [-i] [-e] [-v] [-f type[,type..]].. <input1> <input2> [..<input_n>]. Accepts GFF or GTF format.\n" \
        "-v Print version and exit\n" \
        "-f Comma seperated types of features to merge. Must be terms or accessions from the SOFA sequence ontology or \"ALL\". (Can be provided more than once to specify multiple merge groups)\n" \
        "-i Ignore strand, merge feature regardless of strand\n" \
        "-e Exclude component features from output"

def merge(target_db, features, ignore_strand=False, ignore_featuretype=False):
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
    the newly merged features and a list of the features that compose it.
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
    feature_components = [features[0]]

    # We don't need to check the first one, so start at feature #2.
    for feature in features[1:]:
        # Does this feature start within the currently merged feature?...
        if feature.start <= current_merged_stop + 1 and feature.seqid == seqid and (ignore_strand or feature.strand == strand) and (ignore_featuretype or feature.featuretype == featuretype):
            feature_components.append(feature)
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
            yield target_db._feature_returner(
                seqid=seqid,
                source='',
                featuretype="sequence_feature" if len(features) > 1 else featuretype,
                start=current_merged_start,
                end=current_merged_stop,
                score='.',
                strand='.' if ignore_strand else strand,
                frame='.',
                attributes=''), feature_components

            # and we start a new one, initializing with this feature's
            # start and stop.
            current_merged_start = feature.start
            current_merged_stop = feature.stop
            seqid = feature.seqid
            strand = feature.strand
            featuretype = feature.featuretype
            feature_components = [feature]

    # need to yield the last one.
    if len(features) == 1:
        feature = features[0]
    yield target_db._feature_returner(
        seqid=seqid,
        source='',
        featuretype="sequence_feature" if len(features) > 1 else featuretype,
        start=current_merged_start,
        end=current_merged_stop,
        score='.',
        strand='.' if ignore_strand else strand,
        frame='.',
        attributes=''), feature_components

if __name__ == '__main__':
    ignore_strand = False
    ignore_featuretypes = False
    exclude_components = False
    featuretypes_groups = []

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'vief:')
        for opt, val in opts:
            if opt == '-v':
                import __version
                print(__version.__versionstr__)
                exit(0)
            elif opt == '-i':
                ignore_strand = True
            elif opt == '-e':
                exclude_components = True
            elif opt == '-f':
                ignore_featuretypes = True
                if val != "ALL": featuretypes_groups.append(tuple(filter(None, val.split(','))))

    except getopt.GetoptError as err:
        args = []

    if len(args) < 1:
        print(usage, file=sys.stderr)
        exit(1)

    if not featuretypes_groups:
        featuretypes_groups.append(None)

    if ignore_featuretypes:
        merge_order = ('seqid', 'strand', 'start', 'featuretype')
    else:
        merge_order = ('seqid', 'featuretype', 'strand', 'start')

    #Load input data
    db = gffutils.create_db(args[0], ":memory:", merge_strategy="merge")
    for input in args[2:]:
        db = db.update(input, make_backup=False)

    remaining_featuretypes = set(db.featuretypes())

    #Merge features per featuregroup
    for featuregroup in featuretypes_groups:
        if featuregroup: remaining_featuretypes -= featuregroup
        else: remaining_featuretypes = set()
        for merged, components in merge(db, db.all_features(featuretype=featuregroup, order_by=merge_order), ignore_strand, ignore_featuretypes):
            merged.id = "merged"
            if len(components) > 1:
                for component in components:
                    if not component.id:
                        component.id = hex(hash(component) + 2**63)[2:]
                    merged.id += "-" + component.id

                if len(merged.id) > 32:
                    merged.id = hex(hash(merged.id)+ 2**63)[2:]

                merged.attributes["ID"] = merged.id

                for component in components:
                    if "Parent" not in component.attributes:
                        component.attributes["Parent"] = []
                    component.attributes["Parent"].append(merged.id)

            if len(components) == 1 or not exclude_components:
                for component in components:
                    print(component)

            if len(components) > 1: # Don't output merged record of single record
                # Eliminate duplicate sources before adding to merged.source
                merged.source = ",".join(set(component.source for component in components))
                print(merged)

    #Output any features that may have not been in the -f arguments
    if remaining_featuretypes:
        for feature in db.all_features(featuretype=remaining_featuretypes, order_by=merge_order):
            print(feature)