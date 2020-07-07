"""
Feature Merge - Merge annotations in GFF files.
Can merge overlapping annotations based on criteria.

Run without arguments for options help.

Nolan Woods (nolan_w@sfu.ca) 2019
Brinkman Lab, SFU
"""
import gffutils
import sys
import os
import getopt
from typing import Sequence, Set, Callable

from . import merge_criteria as mc

usage = """
Usage: feature_merge [-i] [-e] [-x] [-v] [-t <number>]  [-m merge|append|error|skip|replace] [-f type[,type..]].. <input1> [<input_n>..]
Accepts GFF or GTF format.
-v Print version and exit
-f Comma seperated types of features to merge. Must be terms or accessions from the SOFA sequence ontology, \"ALL\", or \"NONE\". (Can be provided more than once to specify multiple merge groups)
-i Ignore strand, merge feature regardless of strand
-s Ignore sequence id, merge feature regardless of sequence id
-x Only merge features with identical coordinates
-t Threshold distance between features to merge 
-e Exclude component features from output
-m Merge strategy used to deal with id collisions between input files.
    merge: attributes of all features with the same primary key will be merged
    append: entry will have a unique, autoincremented primary key assigned to it (default)
    error: exception will be raised. This means you will have to edit the file yourself to fix the duplicated IDs
    skip: ignore duplicates, emitting a warning
    replace: keep last duplicate
"""[1:-1]

merge_strategies = {"merge": "merge", "append": "create_unique", "error": "error", "skip": "warning",
                    "replace": "replace"}


# The following two functions are waiting to be integrated into the gffutils library.
# PR https://github.com/daler/gffutils/pull/130 needs to be updated with below implementation
# PR https://github.com/daler/gffutils/pull/131
# PR https://github.com/daler/gffutils/pull/152
# Issue https://github.com/daler/gffutils/issues/119

no_children = tuple()


def _finalize_merge(feature, feature_children):
    if len(feature_children) > 1:
        feature.source = ','.join(set(child.source for child in feature_children))
        feature.children = feature_children
    else:
        feature.children = no_children
    return feature


def merge(self, features, ignore_strand=False,
          merge_criteria: [Callable] = (mc.seqid, mc.overlap_end_inclusive, mc.strand, mc.feature_type),
          multiline: bool = False):
    """
    Merge features matching criteria together

    Criteria callback interface:
    callback(acc: gffutils.Feature, cur: gffutils.Feature, components: [gffutils.Feature]) -> bool
    acc: current accumulated feature
    cur: candidate feature to merge
    components: list of features that compose acc
    return: true to merge cur into acc, false to set cur to acc (start a new merged feature)
    If merge criteria allows different feature types then the merged features feature types should have their
    feature_type property reassigned to a more specific ontology value.

    Returned Features have a special property called 'children' that is a list of the component features.
    This only exists for the lifetime of the Feature instance.

    :param features: Iterable of Feature instances to merge
    :param ignore_strand: DEPRECIATED remove 'strand' from criteria if true
    :param merge_criteria: List of merge criteria callbacks
    :param multiline: True to emit multiple features with the same ID attribute, False otherwise
    :return: generator emitting merged Feature instances
    """
    if not isinstance(merge_criteria, list):
        try:
            merge_criteria = list(merge_criteria)
        except TypeError:
            merge_criteria = [merge_criteria]

    if ignore_strand and mc.strand in merge_criteria:
        merge_criteria.remove(mc.strand)

    # To start, we create a merged feature of just the first feature.
    features = iter(features)
    last_id = None
    current_merged = None
    feature_children = []

    for feature in features:
        if current_merged is None:
            if all(criteria(feature, feature, feature_children) for criteria in merge_criteria):
                current_merged = feature
                feature_children = [feature]
            else:
                yield _finalize_merge(feature, no_children)
                last_id = None
            continue

        if len(feature_children) == 0: # current_merged is last feature and unchecked
            if all(criteria(current_merged, current_merged, feature_children) for criteria in merge_criteria):
                feature_children.append(current_merged)
            else:
                yield _finalize_merge(current_merged, no_children)
                current_merged = feature
                last_id = None
                continue

        if all(criteria(current_merged, feature, feature_children) for criteria in merge_criteria):
            # Criteria satisfied, merge
            # TODO Test multiline records and iron out the following code
            #if multiline and (feature.start > current_merged.end + 1 or feature.end + 1 < current_merged.start):
            #    # Feature is possibly multiline (discontiguous), keep ID but start new record
            #    yield _finalize_merge(current_merged, feature_children)
            #    current_merged = feature
            #    feature_children = [feature]

            if len(feature_children) == 1:
                # Current merged is only child and merge is going to occur, make copy
                current_merged = vars(current_merged).copy()
                del current_merged['attributes']
                del current_merged['extra']
                del current_merged['dialect']
                del current_merged['keep_order']
                del current_merged['sort_attribute_values']
                current_merged = self._feature_returner(**current_merged)
                if not last_id:
                    # Generate unique ID for new Feature
                    # This can be simplified after https://github.com/daler/gffutils/pull/135
                    self._autoincrements[current_merged.featuretype] = \
                        self._autoincrements.get(current_merged.featuretype, 0) + 1
                    last_id = current_merged.featuretype + '_' + str(self._autoincrements[current_merged.featuretype])
                current_merged['ID'] = last_id
                current_merged.id = last_id

            feature_children.append(feature)

            # Merge attributes. Removed as it doesn't make sense to collect attributes in an aggrigate feature when
            # Parent relationships present
            # current_merged.attributes = gffutils.helpers.merge_attributes(feature.attributes, current_merged.attributes)
            # Preserve ID
            # current_merged['ID'] = last_id

            # Set mismatched properties to ambiguous values
            if feature.seqid not in current_merged.seqid.split(','): current_merged.seqid += ',' + feature.seqid
            if feature.strand != current_merged.strand: current_merged.strand = '.'
            if feature.frame != current_merged.frame: current_merged.frame = '.'
            if feature.featuretype != current_merged.featuretype: current_merged.featuretype = "sequence_feature"

            if feature.start < current_merged.start:
                # Extends prior, so set a new start position
                current_merged.start = feature.start

            if feature.end > current_merged.end:
                # Extends further, so set a new stop position
                current_merged.end = feature.end

        else:
            yield _finalize_merge(current_merged, feature_children)
            current_merged = feature
            feature_children = []
            last_id = None

    if current_merged:
        yield _finalize_merge(current_merged, feature_children)


def update(self, data, **kwargs):
    """
    Ripped this out of FeatureDB.update() to deal with a bug.
    Update database with features in `data`.
    If the file is empty, return rather than throw exception.

    data : str, iterable, FeatureDB instance
        If FeatureDB, all data will be used. If string, assume it's
        a filename of a GFF or GTF file.  Otherwise, assume it's an
        iterable of Feature objects.  The classes in gffutils.iterators may
        be helpful in this case.
    """
    from gffutils import create
    from gffutils import iterators

    # Handle all sorts of input
    data = iterators.DataIterator(data)

    if self.dialect['fmt'] == 'gtf':
        if 'id_spec' not in kwargs:
            kwargs['id_spec'] = {'gene': 'gene_id', 'transcript': 'transcript_id'}
        db = create._GTFDBCreator(
            data=data, dbfn=self.dbfn, dialect=self.dialect, **kwargs)
    elif self.dialect['fmt'] == 'gff3':
        if 'id_spec' not in kwargs:
            kwargs['id_spec'] = 'ID'
        db = create._GFFDBCreator(
            data=data, dbfn=self.dbfn, dialect=self.dialect, **kwargs)
    else:
        raise ValueError

    peek, data._iter = iterators.peek(data._iter, 1)
    if len(peek) == 0: return  # If the file is empty then do nothing

    db._autoincrements.update(self._autoincrements)
    db._populate_from_lines(data)
    db._update_relations()
    db._finalize()
    self._autoincrements.update(db._autoincrements)


def assign_child(parent, child):
    """
    Helper for add_relation()
    Sets 'Parent' attribute to parent['ID']
    :param parent: parent Feature
    :param child: child Feature
    :return: child
    """
    child.attributes['Parent'] = parent['ID']
    return child


def merge_all(self,
              merge_order: (str,) = ('seqid', 'featuretype', 'strand', 'start'),
              merge_criteria: '[Callable]' = (mc.seqid, mc.overlap_end_inclusive, mc.strand, mc.feature_type),
              featuretypes_groups: 'Sequence[Set[str]]' = (None,),
              exclude_components: bool = False):
    """
    Merge all features in database according to criteria.
    Merged features will be assigned as children of the merged record.
    The resulting records are added to the database.

    :param merge_order: Ordered list of columns with which to group features before evaluating criteria
    :param merge_criteria: List of merge criteria callbacks. See merge().
    :param featuretypes_groups: iterable of sets of featuretypes to merge together
    :param exclude_components: True: child features will be discarded. False to keep them.
    :return: list of merge features
    """

    if not len(featuretypes_groups):
        # Can't be empty
        featuretypes_groups = (None,)

    result_features = []

    # Merge features per featuregroup
    for featuregroup in featuretypes_groups:
        for merged in merge(self, self.all_features(featuretype=featuregroup, order_by=merge_order),
                            merge_criteria=merge_criteria):
            # If feature is result of merge
            if merged.children:
                self._insert(merged, self.conn.cursor())
                if exclude_components:
                    # Remove child features from DB
                    self.delete(merged.children)
                else:
                    # Add child relations to DB
                    for child in merged.children:
                        self.add_relation(merged, child, 1, child_func=assign_child)

                result_features.append(merged)
            else:
                pass  # Do nothing, feature is already in DB

    return result_features


# -- Begin feature_merge specific code --

def get_args(sysargs):
    ignore_seqid = False
    ignore_strand = False
    ignore_featuretypes = False
    exact_only = False
    threshold = None
    exclude_components = False
    featuretypes_groups = []
    merge_strategy = "create_unique"
    merge_criteria = []
    # Parse arguments
    try:
        opts, args = getopt.gnu_getopt(sysargs, 'visecxf:m:t:')
        for opt, val in opts:
            if opt == '-v':
                from . import __version
                print(__version.__version__)
                exit(0)
            elif opt == '-i':
                ignore_strand = True
            elif opt == '-e':
                exclude_components = True
            elif opt == '-f':
                if val != "NONE": ignore_featuretypes = True
                if val != "ALL": featuretypes_groups.append(set(filter(None, val.split(','))))
            elif opt == '-x':
                exact_only = True
            elif opt == '-m':
                if val in merge_strategies:
                    merge_strategy = merge_strategies[val]
                else:
                    raise getopt.GetoptError("Invalid merge strategy", opt)
            elif opt == '-t':
                threshold = int(val)
            elif opt == '-s':
                ignore_seqid = True

    except getopt.GetoptError as err:
        # TODO raise exception rather than exit
        print("Argument error(", err.opt, "): ", err.msg, file=sys.stderr)
        args = []

    if len(args) < 1:
        print(usage, file=sys.stderr)
        # TODO raise exception rather than exit
        exit(1)

    # Remove any empty files as GFFutils gets angry
    paths = list(filter(os.path.getsize, args))
    if not len(paths):
        # TODO raise exception rather than exit
        exit(0)

    if not ignore_seqid:
        merge_criteria.append(mc.seqid)

    if exact_only:
        merge_criteria.append(mc.exact_coordinates_only)
    elif threshold is not None:
        merge_criteria.append(mc.overlap_any_threshold(threshold))
    else:
        merge_criteria.append(mc.overlap_any_inclusive)

    if not ignore_strand:
        merge_criteria.append(mc.strand)

    if ignore_featuretypes:
        if ignore_strand:
            merge_order = ('seqid', 'start', 'featuretype')
        else:
            merge_order = ('seqid', 'strand', 'start', 'featuretype')
    else:
        merge_criteria.append(mc.feature_type)
        if ignore_strand:
            merge_order = ('seqid', 'featuretype', 'start')
        else:
            merge_order = ('seqid', 'featuretype', 'strand', 'start')

    return paths, merge_strategy, merge_order, merge_criteria, featuretypes_groups, exclude_components


def load_data(paths: [str], merge_strategy: str = "create_unique") -> gffutils.FeatureDB:
    paths = list(filter(lambda f: os.path.getsize(f), paths))
    db = None
    while len(paths):
        try:
            db = gffutils.create_db(paths[0], ":memory:", merge_strategy=merge_strategy)
            break
        except ValueError as e:
            print("Error while parsing ", paths[0], e, file=sys.stderr)
            del paths[0]

    for path in paths[1:]:
        try:
            update(db, path, merge_strategy=merge_strategy)
        except ValueError as e:
            print("Error while parsing ", path, e, file=sys.stderr)

    if len(paths) == 0:
        raise ValueError("No valid input data")

    # Deal with autoincrements being behind by one
    for a in db._autoincrements:
        db._autoincrements[a] += 1

    return db
