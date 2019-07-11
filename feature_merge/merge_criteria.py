"""
Merge criteria used by FeatureDB merge() and merge_all()
"""
from gffutils import Feature


def exact_coordinates_only(acc: Feature, cur: Feature, components: [Feature]):
    return cur.start == acc.start and cur.stop == acc.end and cur.seqid == acc.seqid


def any_overlap_inclusive(acc: Feature, cur: Feature, components: [Feature]):
    return cur.seqid == acc.seqid and (
            acc.start <= cur.start <= acc.end + 1 or
            acc.start <= cur.end + 1 <= acc.end + 1
    )


def strand(acc: Feature, cur: Feature, components: [Feature]):
    return acc.strand == cur.strand


def feature_type(acc: Feature, cur: Feature, components: [Feature]):
    return acc.featuretype == cur.featuretype