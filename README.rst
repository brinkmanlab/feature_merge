.. image:: https://zenodo.org/badge/166121742.svg
    :target: https://zenodo.org/badge/latestdoi/166121742

=============
feature_merge
=============
Merge features in GFF files

Installation::

    pip install feature-merge

or::

    conda install feature_merge

or::

    git clone https://github.com/brinkmanlab/feature_merge.git
    cd feature_merge
    ./setup.py install

Usage::

    feature_merge [-i] [-e] [-x] [-s] [-v] [-m merge|append|error|skip|replace] [-f type[,type..]].. <input1> [<input_n>..]
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
        
See CONTRIBUTING.rst_ for information on contributing to this repo.

.. _CONTRIBUTING.rst: CONTRIBUTING.rst
