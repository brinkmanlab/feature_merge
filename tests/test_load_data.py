from unittest import TestCase, expectedFailure

from feature_merge import load_data, gffutils
from . import paths, num_features

class TestLoad_data(TestCase):
    empty_path = 'test-data/empty_test.gff'
    header_only_path = 'test-data/empty_test_header.gff'

    def test_load_data(self):
        db = load_data(paths) # type: gffutils.FeatureDB
        self.assertEqual(num_features, db.count_features_of_type())

    @expectedFailure
    def test_load_empty(self):
        db = load_data((self.empty_path,))
        self.assertEqual(0, db.count_features_of_type())

    @expectedFailure
    def test_load_header_only(self):
        db = load_data((self.header_only_path,))
        self.assertEqual(0, db.count_features_of_type())