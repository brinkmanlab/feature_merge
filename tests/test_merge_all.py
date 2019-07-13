from feature_merge import merge_all
from . import TestWithSynthDB, num_synthetic_features, num_synthetic_overlap


class TestMerge_all(TestWithSynthDB):
    def test_defaults(self):
        merged = merge_all(self.db)
        dump = self._dump_db()
        self.assertEqual(num_synthetic_features + 1, self.db.count_features_of_type(), dump)
        self.assertEqual(1, len(merged), dump)
        self.assertEqual(num_synthetic_overlap, len(merged[0].children), dump)

    def test_empty(self):
        self.db.delete(self.db.all_features())
        dump = self._dump_db()
        self.assertEqual(0, self.db.count_features_of_type(), dump)
        merged = merge_all(self.db)
        dump = self._dump_db()
        self.assertEqual(0, self.db.count_features_of_type(), dump)
        self.assertEqual(0, len(merged), dump)

    def test_one(self):
        features = self.db.all_features()
        next(features)
        self.db.delete(features)
        dump = self._dump_db()
        self.assertEqual(1, self.db.count_features_of_type(), dump)
        merged = merge_all(self.db)
        dump = self._dump_db()
        self.assertEqual(1, self.db.count_features_of_type(), dump)
        self.assertEqual(0, len(merged), dump)

    def test_no_overlap(self):
        self.db.delete(f for f in self.db.all_features() if f.id not in ('basic1', 'no_overlap1'))
        dump = self._dump_db()
        self.assertEqual(2, self.db.count_features_of_type(), dump)
        merged = merge_all(self.db)
        dump = self._dump_db()
        self.assertEqual(2, self.db.count_features_of_type(), dump)
        self.assertEqual(0, len(merged), dump)

    def test_merge_groups(self):
        merged = merge_all(self.db, featuretypes_groups=({'sequence_feature', 'misc_feature'},))
        dump = self._dump_db()
        self.assertEqual(num_synthetic_features + 1, self.db.count_features_of_type(), dump)
        self.assertEqual(1, len(merged), dump)
        self.assertEqual(num_synthetic_overlap, len(merged[0].children), dump)

    def test_exclude_components(self):
        merged = merge_all(self.db, exclude_components=True)
        dump = self._dump_db()
        self.assertEqual(6, self.db.count_features_of_type(), dump)
        self.assertEqual(1, len(merged), dump)
        self.assertEqual(num_synthetic_overlap, len(merged[0].children), dump)
