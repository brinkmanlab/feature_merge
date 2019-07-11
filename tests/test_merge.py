from copy import deepcopy

from feature_merge import merge, mc
from . import TestWithSynthDB, num_synthetic_features, num_synthetic_overlap


class TestMerge(TestWithSynthDB):
    def _merge_and_compare(self, ids, expected_merge_children, **kwargs):
        features = [self.db[id] for id in ids]
        features_backup = deepcopy(features)
        merged = list(merge(self.db, features, **kwargs))
        dump = self._dump_db()
        m_dump = '\n'.join(map(str, merged))

        # Ensure db wasn't modified
        self.assertEqual(num_synthetic_features, self.db.count_features_of_type(), dump)
        for expected, actual in zip(features_backup, ids):
            self.assertEqual(expected, self.db[actual], '\n'.join(map(str, [expected, self.db[actual]])) + '\n')

        # Ensure input Features not modified
        for expected, actual in zip(features_backup, features):
            self.assertEqual(expected, actual, '\n'.join(map(str, [expected, actual])) + '\n')

        # Ensure output matches expected
        self.assertEqual(len(expected_merge_children), len(merged), m_dump)
        for expected, actual in zip(expected_merge_children, merged):
            self.assertEqual(len(expected), len(actual.children), '\n'.join(map(str, [expected, actual])) + '\n')
            self.assertTrue(all(child.id in expected for child in actual.children))

        return merged

    def test_defaults(self):
        merged = list(merge(self.db, self.db.all_features(order_by=('seqid', 'featuretype', 'strand', 'start'))))
        dump = '\n'.join(map(str, merged))
        self.assertEqual(num_synthetic_features, self.db.count_features_of_type(), dump)
        self.assertEqual(6, len(merged), dump)
        merged_count = 0
        for f in merged:
            if f.children:
                self.assertEqual(num_synthetic_overlap, len(f.children), dump)
                merged_count += 1

        self.assertEqual(1, merged_count)

    def test_none(self):
        merged = list(merge(self.db, []))
        dump = self._dump_db()
        self.assertEqual(num_synthetic_features, self.db.count_features_of_type(), dump)
        self.assertEqual(0, len(merged))

    def test_one(self):
        self._merge_and_compare(['basic1'], [[]])

    def test_no_overlap(self):
        self._merge_and_compare(['basic1', 'no_overlap1'], [[],[]])

    def test_perfect_overlap(self):
        self._merge_and_compare(['basic1', 'perfect_overlap1'], [['basic1', 'perfect_overlap1']])

    def test_overlap(self):
        self._merge_and_compare(['basic1', 'adjacent1'], [['basic1', 'adjacent1']])

    def test_any_overlap(self):
        self._merge_and_compare(['basic1', 'overlap_start1'], [['basic1', 'overlap_start1']])
        self._merge_and_compare(['overlap_start1', 'basic1'], [['basic1', 'overlap_start1']])
        self._merge_and_compare(['adjacent1', 'basic1'], [['basic1', 'adjacent1']])
        self._merge_and_compare(['adjacent4', 'basic1'], [['basic1', 'adjacent4']])

    def test_adjacent_len1(self):
        self._merge_and_compare(['basic1', 'adjacent2'], [['basic1', 'adjacent2']])
        self._merge_and_compare(['adjacent2', 'basic1'], [['basic1', 'adjacent2']])
        self._merge_and_compare(['basic1', 'adjacent3'], [['basic1', 'adjacent3']])
        self._merge_and_compare(['adjacent3', 'basic1'], [['basic1', 'adjacent3']])

    def test_end_overlap(self):
        self._merge_and_compare(['basic1', 'overlap_end1'], [['basic1', 'overlap_end1']])

    def test_ignore_strand(self):
        self._merge_and_compare(['basic1', 'strand_plus1'], [['basic1', 'strand_plus1']], ignore_strand=True)
        self._merge_and_compare(['basic1', 'strand_minus1'], [['basic1', 'strand_minus1']], ignore_strand=True)
        self._merge_and_compare(['basic1', 'strand_plus1', 'strand_minus1'], [['basic1', 'strand_plus1', 'strand_minus1']], ignore_strand=True)
        self._merge_and_compare(['basic1', 'strand_plus1'], [['basic1', 'strand_plus1']], merge_criteria=(mc.any_overlap_inclusive, mc.feature_type))
        self._merge_and_compare(['basic1', 'strand_minus1'], [['basic1', 'strand_minus1']], merge_criteria=(mc.any_overlap_inclusive, mc.feature_type))
        self._merge_and_compare(['basic1', 'strand_plus1', 'strand_minus1'], [['basic1', 'strand_plus1', 'strand_minus1']], merge_criteria=(mc.any_overlap_inclusive, mc.feature_type))

    # TODO test various feature orders

    # TODO test various merge criteria

    # TODO test multiline, this doesn't really work until gffutils better supports multiline features