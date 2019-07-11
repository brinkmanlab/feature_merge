from unittest import TestCase
import gffutils


paths = (
    'test-data/615598_genome.gbk_1.gff3',
    'test-data/615598_genome.gbk_2.gff3',
    'test-data/915596_genome.gbk_1.gff3',
    'test-data/915596_genome.gbk_2.gff3',
    'test-data/1015597_genome.gbk_1.gff3',
    'test-data/1015597_genome.gbk_2.gff3',
)
num_features = 27

synthetic_path = 'test-data/synthetic.gff3'
num_synthetic_features = 18
num_synthetic_overlap = 13


class TestWithSynthDB(TestCase):
    def setUp(self) -> None:
        self.db = gffutils.create_db(synthetic_path, ":memory:", merge_strategy='create_unique')  # type: gffutils.FeatureDB
        self.assertEqual(num_synthetic_features, self.db.count_features_of_type())

    def _dump_db(self):
        return '\n'.join(map(str, self.db.all_features()))