import unittest
import common_objects
from simply_python import data_structures as ds

class MutationTester(unittest.TestCase):

    def test_mutation_maker(self):
        seq= ds.StringBuffer('A'*20)
        seq2= common_objects.mutate_sequence(seq, percent_id=0.5)
        self.assertNotEqual(seq, seq2)
        self.assertNotEqual(ds.StringBuffer('A'*20), seq2)
        self.assertEqual(ds.StringBuffer('A'*20), seq)

if __name__ == '__main__':
    unittest.main()
