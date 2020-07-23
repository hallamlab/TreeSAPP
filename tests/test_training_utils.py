import unittest


class ClassifierTester(unittest.TestCase):
    def test_augment_training_set(self):
        from treesapp.training_utils import augment_training_set
        import numpy as np
        test_array = np.array([0.33025841, 0.33025841, 0.07989843, 0.12523663, 0.13434813])
        result = augment_training_set(test_array)
        self.assertTrue(3, len(result))


if __name__ == '__main__':
    unittest.main()
