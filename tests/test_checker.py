import unittest
from gpcrml_lib import checker


class TestChecker(unittest.TestCase):
    """
    Tests for checker.py
    """
    def test_special_character(self):
        self.assertEqual(checker.special_character('iansdiuhad'), False)
        self.assertEqual(checker.special_character('654654648'), False)
        self.assertEqual(checker.special_character('asd#+-gfdg'), True)
        self.assertEqual(checker.special_character('6546!"§$654'), True)

    def test_len_is_four(self):
        self.assertEqual(checker.len_is_four('4567'), True)
        self.assertEqual(checker.len_is_four('ASDF'), True)
        self.assertEqual(checker.len_is_four('asd'), False)
        self.assertEqual(checker.len_is_four('6546!'), False)


# for running the tests directly in editor instead of command line: $ python3.6 -m unittest tests.test_checker
if __name__ == '__main__':
    unittest.main()
