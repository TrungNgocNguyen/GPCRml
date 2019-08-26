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

    def test_len_is_x(self):
        self.assertEqual(checker.len_is_x('4567', 4), True)
        self.assertEqual(checker.len_is_x('ASDF1', 5), True)
        self.assertEqual(checker.len_is_x('asd', 2), False)
        self.assertEqual(checker.len_is_x('6546!', 4.5), False)

    def test_capital_letter_or_digit(self):
        self.assertEqual(checker.capital_letter_or_digit('45AS'), True)
        self.assertEqual(checker.capital_letter_or_digit('4As3'), False)
        self.assertEqual(checker.capital_letter_or_digit('ASD!'), False)
        self.assertEqual(checker.capital_letter_or_digit('6548'), True)


# for running the tests directly in editor instead of command line: $ python3.6 -m unittest tests.test_checker
if __name__ == '__main__':
    unittest.main()
