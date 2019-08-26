def len_is_x(string, x):
    """
    Check if string length is 4

    :param string: String
    :return: True/False
    """
    if len(string) == x:
        return True
    else:
        return False


special_characters = ['~', ':', "'", '+', '[', '\\', '@', '^', '{', '%', '(', '-', '"', '*', '|', ',',
                      '&', '<', '`', '}', '.', '_', '=', ']', '!', '>', ';', '?', '#', '$', ')', '/']


def special_character(string, special_characters = special_characters):
    """
    Check if string contains special characters
    :param string: String
    :return: True/False
    """

    if any(letter in string for letter in special_characters):
        return True
    else:
        return False


def capital_letter_or_digit(string):
    """
    Check if string contains only capital letters or digits
    :param string: e.g. '3V2Y'
    :return: True/False
    """
    for letter in string:
        if letter.isupper() or letter.isdigit():
            continue
        else:
            return False
    return True
