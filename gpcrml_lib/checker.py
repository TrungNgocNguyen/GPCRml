def len_is_four(string):
    """
    Check if string length is 4

    :param string: String
    :return: True/False
    """
    if len(string) == 4:
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
