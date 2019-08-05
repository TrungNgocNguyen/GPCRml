import os


def dir_exist(directory):
    """
    Check if directory exist, if not make directory
    :param directory: Directory to check
    :return: True/False
    """
    if not os.path.exists(directory):
        os.mkdir(directory)


def pickle_exist(directory, pickle_name):
    """
    Check if pickle file exist. If yes return pickle, if no return empty dict.
    :param directory: Directory of pdb_code_dict
    :param pickle_name: Name of pickle in directory (without .p)
    :return: Pickled dict if exist, if not return empty dict
    """
    if os.path.exists('{}/{}.p'.format(directory, pickle_name)):
        loaded_pickle = pickle.load(open('{}/{}.p'.format(directory, pickle_name), 'rb'))
        return loaded_pickle
    else:
        return {}


special_characters = ['~', ':', "'", '+', '[', '\\', '@', '^', '{', '%', '(', '-', '"', '*', '|', ',', '&', '<', '`', '}', '.', '_', '=', ']', '!', '>', ';', '?', '#', '$', ')', '/']


def special_character_checker(string, special_characters = special_characters):
    """
    Check if string contains special characters
    :param string: String
    :return: True/False
    """

    if any(letter in string for letter in special_characters):
        return True
    else:
        return False
