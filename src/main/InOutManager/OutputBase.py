import abc


class OutputBaseClass(metaclass=abc.ABCMeta):

    """
    Base class for displaying or saving output files in correct format.
    """

    @abc.abstractmethod
    def display(self, result):
        '''
        Validate the content of the input file. Exit with error if the file is corrupted or contains strange character
        :param filepath: path to the input file
        :return: TRUE if the file can be loaded with the correct format
        '''

    @abc.abstractmethod
    def save(self, result, filepath):
        '''
        Load the input file for future uses
        :param filepath: path to the input file
        :return: a dictionary of inputs
        '''