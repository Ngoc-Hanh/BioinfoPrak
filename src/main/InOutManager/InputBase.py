import abc


class InputBaseClass(metaclass=abc.ABCMeta):

    """
    Base class for validating and loading input files.
    """

    @abc.abstractmethod
    def load(self, filepath):
        '''
        Load the input file for future uses
        :param filepath: path to the input file
        :return: a sequencesionary of inputs
        '''