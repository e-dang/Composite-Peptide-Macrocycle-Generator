
from copy import copy
from abc import ABC, abstractmethod
import project_io


class IRecordIterator(ABC):

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __del__(self):
        pass

    @abstractmethod
    def get_next(self):
        pass


class IDIterator:

    ALPHABET = 'a b c d e f g h i j k l m n o p q r s t u v w x y z'.split(' ')

    def __init__(self):

        self.io = project_io.JsonIDIO()
        self.record = list(self.io.load())[0]
        self.count = self.record['count']
        self.prefix = self.record['prefix']

    def __del__(self):

        self.record['count'] = self.count
        self.record['prefix'] = self.prefix
        self.io.save(self.record)

    def get_next(self):
        """
        Method for generating new IDs for parent_side_chains.

        Args:
            count (int): The count stored in the MongoDataBase, which is used for determining the next letter in the ID.
            prefix (str): The prefix to which the next letter in the _id will be appended to.

        Returns:
            str: The new ID.
        """

        try:
            code = self.prefix + self.ALPHABET[self.count]
        except IndexError:
            self.count = 0
            self.prefix = self.set_prefix(copy(self.prefix))
            code = self.prefix + self.ALPHABET[self.count]

        self.count += 1
        return code

    def set_prefix(self, prefix):
        """
        Recursive method for rotating the prefix letter once they reach 'z'. For example, a prefix 'zz' will turn into
        'aaa'.

        Args:
            prefix (str): The prefix to be rotated.

        Returns:
            str: The rotated prefix.
        """

        # initial prefix assignment
        if prefix == '':
            return 'a'

        # increment last letter of prefix and recursively wrap if necessary
        ending = ord(prefix[-1]) + 1
        if ending > ord('z'):
            prefix = self.set_prefix(prefix[:-1]) if len(prefix) > 1 else 'a'
            return prefix + 'a'

        # no recursive wrapping needed
        return prefix[:-1] + str(chr(ending))


class IndexIterator:

    def __init__(self):
        self.io = project_io.JsonIndexIO()
        self.record = list(self.io.load())[0]
        self.index = self.record['index']

    def __del__(self):

        self.record['index'] = self.index
        self.io.save(self.record)

    def get_next(self):

        index = self.index
        self.index += 1
        return index
