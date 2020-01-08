
from abc import ABC, abstractmethod
from copy import copy


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

    @abstractmethod
    def save(self):
        pass


class IDIterator:

    ALPHABET = 'a b c d e f g h i j k l m n o p q r s t u v w x y z'.split(' ')

    def __init__(self, id_io):

        self.id_io = id_io
        self.record = self.id_io.load()
        self.count = self.record['count']
        self.prefix = self.record['prefix']
        self.current = True

    def __del__(self):

        if not self.current:
            self.save()

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
        self.current = False
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

    def save(self):

        self.record['count'] = self.count
        self.record['prefix'] = self.prefix
        self.id_io.save(self.record)
        self.current = True


class IndexIterator:

    def __init__(self, index_io):

        self.index_io = index_io
        self.record = self.index_io.load()
        self.index = self.record['index']
        self.current = True

    def __del__(self):

        if not self.current:
            self.save()

    def get_next(self):

        index = self.index
        self.index += 1
        self.current = False
        return index

    def save(self):

        self.record['index'] = self.index
        self.index_io.save(self.record)
        self.current = True
