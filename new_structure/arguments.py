from abc import ABC, abstractmethod
from itertools import product


class IArgumentProducer(ABC):

    @abstractmethod
    def __call__(self, *data):
        pass


class VanillaArgProducer(IArgumentProducer):

    def __call__(self, data):

        return data


class CartesianProductArgProducer(IArgumentProducer):

    def __init__(self, filter_func=lambda x: True):
        self.filter_func = filter_func

    def __call__(self, *data):
        return filter(self.filter_func, product(*data))


class PeptideGeneratorArgProducer(IArgumentProducer):

    def __call__(self, *data):

        pass
