from abc import ABC, abstractmethod

import project_io


class IProxy(ABC):

    @abstractmethod
    def __init__(self):
        pass

    def __getitem__(self, key):

        if not self.flag:  # pylint: disable=E0203
            self.data = self.hasher()  # pylint: disable=E1103
            self.flag = True

        return self.data[key]


class RegioSQMProxy(IProxy):

    def __init__(self):
        self.flag = False
        self.data = {}
        self.hasher = project_io.get_hashed_regiosqm_predictions


class pKaProxy(IProxy):

    def __init__(self):
        self.flag = False
        self.data = {}
        self.hasher = project_io.get_hashed_pka_predictions
