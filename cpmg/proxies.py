import cpmg.repository as repo


def hash_predictions(predictions):
    hashed_predictions = {}
    for prediction in predictions:
        hashed_predictions[prediction.reacting_mol] = prediction.predictions

    return hashed_predictions


class AbstractProxy():

    def __getitem__(self, key):

        if not self.flag:  # pylint: disable=E0203
            self.data = self.load_data()  # pylint: disable=E1103
            self.flag = True

        return self.data[key]

    def load_data(self):
        return hash_predictions(self.loader.load())  # pylint: disable=E1103


class HashedRegioSQMProxy(AbstractProxy):

    def __init__(self):
        self.loader = repo.create_regiosqm_repository()
        self.flag = False
        self.data = {}


class HashedpKaProxy(AbstractProxy):

    def __init__(self):
        self.loader = repo.create_pka_repository()
        self.flag = False
        self.data = {}
