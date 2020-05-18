import cpmg.repository as repo


def hash_predictions(predictions):
    hashed_predictions = {}
    for prediction in predictions:
        hashed_predictions[prediction.reacting_mol] = prediction.predictions

    return hashed_predictions


class AbstractProxy():

    def __getitem__(self, key):

        if not self.flag:  # pylint: disable=E0203
            self.load_data()
            self.flag = True

        return self.data[key]  # pylint: disable=no-member

    def load_data(self):
        pass


class HashedRegioSQMProxy(AbstractProxy):

    def __init__(self):
        self.flag = False
        self.data = {}

    def load_data(self):
        self.data = hash_predictions(repo.create_regiosqm_repository().load())


class HashedpKaProxy(AbstractProxy):

    def __init__(self):
        self.flag = False
        self.data = {}

    def load_data(self):
        self.data = hash_predictions(repo.create_pka_repository().load())
