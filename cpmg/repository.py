import cpmg.config as config
import cpmg.models as models
from cpmg.hdf5 import HDF5Initializer, HDF5Repository
from cpmg.ranges import WholeRange

HDF5 = 'hdf5'


class RepositoryInitializer:
    def __init__(self, impl):
        self.impl = impl

    def initialize(self):
        self.impl.initialize()


class AbstractRepository:
    TYPE = None
    CATEGORY = None

    def __init__(self, impl):
        self.impl = impl
        self.failed_instances = []

    def load(self, key=WholeRange()):
        for _id, data in self.impl.load(self.CATEGORY, key):
            yield self.TYPE.from_dict(data, _id=_id)

    def save(self, data):
        return self.impl.save(self.CATEGORY, data)

    def get_num_records(self):
        return self.impl.get_num_records(self.CATEGORY)

    def _check_type(self, data):
        for model in data:
            if not isinstance(model, self.TYPE):
                print(
                    f'Type error! Repository of type {self.TYPE} cannot save model of type {type(model)}. The instance has '
                    f'been skipped and saved in instance variable \'self.failed_instances\'.')
                self.failed_instances.append(model)
            else:
                yield model.to_dict()


class BackboneRepository(AbstractRepository):
    TYPE = models.Backbone
    CATEGORY = 'backbone'

    def __init__(self, impl):
        super().__init__(impl)

    def save(self, data):
        return super().save(self._check_type(data))


class ConnectionRepository(AbstractRepository):
    TYPE = models.Connection
    CATEGORY = 'connections'

    def __init__(self, impl):
        super().__init__(impl)

    def save(self, data):
        return super().save(self._check_type(data))


class TemplateRepository(AbstractRepository):
    TYPE = models.Template
    CATEGORY = 'templates'

    def __init__(self, impl):
        super().__init__(impl)

    def save(self, data):
        return super().save(self._check_type(data))


class SidechainRepository(AbstractRepository):
    TYPE = models.Sidechain
    CATEGORY = 'sidechains'

    def __init__(self, impl):
        super().__init__(impl)

    def save(self, data):
        return super().save(self._check_type(data))


class MonomerRepository(AbstractRepository):
    TYPE = models.Monomer
    CATEGORY = 'monomers'

    def __init__(self, impl):
        super().__init__(impl)

    def save(self, data):
        return super().save(self._attach_indices(self._check_type(data)))

    def _attach_indices(self, data):
        current_num_monomers = self.get_num_records()
        for i, monomer in enumerate(data):
            if monomer['index'] is None:
                monomer['index'] = current_num_monomers + i
                yield monomer
            else:
                print(
                    f'Warning - A monomer with an index != None indicates that the monomer has already been saved to '
                    f'the repository! Saving monomer in instance variable \'self.failed_instances\' kekule: '
                    f'{monomer.kekule} index: {monomer.index}')
                self.failed_instances.append(monomer)


class PeptideRepository(AbstractRepository):
    TYPE = models.Peptide
    CATEGORY = 'peptides'

    def __init__(self, impl):
        super().__init__(impl)

    def save(self, data):
        return super().save(self._check_type(data))


class TemplatePeptideRepository(AbstractRepository):
    TYPE = models.TemplatePeptide
    CATEGORY = 'template_peptides'

    def __init__(self, impl):
        super().__init__(impl)

    def save(self, data):
        return super().save(self._check_type(data))


class RegioSQMRepository(AbstractRepository):
    TYPE = models.RegioSQMPrediction
    CATEGORY = 'regiosqm'

    def __init__(self, impl):
        super().__init__(impl)

    def save(self, data):
        return super().save(self._check_type(data))


def repository_impl_from_string(impl):
    if impl == HDF5:
        return HDF5Repository()
    else:
        raise ValueError('Unrecognized repository implementation!')


def create_repository_initializer(impl=None):
    impl = impl or config.DATA_FORMAT
    if impl == HDF5:
        return RepositoryInitializer(HDF5Initializer())
    else:
        raise ValueError('Unrecognized repository implementation!')


def get_repository(repository):
    def repository_closure(impl=None):
        impl = impl or config.DATA_FORMAT
        return repository(repository_impl_from_string(impl))

    return repository_closure


create_backbone_repository = get_repository(BackboneRepository)
create_connection_repository = get_repository(ConnectionRepository)
create_template_repository = get_repository(TemplateRepository)
create_sidechain_repository = get_repository(SidechainRepository)
create_monomer_repository = get_repository(MonomerRepository)
create_peptide_repository = get_repository(PeptideRepository)
create_template_peptide_repository = get_repository(TemplatePeptideRepository)
create_regiosqm_repository = get_repository(RegioSQMRepository)

# class SaveRequest:
#     def __init__(self, data, data_type, peptide_length=None, job_num=None):
#         self.data = data
#         self.data_type = data_type
#         self.peptide_length = peptide_length
#         self.job_num = job_num

#     @classmethod
#     def create_request(cls, data, meta_data):
#         return cls(data, meta_data.data_type, meta_data.peptide_length, meta_data.job_num)


# class LoadRequest:

#     def __init__(self, data_type, peptide_length=None, data_chunk=WholeRange()):
#         self.data_type = data_type
#         self.peptide_length = peptide_length
#         self.data_chunk = data_chunk
#         self.next_request = None

#     def add_request(self, request):
#         if self.next_request is None:
#             self.next_request = request
#         else:
#             self.next_request.add_request(request)

#     def get_next_request(self, result):
#         try:
#             self.next_request.update_data_chunk(result)
#         except AttributeError:
#             pass

#         return self.next_request

#     def update_data_chunk(self, data_chunk):
#         self.data_chunk = data_chunk

#     def format_result(self, func):
#         self.result = func(self.result)
