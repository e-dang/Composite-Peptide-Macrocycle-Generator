import cpmg.config as config
import cpmg.models as models
import cpmg.hdf5 as hdf5
from cpmg.ranges import WholeRange, Key

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

    def __repr__(self):
        return self.impl.__repr__()

    def load(self, key=Key(WholeRange())):
        for _id, data in self.impl.load(key):
            yield self.TYPE.from_dict(data, _id=_id)

    def save(self, data):
        return self.impl.save(data)

    def get_num_records(self):
        return self.impl.get_num_records()

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
        super().__init__(impl.backbone_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class ConnectionRepository(AbstractRepository):
    TYPE = models.Connection
    CATEGORY = 'connections'

    def __init__(self, impl):
        super().__init__(impl.connection_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class TemplateRepository(AbstractRepository):
    TYPE = models.Template
    CATEGORY = 'templates'

    def __init__(self, impl):
        super().__init__(impl.template_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class SidechainRepository(AbstractRepository):
    TYPE = models.Sidechain
    CATEGORY = 'sidechains'

    def __init__(self, impl):
        super().__init__(impl.sidechain_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class MonomerRepository(AbstractRepository):
    TYPE = models.Monomer
    CATEGORY = 'monomers'

    def __init__(self, impl):
        super().__init__(impl.monomer_repo)

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
        super().__init__(impl.peptide_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class TemplatePeptideRepository(AbstractRepository):
    TYPE = models.TemplatePeptide
    CATEGORY = 'template_peptides'

    def __init__(self, impl):
        super().__init__(impl.template_peptide_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class MacrocycleRepository(AbstractRepository):
    TYPE = models.Macrocycle
    CATEGORY = 'macrocycles'

    def __init__(self, impl):
        super().__init__(impl.macrocycle_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class ReactionRepository(AbstractRepository):
    TYPE = models.Reaction
    CATEGORY = 'reactions'

    def __init__(self, impl):
        super().__init__(impl.reaction_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class RegioSQMRepository(AbstractRepository):
    TYPE = models.RegioSQMPrediction
    CATEGORY = 'regiosqm'

    def __init__(self, impl):
        super().__init__(impl.regiosqm_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class pKaRepository(AbstractRepository):
    TYPE = models.pKaPrediction
    CATEGORY = 'pka'

    def __init__(self, impl):
        super().__init__(impl.pka_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class PeptidePlanRepository(AbstractRepository):
    TYPE = models.PeptidePlan
    CATEGORY = 'peptide_plan'

    def __init__(self, impl):
        super().__init__(impl.peptide_plan_repo)

    def save(self, data):
        return super().save(self._check_type(data))

    def load(self, key):
        if key.peptide_length is None:
            raise ValueError('Must specify a peptide length when loading a peptide plan!')

        return self.TYPE.from_array_tuple(key.peptide_length, self.impl.load(key))

    def _check_type(self, data):
        if not isinstance(data, self.TYPE):
            print(
                f'Type error! Repository of type {self.TYPE} cannot save model of type {type(data)}. The instance has '
                f'been skipped and saved in instance variable \'self.failed_instances\'.')
            self.failed_instances.append(data)
        else:
            return data.data()


def repository_impl_from_string(impl=None):
    impl = impl or config.DATA_FORMAT

    if impl == HDF5:
        return hdf5.HDF5Repository.instance()

    raise ValueError('Unrecognized repository implementation!')


# def create_repository_initializer(impl=None):
#     impl = impl or config.DATA_FORMAT
#     if impl == HDF5:
#         return RepositoryInitializer(HDF5Initializer())
#     else:
#         raise ValueError('Unrecognized repository implementation!')


def get_repository(repository):

    def repository_closure(impl=None):
        return repository(repository_impl_from_string(impl or config.DATA_FORMAT))

    return repository_closure


create_backbone_repository = get_repository(BackboneRepository)
create_connection_repository = get_repository(ConnectionRepository)
create_template_repository = get_repository(TemplateRepository)
create_sidechain_repository = get_repository(SidechainRepository)
create_monomer_repository = get_repository(MonomerRepository)
create_peptide_repository = get_repository(PeptideRepository)
create_template_peptide_repository = get_repository(TemplatePeptideRepository)
create_macrocycle_repository = get_repository(MacrocycleRepository)
create_reaction_repository = get_repository(ReactionRepository)
create_regiosqm_repository = get_repository(RegioSQMRepository)
create_pka_repository = get_repository(pKaRepository)
create_peptide_plan_repository = get_repository(PeptidePlanRepository)
