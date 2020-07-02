import cpmg.config as config
import cpmg.hdf5 as hdf5
import cpmg.models as models
import cpmg.mongodb as mongo
import cpmg.utils as utils
from cpmg.ranges import Key, WholeRange


def add_completed_tag(record):
    record['completed'] = False
    return record


class AbstractRepository:
    TYPE = None

    def __init__(self, impl):
        self.impl = impl
        self.failed_instances = []

    def __repr__(self):
        return self.impl.__repr__()

    def load(self, key=Key(WholeRange())):
        for doc in self.impl.load(key):
            yield self.TYPE.from_dict(doc)

    def save(self, data):
        return self.impl.save(data)

    def get_num_records(self):
        return self.impl.get_num_records()

    def remove_records(self, key):
        return self.impl.remove_records(key)

    def deactivate_records(self, key):
        return self.impl.deactivate_records(key)

    def activate_records(self, key):
        return self.impl.activate_records(key)

    def load_inactive_records(self, key):
        for _id, data in self.impl.load_inactivate_records(key):
            yield self.TYPE.from_dict(data, _id=_id)

    def mark_complete(self, ids):
        self.impl.mark_complete(ids)

    def deactivate_completed(self):
        self.impl.deactivate_completed()

    def _check_type(self, data):
        for model in data:
            if not isinstance(model, self.TYPE):
                print(
                    f'Type error! Repository of type {self.TYPE} cannot save model of type {type(model)}. The instance has '
                    f'been skipped and saved in instance variable \'self.failed_instances\'.')
                self.failed_instances.append(model)
            else:
                yield add_completed_tag(model.to_dict())


class BackboneRepository(AbstractRepository):
    TYPE = models.Backbone
    STRING = TYPE.STRING

    def __init__(self, impl):
        super().__init__(impl.backbone_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class ConnectionRepository(AbstractRepository):
    TYPE = models.Connection
    STRING = TYPE.STRING

    def __init__(self, impl):
        super().__init__(impl.connection_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class TemplateRepository(AbstractRepository):
    TYPE = models.Template
    STRING = TYPE.STRING

    def __init__(self, impl):
        super().__init__(impl.template_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class SidechainRepository(AbstractRepository):
    TYPE = models.Sidechain
    STRING = TYPE.STRING

    def __init__(self, impl):
        super().__init__(impl.sidechain_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class MonomerRepository(AbstractRepository):
    TYPE = models.Monomer
    STRING = TYPE.STRING

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
    STRING = TYPE.STRING

    def __init__(self, impl):
        super().__init__(impl.peptide_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class TemplatePeptideRepository(AbstractRepository):
    TYPE = models.TemplatePeptide
    STRING = TYPE.STRING

    def __init__(self, impl):
        super().__init__(impl.template_peptide_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class MacrocycleRepository(AbstractRepository):
    TYPE = models.Macrocycle
    STRING = TYPE.STRING

    def __init__(self, impl):
        super().__init__(impl.macrocycle_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class ConformerRepository(AbstractRepository):
    TYPE = models.Conformer
    STRING = TYPE.STRING

    def __init__(self, impl):
        super().__init__(impl.conformer_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class ReactionRepository(AbstractRepository):
    TYPE = models.Reaction
    STRING = TYPE.STRING

    def __init__(self, impl):
        super().__init__(impl.reaction_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class RegioSQMRepository(AbstractRepository):
    TYPE = models.RegioSQMPrediction
    STRING = TYPE.STRING

    def __init__(self, impl):
        super().__init__(impl.regiosqm_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class pKaRepository(AbstractRepository):
    TYPE = models.pKaPrediction
    STRING = TYPE.STRING

    def __init__(self, impl):
        super().__init__(impl.pka_repo)

    def save(self, data):
        return super().save(self._check_type(data))


class PeptidePlanRepository(AbstractRepository):
    TYPE = models.PeptidePlan
    STRING = TYPE.STRING

    def __init__(self, impl):
        super().__init__(impl.peptide_plan_repo)

    def save(self, data):
        return super().save(self._check_type(data))

    def load(self, key):
        if key.peptide_length is None:
            return None

        for record in self.impl.load(key):
            yield (record['_id'], tuple(record['combination']))

    def _check_type(self, data):
        for record in data:
            if not isinstance(record, self.TYPE):
                print(
                    f'Type error! Repository of type {self.TYPE} cannot save model of type {type(record)}. The instance has '
                    f'been skipped and saved in instance variable \'self.failed_instances\'.')
                self.failed_instances.append(record)
            else:
                for combo in record.data():
                    yield add_completed_tag(combo)


class InactivesRepository(AbstractRepository):
    STRING = 'inactives'

    def __init__(self, impl):
        super().__init__(impl.inactives_repo)


def repository_impl_from_string(impl=None):
    impl = impl or config.DATA_FORMAT

    if impl == config.HDF5:
        return hdf5.HDF5Repository.instance()

    if impl == config.MONGO:
        return mongo.MongoRepository.instance()

    raise ValueError('Unrecognized repository implementation!')


def get_repository(repository):

    def repository_closure(impl=None):
        impl = impl or config.DATA_FORMAT
        if isinstance(impl, str):
            return repository(repository_impl_from_string(impl))

        return repository(impl)

    return repository_closure


create_backbone_repository = get_repository(BackboneRepository)
create_connection_repository = get_repository(ConnectionRepository)
create_template_repository = get_repository(TemplateRepository)
create_sidechain_repository = get_repository(SidechainRepository)
create_monomer_repository = get_repository(MonomerRepository)
create_peptide_repository = get_repository(PeptideRepository)
create_template_peptide_repository = get_repository(TemplatePeptideRepository)
create_macrocycle_repository = get_repository(MacrocycleRepository)
create_conformer_repository = get_repository(ConformerRepository)
create_reaction_repository = get_repository(ReactionRepository)
create_regiosqm_repository = get_repository(RegioSQMRepository)
create_pka_repository = get_repository(pKaRepository)
create_peptide_plan_repository = get_repository(PeptidePlanRepository)
create_inactives_repository = get_repository(InactivesRepository)


class CPMGRepository:
    STRING = 'all'

    def __init__(self, impl=None):
        self.connection_repo = create_connection_repository(impl)
        self.backbone_repo = create_backbone_repository(impl)
        self.template_repo = create_template_repository(impl)
        self.sidechain_repo = create_sidechain_repository(impl)
        self.monomer_repo = create_monomer_repository(impl)
        self.peptide_repo = create_peptide_repository(impl)
        self.template_peptide_repo = create_template_peptide_repository(impl)
        self.macrocycle_repo = create_macrocycle_repository(impl)
        self.conformer_repo = create_conformer_repository(impl)
        self.reaction_repo = create_reaction_repository(impl)
        self.regiosqm_repo = create_regiosqm_repository(impl)
        self.pka_repo = create_pka_repository(impl)
        self.peptide_plan_repo = create_peptide_plan_repository(impl)

    def __repr__(self):
        string = '/'
        for repo in self.__dict__.values():
            string += '\n\t' + repo.__repr__()

        return string

    def load(self, key):
        for repo in self.__dict__.values():
            try:
                for record in repo.load(key):
                    yield record
            except (KeyError, TypeError):
                continue


get_all_repository_strings = utils.get_module_strings(__name__)


def create_repository_from_string(string, impl=None):
    for _, member in utils.get_classmembers(__name__):
        try:
            if string == member.STRING:
                return member(repository_impl_from_string(impl))
        except AttributeError:
            pass

    raise ValueError(f'Unrecognized repository string: {string}')
