from argparse import ArgumentParser
import cpmg.repository as r
import cpmg.models as m

ROOT = 'all'


def get_print_options():
    model_strings = m.get_all_model_strings()
    model_strings.append(ROOT)
    return model_strings


def print_database(model_type):
    if model_type == ROOT:
        repo = r.repository_impl_from_string()
    elif model_type == m.Connection.STRING:
        repo = r.create_connection_repository()
    elif model_type == m.Backbone.STRING:
        repo = r.create_backbone_repository()
    elif model_type == m.Template.STRING:
        repo = r.create_template_repository()
    elif model_type == m.Sidechain.STRING:
        repo = r.create_sidechain_repository()
    elif model_type == m.Monomer.STRING:
        repo = r.create_monomer_repository()
    elif model_type == m.Peptide.STRING:
        repo = r.create_peptide_repository()
    elif model_type == m.TemplatePeptide.STRING:
        repo = r.create_template_peptide_repository()
    elif model_type == m.Macrocycle.STRING:
        repo = r.create_macrocycle_repository()
    elif model_type == m.RegioSQMPrediction.STRING:
        repo = r.create_regiosqm_repository()
    elif model_type == m.pKaPrediction.STRING:
        repo = r.create_pka_repository()
    elif model_type == m.PeptidePlan.STRING:
        repo = r.create_peptide_plan_repository()

    print(repo)


class QueryArgParser:
    def __init__(self):
        parser = ArgumentParser()
        parser.add_argument('-p', '--print', choices=get_print_options(), default=ROOT, nargs='?', const=ROOT,
                            help='Prints the number of records and structure if applicable in the repository for the specified type. Default type is all.')

        self.args = parser.parse_args()

    def execute(self):
        if self.args.print is not None:
            print_database(self.args.print)


if __name__ == "__main__":
    query_parser = QueryArgParser()
    query_parser.execute()
