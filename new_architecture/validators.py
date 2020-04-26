import new_architecture.utils as utils
import macrocycles.exceptions as exceptions


def validate_backbone(backbone):
    try:
        _, map_nums = zip(*utils.get_atom_map_nums(backbone))
        if not all(map(lambda x: x == 1, map_nums)):
            raise ValueError
    except ValueError:
        raise exceptions.InvalidMolecule(
            f'Backbone molecule is missing an atom map number or the atom map number is not equal to 1')

    return True
