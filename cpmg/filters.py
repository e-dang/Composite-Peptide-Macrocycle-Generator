from functools import wraps

import cpmg.proxies as proxies
import cpmg.reactions as rxns
import cpmg.config as config


def regiosqm_filter(original_func):
    predictions = proxies.HashedRegioSQMProxy()
    types = (rxns.FriedelCrafts.TYPE, rxns.PictetSpangler.TYPE)

    @wraps(original_func)
    def regiosqm_filter_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            if original_result.type in types:
                if original_result.rxn_atom_idx in predictions[original_result.reacting_mol['kekule']]:
                    data.append(original_result)
            else:  # filter is not applicable to this type of reaction
                data.append(original_result)

        return data

    return regiosqm_filter_wrapper


def pka_filter(original_func):

    predictions = proxies.HashedpKaProxy()
    cutoff = config.PKA_CUTOFF
    types = (rxns.TsujiTrost.TYPE, )

    @wraps(original_func)
    def pka_filter_wrapper(*args, **kwargs):

        data = []
        for original_result in original_func(*args, **kwargs):
            if original_result.type in types:
                pkas = predictions[original_result.reacting_mol['kekule']]
                rxn_atom = original_result.rxn_atom_idx
                if pkas[str(rxn_atom)] < cutoff and pkas[str(rxn_atom)] > 0:
                    data.append(original_result)
            else:   # filter is not applicable to this type of reaction
                data.append(original_result)

        return data

    return pka_filter_wrapper
