import pytest

import cpmg.parallelizers as parallelizers


def test_get_all_parallelizer_strings():
    parallelizer_strings = set(parallelizers.get_all_parallelizer_strings())

    assert parallelizer_strings == {parallelizers.SingleProcess.STRING,
                                    parallelizers.Multiprocess.STRING}


@pytest.mark.parametrize('parallelizer', [
    (parallelizers.SingleProcess),
    (parallelizers.Multiprocess)
])
def test_create_parallelizer_from_string(parallelizer):
    produced_parallelizer = parallelizers.create_parallelizer_from_string(parallelizer.STRING)

    assert isinstance(produced_parallelizer, parallelizer)
