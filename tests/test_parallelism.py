import importlib

import mock
import pytest

import cpmg.parallelism as p


@pytest.fixture()
def reset_parallelism():
    importlib.reload(p)


@pytest.fixture()
def mock_mpi(reset_parallelism):
    with mock.patch.object(p, 'MPI') as mock_mpi_obj:
        yield mock_mpi_obj


@pytest.fixture()
def mock_atexit():
    with mock.patch.object(p, 'atexit') as mock_atexit_obj:
        yield mock_atexit_obj


def test_get_parallelism_strings():
    assert set(p.get_parallelism_strings()) == {p.LEVEL_0,
                                                p.LEVEL_1,
                                                p.LEVEL_2}


def test_set_get_level(reset_parallelism):
    assert p.Parallelism.get_level() == None

    p.Parallelism.set_level(p.LEVEL_0)

    assert p.Parallelism.get_level() == p.LEVEL_0


def test_set_get_level_disributed(mock_mpi, mock_atexit):
    p.Parallelism.set_level(p.LEVEL_2)

    mock_atexit.register.assert_called_once_with(mock_mpi.Finalize)
    mock_mpi.Init_thread.assert_called_once()

    assert p.Parallelism.get_level() == p.LEVEL_2


def test_set_level_fail(reset_parallelism):
    p.Parallelism.set_level(p.LEVEL_0)

    with pytest.raises(p.ParallelismAlreadySet):
        p.Parallelism.set_level(p.LEVEL_1)


def test_is_distributed(mock_mpi, mock_atexit):
    p.Parallelism.set_level(p.LEVEL_2)

    assert p.Parallelism.is_distributed()


def test_is_distributed_fail(mock_mpi, mock_atexit):
    p.Parallelism.set_level(p.LEVEL_1)

    assert not p.Parallelism.is_distributed()
