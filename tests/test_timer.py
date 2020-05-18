import importlib
from time import sleep

import pytest

import cpmg.timer as t


@pytest.fixture()
def reset_timer():
    importlib.reload(t)


@pytest.fixture()
def started_timer(reset_timer):
    timer = t.GlobalTimer.instance(time_allocation=2, buffer_time=1)
    timer.start()
    yield timer


def test_instance(reset_timer):
    timer1 = t.GlobalTimer.instance()
    timer2 = t.GlobalTimer.instance()

    assert timer1 is timer2
    assert not timer1.has_started()


def test_start(started_timer):
    assert started_timer.has_started()


def test_start_fail(started_timer):
    with pytest.raises(t.TimerAlreadyStarted):
        started_timer.start()


def test_elapsed_time(started_timer):
    sleep(0.5)

    time = started_timer.elapsed_time()

    assert abs(time - 0.5) < 0.01


def test_elapsed_time_fail(reset_timer):
    timer = t.GlobalTimer()

    with pytest.raises(t.TimerNotStarted):
        time = timer.elapsed_time()


def test_is_near_complete_true(started_timer):
    sleep(1)
    assert started_timer.is_near_complete()


def test_is_near_complete_false(started_timer):
    assert not started_timer.is_near_complete()
