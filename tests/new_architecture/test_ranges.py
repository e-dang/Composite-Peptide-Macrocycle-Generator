import pytest

import new_architecture.ranges as ranges


@pytest.fixture
def static_range():
    yield ranges.Range(5, 10)


@pytest.mark.parametrize('test_min,test_max,static_range,expected', [(0, 2, '', False), (0, 5, '', False), (0, 6, '', True), (5, 7, '', True), (6, 9, '', True), (6, 15, '', True), (10, 15, '', False), (11, 15, '', False)], indirect=['static_range'])
def test_range(test_min, test_max, static_range, expected):
    dynamic_range = ranges.Range(test_min, test_max)
    assert((dynamic_range in static_range) == expected)
