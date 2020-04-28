import pytest
import new_architecture.ranges as ranges


@pytest.fixture
def static_range():
    yield ranges.Range(5, 10)


@pytest.mark.parametrize('test_min,test_max,static_range,expected', [(0, 2, '', False), (0, 5, '', False), (0, 6, '', True), (5, 7, '', True), (6, 9, '', True), (6, 15, '', True), (10, 15, '', False), (11, 15, '', False)], indirect=['static_range'])
def test_range(test_min, test_max, static_range, expected):
    dynamic_range = ranges.Range(test_min, test_max)
    assert((dynamic_range in static_range) == expected)
    assert((static_range in dynamic_range) == expected)


def test_whole_range(static_range):
    whole_range = ranges.WholeRange()
    assert(static_range in whole_range)
    assert(whole_range in static_range)


@pytest.mark.parametrize('test_min,test_max,static_range,expected', [(0, 2, '', False), (0, 5, '', False), (0, 6, '', True), (5, 7, '', True), (6, 9, '', True), (6, 15, '', True), (10, 15, '', False), (11, 15, '', False)], indirect=['static_range'])
def test_discrete_data_chunk(test_min, test_max, static_range, expected):
    int_avg = int((test_min + test_max) / 2)
    data_chunk = ranges.DiscreteDataChunk([test_min, int_avg, test_max])
    assert((data_chunk in static_range) == expected)
    assert(int_avg in data_chunk)
