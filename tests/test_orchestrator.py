import pytest
import mock
import cpmg.orchestrator as o


@pytest.fixture()
def single_process_sidechain_orchestrator():
    yield o.Orchestrator.from_strings(o.LEVEL_0, o.g.SidechainModifier.STRING)


@pytest.fixture(params=[5, 100])
def result_buffer(request):
    chunk_size = request.param
    mock_saver = mock.MagicMock(spec=o.h.SidechainDataHandler)
    # mock_saver.save.return_value = [str(uuid.uuid4()) for _ in range(len())]
    mock_saver.save.side_effect = lambda x: [f'test_id{i}' for i in range(len(x))]
    buffer = o.ResultBuffer(mock_saver, chunk_size=chunk_size)
    yield buffer, mock_saver, chunk_size or o.config.CAPACITY


@pytest.mark.parametrize('parallelism, expected', [
    (o.LEVEL_0, o.SingleProcessOrchestrator),
    (o.LEVEL_1, o.MultiProcessOrchestrator),
    (o.LEVEL_2, o.DistributedOrchestrator),
])
@pytest.mark.parametrize('operation', [
    (o.g.SidechainModifier.STRING),
    (o.g.MonomerGenerator.STRING),
    (o.g.PeptidePlanGenerator.STRING),
    (o.g.PeptideGenerator.STRING),
    (o.g.TemplatePeptideGenerator.STRING),
    (o.g.MacrocycleGenerator.STRING),
    (o.g.ConformerGenerator.STRING)
])
def test_orchestrator_from_strings(parallelism, operation, expected):
    orchestrator = o.Orchestrator.from_strings(parallelism, operation)

    assert isinstance(orchestrator.impl, expected)


def test_orchestrator_execute(single_process_sidechain_orchestrator):
    params = {'a': 1, 'b': 2}
    with mock.patch.object(o.SingleProcessOrchestrator, 'execute') as mock_execute:
        _ = single_process_sidechain_orchestrator.execute(**params)
        mock_execute.assert_called_once_with(**params)


def test_result_buffer():
    buffer = o.ResultBuffer(mock.MagicMock())

    assert buffer.chunk_size == o.config.CAPACITY


@pytest.mark.parametrize('chunk_size', [
    (-1),
    ('str'),
    (1.0)
])
def test_result_buffer_fail(chunk_size):
    with pytest.raises(o.InvalidChunkSize):
        _ = o.ResultBuffer(mock.MagicMock(), chunk_size=chunk_size)


def test_result_buffer_len(result_buffer):
    buffer, _, _ = result_buffer
    assert len(buffer) == 0


def test_result_buffer_add(result_buffer):
    buffer, _, chunk_size = result_buffer

    with mock.patch.object(buffer, 'flush') as mock_flush:
        for i in range(chunk_size - 1):
            buffer.add(i)
            assert len(buffer) == i + 1
            mock_flush.assert_not_called()

        buffer.add(0)
        mock_flush.assert_called_once()


def test_result_buffer_flush(result_buffer):
    buffer, mock_saver, _ = result_buffer

    buffer.add(0)
    buffer.flush()

    mock_saver.save.assert_called_once_with([0])
    assert buffer.buffer == []
    assert buffer.ids == mock_saver.save.side_effect([0])
