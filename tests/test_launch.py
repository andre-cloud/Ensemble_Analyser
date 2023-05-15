from src.launch import launch
import mock
from ase.calculators.calculator import CalculationFailed
import pytest

@mock.patch("builtins.open", mock.mock_open())
def test_launch():
    conf = mock.MagicMock()
    conf.get_ase_atoms = mock.Mock(side_effect=CalculationFailed("test"))
    log = mock.MagicMock()
    log.error = mock.MagicMock()
    log.critical = mock.MagicMock()
    protocol = mock.MagicMock()
    protocol.get_calculator = mock.Mock(return_value=(mock.MagicMock(), "test label"))

    with pytest.raises(RuntimeError):
        launch(
        1, conf, protocol, mock.MagicMock(), log, 0, mock.MagicMock()
        )
    log.critical.assert_called_once()
    log.error.assert_called_once()

