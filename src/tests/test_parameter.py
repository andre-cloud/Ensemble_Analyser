from src.parser_parameter import get_conf_parameters

import pytest
import mock


def test_get_conf_parameters():
    conf = mock.MagicMock()
    conf.mult = 1
    conf.weight_mass = 108.14
    conf.folder = "tests/assets"

    p = mock.MagicMock()
    p.number = 1
    p.calculator = "orca"
    p.freq_fact = 1
    p.freq = False

    log = mock.MagicMock()

    output = get_conf_parameters(
        conf=conf, number=1, p=p, time="1", temp=298.15, log=log
    )
    assert output is True

    p.number = 3
    p.freq = True
    output = get_conf_parameters(
        conf=conf, number=3, p=p, time="1", temp=298.15, log=log
    )
    assert output is True
    p.number = 4
    p.freq = True

    with pytest.raises(IOError) as exc_info:
        output = get_conf_parameters(
            conf=conf, number=4, p=p, time="1", temp=298.15, log=log
        )
    assert str(exc_info.value) == "No frequency in the output file"
