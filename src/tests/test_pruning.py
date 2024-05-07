import pytest
import ase
import json
import numpy as np
from mock import Mock, patch


from src.pruning import (
    cut_over_thr_max,
    rmsd,
    check_ensemble,
    calculate_rel_energies,
    check,
    refactor_dict,
)


class MockConformer:
    def __init__(self, number, energy):
        self.active = True
        self.number = number
        self.get_energy = energy
        self._last_energy = {}
        self.rotatory = 1
        self.moment = 1

    def get_ase_atoms(self):
        return ase.Atoms(["H", "H"], np.array([[0, 0, 0], [0, 1, 0]]))

    def write_xyz(self):
        return "a"


def test_cut_over_thr_max():
    conf1 = MockConformer(1, 10)
    conf2 = MockConformer(2, 5)
    conf3 = MockConformer(3, 15)
    confs = [conf1, conf2, conf3]
    thrGMAX = 8.0
    log = Mock()

    cut_over_thr_max(confs, thrGMAX, log)

    assert conf1.active is True
    assert conf2.active is True
    assert conf3.active is False


def test_check_ensemble():
    conf1 = MockConformer(1, 10)
    conf2 = MockConformer(2, 20)
    conf3 = MockConformer(3, 30)

    confs = [conf1, conf2, conf3]
    protocol = Mock(graph=False, thrGMAX=25, thrG=1, thrB=2)
    log = Mock()

    with patch("src.pruning.cut_over_thr_max") as mock_cut_over_thr_max, patch(
        "src.pruning.rmsd"
    ) as mock_rmsd:
        result = check_ensemble(confs, protocol, log)

    assert result == confs  # The function should return the original `confs` list
    # The `protocol.graph` attribute should not be changed
    assert protocol.graph is False

    mock_cut_over_thr_max.assert_called_once_with(confs, 25, log)
    assert mock_rmsd.call_count == 3
    mock_rmsd.assert_any_call(conf1.get_ase_atoms(), conf2.get_ase_atoms())
    mock_rmsd.assert_any_call(conf1.get_ase_atoms(), conf3.get_ase_atoms())
    mock_rmsd.assert_any_call(conf2.get_ase_atoms(), conf3.get_ase_atoms())

    confs = [conf1, conf2, conf3]
    protocol.graph = True
    protocol.number = 1

    with patch("src.pruning.cut_over_thr_max") as mock_cut_over_thr_max, patch(
        "src.pruning.rmsd"
    ) as mock_rmsd:
        result = check_ensemble(confs, protocol, log)
    assert result == confs

    conf2.active = False
    with patch("src.pruning.cut_over_thr_max") as mock_cut_over_thr_max, patch(
        "src.pruning.rmsd"
    ) as mock_rmsd:
        result = check_ensemble(confs, protocol, log)

    confs_tmp = confs[:]
    confs.pop(1)
    assert result == confs

    confs = confs_tmp[:]


def test_calculate_rel_energies():
    c = [MockConformer(1, 10.0), MockConformer(2, 20.0), MockConformer(3, 15.0)]
    T = 298.15

    calculate_rel_energies(c, T)

    assert np.isclose(c[0]._last_energy["Erel"], 0.0)
    assert np.isclose(c[0]._last_energy["Pop"], 99.97846114451659)

    assert np.isclose(c[1]._last_energy["Erel"], 10.0)
    assert np.isclose(c[1]._last_energy["Pop"], 4.638224150092233e-06)

    assert np.isclose(c[2]._last_energy["Erel"], 5.0)
    assert np.isclose(c[2]._last_energy["Pop"], 0.02153421725927262)


def test_rmsd():
    # check_positions = np.array([[0, 0, 0], [1, 1, 1]])
    # ref_positions = np.array([[0, 0, 0], [1, 0, 0]])

    check_conformer = MockConformer(1, 10)
    ref_conformer = MockConformer(1, 10)

    # expected_rmsd = np.sqrt(1 / 3) * np.linalg.norm(check_positions - ref_positions)

    assert rmsd(check_conformer.get_ase_atoms(), ref_conformer.get_ase_atoms()) == 0


def test_check():
    check_conformer = MockConformer(2, 10)
    check_conformer.rotatory = 5
    check_conformer.moment = 2
    ref_conformer = MockConformer(1, 5)
    ref_conformer.rotatory = 3
    ref_conformer.moment = 1
    protocol = Mock(graph=False, thrGMAX=25, thrG=1, thrB=2)

    controller = {}

    assert check(check_conformer, ref_conformer, protocol, controller) == False
    assert check_conformer.active == True
    assert controller == {}

    ref_conformer.active = False

    assert check(check_conformer, ref_conformer, protocol, controller) == False
    assert check_conformer.active == True
    assert controller == {}

    check_conformer = MockConformer(3, 3)
    check_conformer.rotatory = 5
    check_conformer.moment = 2

    ref_conformer = MockConformer(1, 5)
    ref_conformer.rotatory = 3
    ref_conformer.moment = 1

    protocol = Mock(graph=False, thrGMAX=25, thrG=3, thrB=2.1)
    controller = {}

    assert check(check_conformer, ref_conformer, protocol, controller) == True
    assert check_conformer.active == False
    assert check_conformer.diactivated_by == 1
    assert controller[0]["∆E [kcal/mol]"] == -2
    assert controller[0]["∆B [e-3 cm-1]"] == 2000
    assert controller[0]["Deactivate"] == True


def test_refactor_ditc():
    input_dict = {0: {"a": 1, "b": 2}, 2: {"a": 3, "b": 4}}
    expected_dict = json.dumps({"a": [1, 3], "b": [2, 4]})

    assert json.dumps(refactor_dict(input_dict)) == expected_dict

    assert refactor_dict({}) == {}
