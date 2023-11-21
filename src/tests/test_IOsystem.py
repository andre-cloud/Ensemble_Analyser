import numpy as np
import pytest, os, json

from src.IOsystem import _parse_xyz_str, mkdir, SerialiseEncoder

def test_parse_xyz_str():
    with open('tests/assets/geom.xyz') as f:
        inp = f.readlines()

    exp_atom = np.array([ "N","N","C","C","O","O","C","O","O","C","H","H","H","H","H","H"])
    exp_geom = np.array([
        [-0.4327    ,   -1.0281    ,   -0.3534],
        [0.3211     ,   -1.9808    ,    0.3099],
        [-1.6507    ,   -0.6978    ,   0.1984],
        [1.3562     ,   -2.5953    ,    -0.3400],
        [1.3760     ,   -2.2144    ,    -1.6218],
        [2.0938     ,   -3.3942    ,    0.1845],
        [2.5341     ,   -2.4677    ,    -2.4541],
        [-2.5896    ,   -0.6011    ,   -0.7422],
        [-1.8017    ,   -0.5880    ,   1.3898],
        [-3.9639    ,   -0.3032    ,   -0.3819],
        [-4.0058    ,   0.0416     ,    0.7039],
        [-4.3636    ,   0.5227     ,    -1.0585],
        [-4.6014    ,   -1.2397    ,     -0.5081],
        [3.2659     ,   -1.5986    ,    -2.3616],
        [2.2054     ,   -2.5738    ,    -3.5405],
        [3.0424     ,   -3.4302    ,    -2.1152],
    ])

    out = _parse_xyz_str(inp)
    assert list(out[0]) == list(exp_atom)
    assert np.allclose(out[1], exp_geom)


def test_mkdir(tmpdir):
    temp_dir = str(tmpdir.mkdir("test_directory"))

    directory = os.path.join(temp_dir, "new_directory")

    mkdir(directory)
    assert os.path.exists(directory)

    with pytest.raises(IOError):
        mkdir(directory)

    assert os.path.exists(directory)


def test_serialise_encoder():
    encoder = SerialiseEncoder()

    arr = np.array([1, 2, 3])
    expected_output = "[1, 2, 3]"
    assert encoder.encode(arr) == expected_output

    obj = {'a': 1, 'b': 2}
    expected_output = '{"a": 1, "b": 2}'
    assert encoder.encode(obj) == expected_output