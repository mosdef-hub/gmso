import numpy as np

import topology


def test_element():
    carbon = topology.core.element.carbon

    assert carbon.name == 'carbon'
    assert carbon.symbol == 'C'
    assert np.isclose(carbon.mass, 12.011)