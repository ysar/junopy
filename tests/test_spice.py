
import numpy as np
import datetime as dt

from JunoPy.spice import InitSpiceKernels, get_juno_position

def test_spice():
    """
    """
    InitSpiceKernels()
    times = [
        dt.datetime(2016, 12, 1),
        dt.datetime(2017, 12, 1),
        dt.datetime(2018, 12, 1),
    ]
    pos, _ = get_juno_position(times, coordsystem='JUNO_JSS')
    print(pos)
