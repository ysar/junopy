
import datetime as dt

from junopy import spice
from junopy import config

def test_spice():
    """
    This test checks if SPICE kernels are loaded and Juno position can be retrieved.
    """
    times = [
        dt.datetime(2016, 12, 1),
        dt.datetime(2017, 12, 1),
        dt.datetime(2018, 12, 1),
    ]
    pos, _ = spice.get_juno_position(times, coordsystem='JUNO_JSS')
    print(pos)
