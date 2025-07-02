
import datetime as dt

from junopy import spice
from junopy import config

def test_spice():
    """
    """
    times = [
        dt.datetime(2016, 12, 1),
        dt.datetime(2017, 12, 1),
        dt.datetime(2018, 12, 1),
    ]
    pos, _ = spice.get_juno_position(times, coordsystem='JUNO_JSS')
    print(pos)


# print(config.ConfigVars.get("SPICE_DIR"))
# kernels = spice.identify_juno_kernels()
# spice._remake_meta_kernel_file(kernels)

# times = [
#         dt.datetime(2016, 12, 1),
#         dt.datetime(2017, 12, 1),
#         dt.datetime(2018, 12, 1),
#     ]
# pos, _ = spice.get_juno_position(times, coordsystem='JUNO_JSS')
# print(pos)