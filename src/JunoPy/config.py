import os

class JunoPyConfig:
    
    __config = {
        "DATA_DIR": "data/",
        "SPICE_DIR": "spice/",
        "INDEX_DIR": "index/",
        "JADE_DATA_DIR": "data/jade/",
        "JEDI_DATA_DIR": "data/jedi/",
        "WAVES_DATA_DIR": "data/waves/",
        "MAG_DATA_DIR": "data/mag/"
    }

    @staticmethod
    def config(name):
        return JunoPyConfig.__config[name]
    
    @staticmethod
    def set(name, value):
        JunoPyConfig.__config[name] = value

    @staticmethod
    def _update_dependent_configs(name):
        pass