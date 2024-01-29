import os
import logging

class ConfigVars:
    
    __config = {"ROOT_DIR": "./"}

    @classmethod
    def get(cls, name):
        return cls.__config[name]
    
    @classmethod
    def set_root(cls, value):
        cls.__config["ROOT_DIR"] = value
        logging.info(f"Root directory changed to {cls.get("ROOT_DIR")}")

        cls._set_dependent_vars()
        

    @classmethod
    def _set_dependent_vars(cls):
        rootDir = cls.__config["ROOT_DIR"]
        dataDir = os.path.join(rootDir, "data")
        cls.__config["DATA_DIR"] = dataDir
        cls.__config["SPICE_DIR"] = os.path.join(rootDir, "spice")
        cls.__config["INDEX_DIR"] = os.path.join(rootDir, "index")
        cls.__config["MAG_DATA_DIR"] = os.path.join(dataDir, "mag")
        cls.__config["JADE_DATA_DIR"] = os.path.join(dataDir, "jade")
        cls.__config["JEDI_DATA_DIR"] = os.path.join(dataDir, "jedi")
        cls.__config["WAVES_SURVEY_DATA_DIR"] = \
            os.path.join(dataDir, "waves-survey")
        cls.__config["TMP_DIR"] = os.path.join(rootDir, "_tmp")

        for key in cls.__config.keys():
            if "DIR" in key:
                os.makedirs(cls.__config[key], exist_ok=True)

        logging.info(f"Data directory changed to {cls.get("DATA_DIR")}")
        logging.info(f"SPICE directory changed to {cls.get("SPICE_DIR")}")
        logging.info(f"Index directory changed to {cls.get("INDEX_DIR")}")

    @classmethod
    def get_all(cls):
        return cls.__config

ConfigVars._set_dependent_vars()