from . import config

# Ensure configuration is updated whenever junopy is imported
config.ConfigVars.update_dependent_vars()
