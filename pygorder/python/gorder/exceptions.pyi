from builtins import Exception

class AnalysisError(Exception):
    """Exception that can be raised when analyzing the trajectory."""
    ...

class WriteError(Exception):
    """Exception that can be raised when writing the results into output files."""
    ...

class ConfigError(Exception):
    """Exception that can be raised when constructing the config Analysis class."""
    ...

class APIError(Exception):
    """Exception that can be raised when accessing the results programmatically."""
    ...