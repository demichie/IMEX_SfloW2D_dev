"""Execute the tracked ETNA geometry script with a numeric ASCII reader only.

Pandas is used by the original script solely to parse the ESRI numeric table.
This adapter leaves every geometric expression and output formatter unchanged.
Its substitution is recorded in N9 evidence and generated rasters are checked
against the user's previously generated example, when available.
"""
import runpy
import sys
import types
import numpy as np


class NumericTable:
    def __init__(self, values):
        self.values = values

    def astype(self, dtype):
        return NumericTable(self.values.astype(dtype))


def read_table(path, *, delim_whitespace, header, skiprows, dtype):
    assert delim_whitespace and header is None and dtype == 'unicode'
    return NumericTable(np.loadtxt(path, skiprows=skiprows, dtype=str))


parser = types.ModuleType('pandas')
parser.read_table = read_table
sys.modules['pandas'] = parser
sys.path.insert(0, '.')
runpy.run_path('create_input_ellipsoid.py', run_name='__main__')
