"""Shared physical constants.

Only project-wide constants should live here. This keeps numerical conversion
factors out of the scientific logic and makes later auditing easier.
"""

# Conversion used when adding a pV contribution to thermodynamic quantities.
PRESSURE_GPA_TO_EV_A3 = 0.0062415091
