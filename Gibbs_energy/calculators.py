"""Calculator factory helpers.

This module centralizes the creation of the atomistic calculator used by the
rest of the project. Keeping the calculator setup in one place makes it easier
to:

- switch backend in the future;
- add common defaults once instead of repeating them in every workflow;
- document clearly which runtime knobs are exposed to the user.
"""

from mace.calculators import MACECalculator


def setup_calculator(model_path: str, device: str = "cuda"):
    """Create and return the MACE calculator used by all workflows."""
    return MACECalculator(
        model_paths=model_path,
        device=device,
    )
