"""Matched-total-capacity presets for H3 lifecycle experiments."""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from causal_model.network_h3_lifecycle import isolated_kernel


@dataclass(frozen=True)
class LandscapePreset:
    """A capacity partition and source-to-destination transport kernel."""

    name: str
    capacities: tuple[int, ...]
    kernel: tuple[tuple[float, ...], ...]

    @property
    def total_capacity(self) -> int:
        return sum(self.capacities)


def _equal(total_capacity: int, patches: int) -> tuple[int, ...]:
    if not isinstance(total_capacity, int) or total_capacity < 1:
        raise ValueError("total_capacity must be a positive integer")
    if not isinstance(patches, int) or patches < 1:
        raise ValueError("patches must be a positive integer")
    if total_capacity % patches:
        raise ValueError("total_capacity must divide exactly among patches")
    return tuple(total_capacity // patches for _ in range(patches))


def one_large(total_capacity: int) -> LandscapePreset:
    """One patch holding all capacity."""
    return LandscapePreset("one_large", _equal(total_capacity, 1), ((1.0,),))


def equal_isolated(total_capacity: int, patches: int) -> LandscapePreset:
    """Equal capacity patches with no between-patch transport."""
    capacities = _equal(total_capacity, patches)
    return LandscapePreset("equal_isolated", capacities, isolated_kernel(patches))


def equal_complete_network(
    total_capacity: int,
    patches: int,
    *,
    self_destination_probability: float = 0.0,
) -> LandscapePreset:
    """Equal patches connected symmetrically after emigration."""
    capacities = _equal(total_capacity, patches)
    value = float(self_destination_probability)
    if not isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError("self_destination_probability must lie in [0, 1]")
    if patches == 1:
        kernel = ((1.0,),)
    else:
        off_diagonal = (1.0 - value) / (patches - 1)
        kernel = tuple(
            tuple(value if source == destination else off_diagonal for destination in range(patches))
            for source in range(patches)
        )
    return LandscapePreset("equal_complete_network", capacities, kernel)
