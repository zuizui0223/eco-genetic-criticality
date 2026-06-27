"""Column-level readiness checks for empirical H1-H3 parameterisation.

The theoretical and simulation layers do not make an empirical claim unless the
necessary field measurements are available.  This module turns that boundary
into an explicit, testable data contract without inventing parameter values.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class HypothesisMeasurementAudit:
    """Available and missing requirements for one empirical hypothesis layer."""

    hypothesis: str
    required_columns: tuple[str, ...]
    available_columns: tuple[str, ...]
    missing_columns: tuple[str, ...]

    @property
    def ready(self) -> bool:
        return not self.missing_columns


@dataclass(frozen=True)
class EmpiricalReadinessAudit:
    """Audits for H1, H2, H3 and their joint readiness."""

    h1: HypothesisMeasurementAudit
    h2: HypothesisMeasurementAudit
    h3: HypothesisMeasurementAudit

    @property
    def all_ready(self) -> bool:
        return self.h1.ready and self.h2.ready and self.h3.ready


H1_REQUIRED_COLUMNS = (
    "patch_id",
    "patch_area",
    "interaction_state",
    "trait_value",
    "performance",
)
H2_REQUIRED_COLUMNS = (
    "patch_id",
    "time",
    "realised_high_trait_abundance",
    "sample_size",
    "high_allele_copies",
)
H3_REQUIRED_COLUMNS = (
    "patch_id",
    "time",
    "census_population",
    "realised_high_trait_abundance",
    "high_allele_copies",
    "sample_size",
    "source_patch_id",
    "destination_patch_id",
    "dispersal_count",
)


def _audit(hypothesis: str, required: tuple[str, ...], available: set[str]) -> HypothesisMeasurementAudit:
    missing = tuple(column for column in required if column not in available)
    return HypothesisMeasurementAudit(
        hypothesis=hypothesis,
        required_columns=required,
        available_columns=tuple(sorted(available)),
        missing_columns=missing,
    )


def audit_empirical_columns(columns: Iterable[str]) -> EmpiricalReadinessAudit:
    """Audit a declared measurement table without assuming its observations.

    ``high_allele_copies`` and ``sample_size`` deliberately require observed
    genotype counts rather than a precomputed allele frequency, so downstream
    finite-sampling uncertainty remains identifiable.
    """
    available = {str(column).strip() for column in columns if str(column).strip()}
    return EmpiricalReadinessAudit(
        h1=_audit("H1", H1_REQUIRED_COLUMNS, available),
        h2=_audit("H2", H2_REQUIRED_COLUMNS, available),
        h3=_audit("H3", H3_REQUIRED_COLUMNS, available),
    )
