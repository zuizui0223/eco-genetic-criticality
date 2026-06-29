"""Causal-model package compatibility helpers."""

from causal_model.multipatch_criticality_dynamics import SimulationSnapshot


# The simulation state tracks the frequency of the high allele.  Keep the
# generic spelling as a read-only alias for consumers that summarise that same
# state, while retaining ``high_allele_frequency`` as the explicit field name.
if not hasattr(SimulationSnapshot, "allele_frequency"):
    SimulationSnapshot.allele_frequency = property(lambda snapshot: snapshot.high_allele_frequency)
