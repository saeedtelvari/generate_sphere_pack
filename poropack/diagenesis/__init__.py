"""
poropack.diagenesis
===================
Post-depositional geological diagenesis modeling:
- apply_compaction: Mechanical uniaxial vertical strain squashing
- apply_cementation: EDT-based syntaxial quartz/calcite pore-throat cementation
- apply_dissolution: Secondary leaching, micro-porosity, and vug creation
"""

from poropack.diagenesis.compaction import apply_compaction
from poropack.diagenesis.cementation import apply_cementation, apply_dissolution, dissolve

__all__ = [
    "apply_compaction",
    "apply_cementation",
    "apply_dissolution",
    "dissolve",
]
