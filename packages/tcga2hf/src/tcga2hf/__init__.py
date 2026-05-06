"""Typed pydantic models + pyarrow schemas for the
`gabrielaltay/tcga-{patients,tabular}-open` HuggingFace datasets.

The pyarrow `*_FIELDS` lists in `tcga2hf.schema` are the single source of
truth for the dataset shape (regenerated from gdcdictionary YAMLs); the
pydantic models in `tcga2hf.models` are derived from those same lists, so
the two stay in sync by construction. See `TcgaHfPatient` for the typed
patient row + helper joins (tumor/normal pairs, mutations-by-gene,
expression lookup, longitudinal timeline).
"""

from tcga2hf.models import TcgaHfPatient

__version__ = "0.1.0"

__all__ = ["TcgaHfPatient", "__version__"]
