"""Reference pydantic implementation of one TCGA patient row.

`TcgaHfPatient` is the canonical, fully-typed Python representation of one row
from this dataset. It mirrors the parquet schema entity-for-entity and adds
convenience methods for the common joins / transforms that a user would
otherwise have to write themselves: flattening the GDC biospecimen tree,
building tumor/normal sample pairs, indexing mutations by gene, looking up
expression by gene name, sorting events into a longitudinal patient timeline.

## Source of truth

The pyarrow `*_FIELDS` lists in `tcga2hf.schema` (regenerated from the
gdcdictionary YAMLs) are the single source of truth. Every entity model below
is generated from one of those lists via `_make_entity`, so the pydantic field
set is — by construction — identical to the parquet schema. Convenience
methods are added by subclassing each auto-generated base.

## Strictness

Every model inherits `extra="forbid"`. A row with a key not declared in the
GDC dictionary will fail to validate rather than be silently dropped. This
makes the model a precise contract: if a user reaches for `patient.foo` and
`foo` doesn't exist as a typed field, the row would never have validated.

This module is deliberately **slow and explicit** — it's a reference for what
the data means, not the fastest way to operate on it. For high-throughput
workflows, work directly off pyarrow / polars and use these models as the
specification.

Usage:

    from datasets import load_dataset
    from tcga2hf.models import TcgaHfPatient

    ds = load_dataset("gabrielaltay/tcga-patients-open", "TCGA-CHOL", split="train")
    for row in ds:
        patient = TcgaHfPatient.model_validate(row)
        for tx_sample, normal_sample in patient.tumor_normal_pairs():
            print(tx_sample.submitter_id, "vs", normal_sample.submitter_id)
        tp53 = patient.mutations_by_gene().get("TP53", [])
        for variant in tp53:
            print(variant.HGVSp_Short, variant.t_alt_count, variant.t_depth)
        for aliquot_id, expr in patient.expression_for_gene("ALB").items():
            print(aliquot_id, "ALB TPM:", expr["tpm_unstranded"])
        for event in patient.timeline():
            print(event.day, event.category)
"""

from __future__ import annotations

from typing import Any, Literal

import pyarrow as pa
from pydantic import BaseModel, ConfigDict, Field, create_model

from tcga2hf.schema import (
    ALIQUOT_FIELDS,
    ANALYTE_FIELDS,
    DEMOGRAPHIC_FIELDS,
    DIAGNOSIS_FIELDS,
    EXPOSURE_FIELDS,
    EXPRESSION_FIELDS,
    FAMILY_HISTORY_FIELDS,
    FOLLOW_UP_FIELDS,
    MUTATION_FIELDS,
    PATIENT_FIELDS,
    PORTION_FIELDS,
    SAMPLE_FIELDS,
    TREATMENT_FIELDS,
)

# ---------------------------------------------------------------------------
# Base config + pyarrow-to-Python type translation
# ---------------------------------------------------------------------------


class _BaseEntity(BaseModel):
    """Strict base for every TCGA entity model.

    `extra="forbid"` makes the GDC-dictionary field set the contract: any
    unknown key raises a ValidationError rather than being silently dropped.
    """

    model_config = ConfigDict(extra="forbid", frozen=False)


def _pa_to_py_scalar(t: pa.DataType) -> type:
    """Map a pyarrow scalar type to its Python primitive."""
    if pa.types.is_integer(t):
        return int
    if pa.types.is_floating(t):
        return float
    if pa.types.is_boolean(t):
        return bool
    if pa.types.is_string(t):
        return str
    raise TypeError(f"unhandled scalar pyarrow type {t!r}")


def _pa_to_py(t: pa.DataType) -> Any:
    """Translate a pyarrow type to a pydantic-compatible annotation.

    Struct / list-of-struct types must be supplied via `_make_entity` overrides
    (we need the actual entity class, which the type alone can't name).
    """
    if pa.types.is_struct(t):
        raise TypeError("struct types must be supplied via overrides")
    if pa.types.is_list(t):
        if pa.types.is_struct(t.value_type):
            raise TypeError("list-of-struct types must be supplied via overrides")
        # Parquet allows null elements inside list columns; reflect that.
        return list[_pa_to_py_scalar(t.value_type) | None]
    return _pa_to_py_scalar(t)


def _make_entity(
    name: str,
    fields: list[pa.Field],
    overrides: dict[str, tuple[Any, Any]] | None = None,
    required: tuple[str, ...] = (),
    base: type[BaseModel] = _BaseEntity,
) -> type[BaseModel]:
    """Build a pydantic model from a pyarrow `*_FIELDS` list.

    For each pyarrow field:
      - if `overrides[name]` is set, use that (used for nested struct fields
        whose Python type is another generated entity);
      - if `name` is in `required`, the field is required (no default);
      - if it's a list type, it defaults to an empty list;
      - otherwise it's `T | None = None`.
    """
    field_defs: dict[str, tuple[Any, Any]] = {}
    for f in fields:
        if overrides and f.name in overrides:
            field_defs[f.name] = overrides[f.name]
            continue
        py_type = _pa_to_py(f.type)
        if f.name in required:
            field_defs[f.name] = (py_type, ...)
        else:
            # Auto-generated fields are all `T | None = None`. Sparsely
            # populated list-of-scalar fields (e.g.
            # follow_up.imaging_anatomic_site) round-trip from parquet as
            # None when absent, which is distinct from an empty list — so we
            # treat them the same as scalar fields. List-of-struct container
            # fields are passed via `overrides` and default to [] there,
            # because the row builder (`clinical._patient_row`) always emits
            # a list for them.
            field_defs[f.name] = (py_type | None, None)
    return create_model(name, __base__=base, **field_defs)


# ---------------------------------------------------------------------------
# Biospecimen hierarchy: case → samples → portions → analytes → aliquots
# Generated bottom-up so each parent can reference the typed child class.
# ---------------------------------------------------------------------------


Aliquot = _make_entity("Aliquot", ALIQUOT_FIELDS, required=("aliquot_id",))

Analyte = _make_entity(
    "Analyte",
    ANALYTE_FIELDS,
    overrides={"aliquots": (list[Aliquot], Field(default_factory=list))},
    required=("analyte_id",),
)

Portion = _make_entity(
    "Portion",
    PORTION_FIELDS,
    overrides={"analytes": (list[Analyte], Field(default_factory=list))},
    required=("portion_id",),
)

# Sample needs convenience methods, so we generate the field-only base then
# attach helpers via a subclass.
_SampleBase = _make_entity(
    "_SampleBase",
    SAMPLE_FIELDS,
    overrides={"portions": (list[Portion], Field(default_factory=list))},
    required=("sample_id",),
)


class Sample(_SampleBase):
    """One biospecimen sample. Fields are the full GDC dictionary set, inherited
    from `_SampleBase`; the methods below walk the embedded biospecimen tree."""

    def all_analytes(self) -> list[Analyte]:
        return [a for p in self.portions for a in p.analytes]

    def all_aliquots(self) -> list[Aliquot]:
        return [a for p in self.portions for an in p.analytes for a in an.aliquots]

    def aliquots_by_analyte_type(self, analyte_type: str) -> list[Aliquot]:
        """Filter aliquots by their parent analyte's analyte_type ('DNA' / 'RNA')."""
        return [
            a
            for p in self.portions
            for an in p.analytes
            if an.analyte_type == analyte_type
            for a in an.aliquots
        ]

    @property
    def is_tumor(self) -> bool:
        return self.tissue_type == "Tumor"

    @property
    def is_normal(self) -> bool:
        return self.tissue_type == "Normal"


# ---------------------------------------------------------------------------
# Clinical entities
# ---------------------------------------------------------------------------


Demographic = _make_entity("Demographic", DEMOGRAPHIC_FIELDS)
Treatment = _make_entity("Treatment", TREATMENT_FIELDS)
Diagnosis = _make_entity(
    "Diagnosis",
    DIAGNOSIS_FIELDS,
    overrides={"treatments": (list[Treatment], Field(default_factory=list))},
)
FollowUp = _make_entity("FollowUp", FOLLOW_UP_FIELDS)
Exposure = _make_entity("Exposure", EXPOSURE_FIELDS)
FamilyHistory = _make_entity("FamilyHistory", FAMILY_HISTORY_FIELDS)


# ---------------------------------------------------------------------------
# Molecular entities
# ---------------------------------------------------------------------------


# Mutation: one row per MAF variant. 143 fields (3 FK + 140 MAF columns).
# Field types follow the parquet schema exactly; nothing is required here
# because the MAF spec doesn't mandate any single column.
Mutation = _make_entity("Mutation", MUTATION_FIELDS)


_GeneExpressionBase = _make_entity(
    "_GeneExpressionBase",
    EXPRESSION_FIELDS,
    required=("aliquot_id", "source_file_id"),
)


class GeneExpression(_GeneExpressionBase):
    """Per-aliquot RNA-Seq STAR gene quantifications (60,660 genes per record).

    `gene_id` / `gene_name` / `gene_type` and the value arrays are index-aligned;
    use `get_gene` to look up by symbol.
    """

    def get_gene(self, gene_name: str) -> dict[str, Any] | None:
        """Return all per-gene values for `gene_name`, or None if absent.

        Genes are indexed by the order they appeared in the source TSV
        (consistent within a GDC release because the STAR pipeline uses a
        fixed GENCODE v36 reference).
        """
        try:
            i = self.gene_name.index(gene_name)
        except ValueError:
            return None
        return {
            "gene_id": self.gene_id[i],
            "gene_name": self.gene_name[i],
            "gene_type": self.gene_type[i],
            "unstranded": self.unstranded[i],
            "tpm_unstranded": self.tpm_unstranded[i],
            "fpkm_unstranded": self.fpkm_unstranded[i],
            "fpkm_uq_unstranded": self.fpkm_uq_unstranded[i],
        }


# ---------------------------------------------------------------------------
# Patient: top-level row entity with cross-modality joins
# ---------------------------------------------------------------------------


# Every TimelineEvent uses the case's `index_date` as zero. The GDC data
# dictionary defines every `days_to_*` field as "the number of days from the
# index date to the date of <event>". For TCGA the index defaults to date of
# initial pathologic diagnosis. Sources:
#   - sample.days_to_collection / sample.days_to_sample_procurement:
#     https://github.com/NCI-GDC/gdcdictionary src/gdcdictionary/schemas/sample.yaml
#   - case.index_date enum (incl. "Diagnosis"):
#     https://github.com/NCI-GDC/gdcdictionary src/gdcdictionary/schemas/case.yaml
#
# Two sample-related categories per the dictionary's definitions:
#   - `bcr_receipt` ← `sample.days_to_collection` — the GDC dictionary's prose
#     description: "received by the Biospecimen Core Resource (BCR) or other
#     center for processing." Note that the GDC dictionary's NCIt term for the
#     same field is "Biospecimen Collection Date Less Initial Pathologic
#     Diagnosis Date Calculated Day Value" — i.e. the dictionary itself uses
#     two slightly different framings for what event this measures. Some TCGA
#     cases have values here that exceed `days_to_death`; that's surprising but
#     it's what GDC ships, and we don't try to reinterpret it.
#   - `sample_procurement` ← `sample.days_to_sample_procurement` — "the date a
#     patient underwent a procedure (e.g. surgical resection) to yield or
#     remove from the patient a sample that was eventually used for research."
#     Sparsely populated in TCGA.
TimelineCategory = Literal[
    "consent",
    "diagnosis",
    "treatment_start",
    "treatment_end",
    "follow_up",
    "sample_procurement",
    "bcr_receipt",
    "lost_to_followup",
    "death",
]


class TimelineEvent(_BaseEntity):
    """One event on a patient's longitudinal timeline.

    `day` is the GDC `days_to_*` value verbatim, anchored to the case's
    `index_date` (TCGA: date of initial pathologic diagnosis). Every event on
    a `TcgaHfPatient.timeline()` shares this anchor, so they can be sorted
    and differenced directly.
    """

    day: float
    category: TimelineCategory
    label: str | None = None
    sample_id: str | None = None
    diagnosis_id: str | None = None
    treatment_id: str | None = None
    follow_up_id: str | None = None


_TcgaHfPatientBase = _make_entity(
    "_TcgaHfPatientBase",
    PATIENT_FIELDS,
    overrides={
        "demographic": (Demographic | None, None),
        "diagnoses": (list[Diagnosis], Field(default_factory=list)),
        "follow_ups": (list[FollowUp], Field(default_factory=list)),
        "exposures": (list[Exposure], Field(default_factory=list)),
        "family_histories": (list[FamilyHistory], Field(default_factory=list)),
        "samples": (list[Sample], Field(default_factory=list)),
        "samples_masked_somatic_mutation": (list[Mutation], Field(default_factory=list)),
        "samples_gene_expression_quantification": (
            list[GeneExpression],
            Field(default_factory=list),
        ),
    },
    required=("case_id", "case_submitter_id", "project_id"),
)


class TcgaHfPatient(_TcgaHfPatientBase):
    """One TCGA patient row from the `gabrielaltay/tcga-patients-open` dataset.

    Fields are the full GDC `case` set (auto-generated from PATIENT_FIELDS).
    Methods below provide common joins so users don't have to walk the
    nested biospecimen tree by hand.
    """

    # ---- biospecimen joins ----

    def all_samples_by_id(self) -> dict[str, Sample]:
        return {s.sample_id: s for s in self.samples}

    def aliquot_to_sample(self) -> dict[str, str]:
        """Reverse lookup: aliquot_id → sample_id, walking the full GDC tree."""
        out: dict[str, str] = {}
        for s in self.samples:
            for p in s.portions:
                for an in p.analytes:
                    for a in an.aliquots:
                        out[a.aliquot_id] = s.sample_id
        return out

    def aliquot_lookup(self) -> dict[str, tuple[Sample, Portion, Analyte, Aliquot]]:
        """aliquot_id → (sample, portion, analyte, aliquot) — full provenance chain."""
        out: dict[str, tuple[Sample, Portion, Analyte, Aliquot]] = {}
        for s in self.samples:
            for p in s.portions:
                for an in p.analytes:
                    for a in an.aliquots:
                        out[a.aliquot_id] = (s, p, an, a)
        return out

    def tumor_samples(self) -> list[Sample]:
        return [s for s in self.samples if s.is_tumor]

    def normal_samples(self) -> list[Sample]:
        return [s for s in self.samples if s.is_normal]

    def tumor_normal_pairs(self) -> list[tuple[Sample, Sample]]:
        """Distinct (tumor, matched_normal) sample pairs found in mutation calls."""
        by_id = self.all_samples_by_id()
        seen: set[tuple[str, str]] = set()
        pairs: list[tuple[Sample, Sample]] = []
        for v in self.samples_masked_somatic_mutation:
            t = getattr(v, "tumor_sample_id", None)
            n = getattr(v, "matched_normal_sample_id", None)
            if not t or not n or (t, n) in seen:
                continue
            seen.add((t, n))
            tumor = by_id.get(t)
            normal = by_id.get(n)
            if tumor and normal:
                pairs.append((tumor, normal))
        return pairs

    # ---- mutation joins ----

    def mutations_by_gene(self) -> dict[str, list[Mutation]]:
        out: dict[str, list[Mutation]] = {}
        for v in self.samples_masked_somatic_mutation:
            sym = getattr(v, "Hugo_Symbol", None)
            if sym:
                out.setdefault(sym, []).append(v)
        return out

    def mutations_by_consequence(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for v in self.samples_masked_somatic_mutation:
            cls = getattr(v, "Variant_Classification", None) or "Unknown"
            out[cls] = out.get(cls, 0) + 1
        return out

    # ---- expression joins ----

    def expression_for_gene(self, gene_name: str) -> dict[str, dict[str, Any]]:
        """{aliquot_id: per-gene value dict} for `gene_name` across RNA-Seq aliquots."""
        out: dict[str, dict[str, Any]] = {}
        for er in self.samples_gene_expression_quantification:
            entry = er.get_gene(gene_name)
            if entry is not None:
                out[er.aliquot_id] = entry
        return out

    # ---- timeline consistency ----

    def consistency_check(self) -> dict[str, int]:
        """Counts of dated-event patterns worth surfacing for downstream consumers.

        These are not assertions of correctness — they're observations that
        consumers might want to filter on or investigate:

          - `bcr_receipts_after_death`: count of samples where
            `days_to_collection` exceeds `demographic.days_to_death`. Per the
            GDC dictionary both are anchored to the case's index date; values
            here mean GDC reports a sample-level date later than the patient's
            recorded death. We don't reinterpret it; flag and let consumers
            (or a human spot-check on the GDC Data Portal / cBioPortal) decide.
          - `samples_with_no_temporal_data`: samples missing both
            `days_to_sample_procurement` and `days_to_collection`. Can't be
            placed on the timeline at all.
          - `pre_index_treatments`: treatments with negative
            `days_to_treatment_start`.
          - `non_diagnosis_index`: 1 if `index_date` is set to something other
            than "Diagnosis" (the TCGA default). Affects what zero means on
            the timeline for this patient.
        """
        counts = {
            "bcr_receipts_after_death": 0,
            "samples_with_no_temporal_data": 0,
            "pre_index_treatments": 0,
            "non_diagnosis_index": 0,
        }
        days_to_death = self.demographic.days_to_death if self.demographic else None
        for s in self.samples:
            if s.days_to_sample_procurement is None and s.days_to_collection is None:
                counts["samples_with_no_temporal_data"] += 1
            if (
                s.days_to_collection is not None
                and days_to_death is not None
                and s.days_to_collection > days_to_death
            ):
                counts["bcr_receipts_after_death"] += 1
        for dx in self.diagnoses:
            for tx in dx.treatments:
                if tx.days_to_treatment_start is not None and tx.days_to_treatment_start < 0:
                    counts["pre_index_treatments"] += 1
        if self.index_date is not None and self.index_date != "Diagnosis":
            counts["non_diagnosis_index"] = 1
        return counts

    # ---- longitudinal timeline ----

    def timeline(self) -> list[TimelineEvent]:
        """Every dated event for this patient, sorted ascending by `day`.

        All events use the same anchor: the case's `index_date` (TCGA: date
        of initial pathologic diagnosis). Each event's `day` is the GDC
        `days_to_*` value verbatim, so positive = after index, negative =
        before. Categories included:

        Clinical (events documented as happening to the patient):
          - `consent`           — `case.days_to_consent`
          - `diagnosis`         — `diagnosis.days_to_diagnosis`
          - `treatment_start`   — `treatment.days_to_treatment_start`
          - `treatment_end`     — `treatment.days_to_treatment_end`
          - `follow_up`         — `follow_up.days_to_follow_up`
          - `lost_to_followup`  — `case.days_to_lost_to_followup`
          - `death`             — `demographic.days_to_death`

        Biospecimen (events documented as involving the sample):
          - `sample_procurement` ← `sample.days_to_sample_procurement`
          - `bcr_receipt`        ← `sample.days_to_collection`

        See the module-level comment block for the dictionary definitions of
        these last two and a note about the `days_to_collection` semantics.
        """
        events: list[TimelineEvent] = []

        if self.days_to_consent is not None:
            events.append(
                TimelineEvent(
                    day=float(self.days_to_consent),
                    category="consent",
                    label=self.consent_type,
                )
            )

        for dx in self.diagnoses:
            if dx.days_to_diagnosis is not None:
                events.append(
                    TimelineEvent(
                        day=float(dx.days_to_diagnosis),
                        category="diagnosis",
                        label=dx.primary_diagnosis,
                        diagnosis_id=dx.diagnosis_id,
                    )
                )
            for tx in dx.treatments:
                if tx.days_to_treatment_start is not None:
                    events.append(
                        TimelineEvent(
                            day=float(tx.days_to_treatment_start),
                            category="treatment_start",
                            label=tx.treatment_type,
                            diagnosis_id=dx.diagnosis_id,
                            treatment_id=tx.treatment_id,
                        )
                    )
                if tx.days_to_treatment_end is not None:
                    events.append(
                        TimelineEvent(
                            day=float(tx.days_to_treatment_end),
                            category="treatment_end",
                            label=tx.treatment_type,
                            diagnosis_id=dx.diagnosis_id,
                            treatment_id=tx.treatment_id,
                        )
                    )

        for fu in self.follow_ups:
            if fu.days_to_follow_up is not None:
                events.append(
                    TimelineEvent(
                        day=float(fu.days_to_follow_up),
                        category="follow_up",
                        label=fu.disease_response,
                        follow_up_id=fu.follow_up_id,
                    )
                )

        for s in self.samples:
            if s.days_to_sample_procurement is not None:
                events.append(
                    TimelineEvent(
                        day=float(s.days_to_sample_procurement),
                        category="sample_procurement",
                        label=s.sample_type,
                        sample_id=s.sample_id,
                    )
                )
            if s.days_to_collection is not None:
                events.append(
                    TimelineEvent(
                        day=float(s.days_to_collection),
                        category="bcr_receipt",
                        label=s.sample_type,
                        sample_id=s.sample_id,
                    )
                )

        if self.days_to_lost_to_followup is not None:
            events.append(
                TimelineEvent(
                    day=float(self.days_to_lost_to_followup),
                    category="lost_to_followup",
                    label=self.lost_to_followup,
                )
            )

        if self.demographic and self.demographic.days_to_death is not None:
            events.append(
                TimelineEvent(
                    day=float(self.demographic.days_to_death),
                    category="death",
                    label=self.demographic.vital_status,
                )
            )

        events.sort(key=lambda e: e.day)
        return events
