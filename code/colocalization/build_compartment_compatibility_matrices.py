#!/usr/bin/env python3
"""Build reviewable compartment-compatibility (C) matrices.

The matrices encode *capacity for physical interaction*, not a posterior
probability and not proof that an interaction occurs.  Native location labels
from HPA and COMPARTMENTS are projected into an explicit canonical physical
topology.  The mpkCCD five-fraction assay is handled separately because its
columns are centrifugation fractions rather than biological compartments.

All weights and mappings are provisional and are written as flat audit tables
so a domain scientist can accept, reject, or revise them before the matrices
are used in the Bayesian colocalization graph.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


PROJECT_DEFAULT = Path(__file__).resolve().parents[2]
GO_RELEASE = "2026-06-19"
GO_URL = (
    "https://release.geneontology.org/"
    f"{GO_RELEASE}/ontology/go-basic.obo"
)
HPA_URL = "https://www.proteinatlas.org/humanproteome/subcellular"
HPA_DATA_URL = "https://www.proteinatlas.org/humanproteome/subcellular/data"
COMPARTMENTS_URL = "https://compartments.jensenlab.org/Downloads"
MPKCCD_URL = (
    "https://esbl.nhlbi.nih.gov/Databases/mpkFractions/"
    "proteomic_fractions_log_files/Proteomics_of_subcellular_fractions.xlsx"
)
TOPOLOGY_URL = "https://www.ncbi.nlm.nih.gov/books/NBK26907/"
NUCLEAR_URL = "https://www.ncbi.nlm.nih.gov/books/NBK26932/"
MITOCHONDRIA_URL = "https://www.ncbi.nlm.nih.gov/books/NBK26894/"
ER_URL = "https://www.ncbi.nlm.nih.gov/books/NBK26841/"


CANONICAL_ZONES = [
    "cytosol_bulk",
    "cytoskeleton",
    "plasma_membrane",
    "extracellular_space",
    "cell_surface_junction",
    "nucleus_interior",
    "nuclear_envelope",
    "er_lumen",
    "er_membrane",
    "golgi_lumen",
    "golgi_membrane",
    "endosome_lumen",
    "endosome_membrane",
    "lysosome_lumen",
    "lysosome_membrane",
    "vesicle_lumen",
    "vesicle_membrane",
    "mitochondrial_matrix",
    "mitochondrial_inner_membrane",
    "mitochondrial_intermembrane_space",
    "mitochondrial_outer_membrane",
    "peroxisome_matrix",
    "peroxisome_membrane",
    "cilium_interior",
    "ciliary_membrane",
    "lipid_droplet_surface",
    "cytoplasm_general",
    "generic_membrane",
    "intracellular_unspecified",
    "sperm_specialized_structure",
]


ZONE_DEFINITIONS = {
    "cytosol_bulk": "Aqueous cytosolic phase.",
    "cytoskeleton": "Cytosolic structural networks and associated assemblies.",
    "plasma_membrane": "Plasma-membrane bilayer; protein orientation is unresolved.",
    "extracellular_space": "Extracellular fluid and extracellular matrix.",
    "cell_surface_junction": "Cell junction, focal-adhesion, and cortical membrane interface.",
    "nucleus_interior": "Nucleoplasm and non-membrane-bounded nuclear substructures.",
    "nuclear_envelope": "Inner/outer nuclear membranes and pore-associated boundary.",
    "er_lumen": "Endoplasmic-reticulum lumen.",
    "er_membrane": "Endoplasmic-reticulum membrane; face is unresolved.",
    "golgi_lumen": "Golgi cisternal lumen.",
    "golgi_membrane": "Golgi membrane; face is unresolved.",
    "endosome_lumen": "Endosomal lumen.",
    "endosome_membrane": "Endosomal membrane; face is unresolved.",
    "lysosome_lumen": "Lysosomal lumen.",
    "lysosome_membrane": "Lysosomal membrane; face is unresolved.",
    "vesicle_lumen": "Generic vesicle, granule, autophagosome, or phagosome lumen.",
    "vesicle_membrane": "Generic vesicle, granule, autophagosome, or phagosome membrane.",
    "mitochondrial_matrix": "Mitochondrial matrix.",
    "mitochondrial_inner_membrane": "Mitochondrial inner membrane.",
    "mitochondrial_intermembrane_space": "Mitochondrial intermembrane space.",
    "mitochondrial_outer_membrane": "Mitochondrial outer membrane.",
    "peroxisome_matrix": "Peroxisomal lumen/matrix.",
    "peroxisome_membrane": "Peroxisomal membrane.",
    "cilium_interior": "Ciliary interior, axoneme, base, or transition-zone interior.",
    "ciliary_membrane": "Ciliary membrane.",
    "lipid_droplet_surface": "Lipid-droplet phospholipid monolayer/cytosolic surface.",
    "cytoplasm_general": "Broad cytoplasm annotation without resolved subcompartment.",
    "generic_membrane": "Membrane annotation without an identified organelle.",
    "intracellular_unspecified": "Broad intracellular annotation without resolved topology.",
    "sperm_specialized_structure": "Sperm-specific structure with insufficient general topology.",
}


@dataclass(frozen=True)
class GoTerm:
    go_id: str
    name: str
    namespace: str
    parents: tuple[str, ...]
    alt_ids: tuple[str, ...]
    obsolete: bool


@dataclass(frozen=True)
class Anchor:
    go_id: str
    zones: dict[str, float]
    specificity: str
    description: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_DEFAULT)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize(weights: dict[str, float]) -> dict[str, float]:
    clean = {
        key: float(value)
        for key, value in weights.items()
        if key in CANONICAL_ZONES and float(value) > 0
    }
    total = sum(clean.values())
    if total <= 0:
        return {}
    return {key: value / total for key, value in clean.items()}


def parse_obo(path: Path) -> tuple[dict[str, GoTerm], dict[str, str]]:
    terms: dict[str, GoTerm] = {}
    alt_to_primary: dict[str, str] = {}
    current: dict[str, object] | None = None

    def finish(block: dict[str, object] | None) -> None:
        if not block or not block.get("id"):
            return
        term = GoTerm(
            go_id=str(block["id"]),
            name=str(block.get("name", "")),
            namespace=str(block.get("namespace", "")),
            parents=tuple(block.get("parents", [])),
            alt_ids=tuple(block.get("alt_ids", [])),
            obsolete=bool(block.get("obsolete", False)),
        )
        terms[term.go_id] = term
        for alt_id in term.alt_ids:
            alt_to_primary[alt_id] = term.go_id

    with path.open("rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "[Term]":
                finish(current)
                current = {
                    "parents": [],
                    "alt_ids": [],
                    "obsolete": False,
                }
                continue
            if line.startswith("[") and line.endswith("]"):
                finish(current)
                current = None
                continue
            if current is None or not line:
                continue
            if line.startswith("id: "):
                current["id"] = line[4:].strip()
            elif line.startswith("name: "):
                current["name"] = line[6:].strip()
            elif line.startswith("namespace: "):
                current["namespace"] = line[11:].strip()
            elif line.startswith("alt_id: "):
                current["alt_ids"].append(line[8:].strip())
            elif line.startswith("is_a: "):
                current["parents"].append(line[6:].split(" ! ", 1)[0].strip())
            elif line.startswith("relationship: part_of "):
                parent = line[len("relationship: part_of ") :].split(" ! ", 1)[0]
                current["parents"].append(parent.strip())
            elif line == "is_obsolete: true":
                current["obsolete"] = True
    finish(current)
    return terms, alt_to_primary


def make_anchors() -> dict[str, Anchor]:
    raw = [
        ("GO:0005829", {"cytosol_bulk": 1}, "specific", "cytosol"),
        ("GO:0005856", {"cytoskeleton": 1}, "specific", "cytoskeleton"),
        ("GO:0005813", {"cytoskeleton": 1}, "specific", "centrosome"),
        ("GO:0030054", {"cell_surface_junction": 1}, "specific", "cell junction"),
        ("GO:0005925", {"cell_surface_junction": 0.6, "cytoskeleton": 0.4}, "specific", "focal adhesion"),
        ("GO:0005886", {"plasma_membrane": 1}, "specific", "plasma membrane"),
        ("GO:0005576", {"extracellular_space": 1}, "specific", "extracellular region"),
        ("GO:0031012", {"extracellular_space": 1}, "specific", "extracellular matrix"),
        ("GO:0005654", {"nucleus_interior": 1}, "specific", "nucleoplasm"),
        ("GO:0005635", {"nuclear_envelope": 1}, "specific", "nuclear envelope"),
        ("GO:0031965", {"nuclear_envelope": 1}, "specific", "nuclear membrane"),
        ("GO:0005634", {"nucleus_interior": 0.85, "nuclear_envelope": 0.15}, "coarse", "nucleus"),
        ("GO:0005788", {"er_lumen": 1}, "specific", "ER lumen"),
        ("GO:0005789", {"er_membrane": 1}, "specific", "ER membrane"),
        ("GO:0005783", {"er_lumen": 0.5, "er_membrane": 0.5}, "coarse", "ER"),
        ("GO:0005796", {"golgi_lumen": 1}, "specific", "Golgi lumen"),
        ("GO:0000139", {"golgi_membrane": 1}, "specific", "Golgi membrane"),
        ("GO:0005794", {"golgi_lumen": 0.5, "golgi_membrane": 0.5}, "coarse", "Golgi apparatus"),
        ("GO:0031904", {"endosome_lumen": 1}, "specific", "endosome lumen"),
        ("GO:0010008", {"endosome_membrane": 1}, "specific", "endosome membrane"),
        ("GO:0005768", {"endosome_lumen": 0.5, "endosome_membrane": 0.5}, "coarse", "endosome"),
        ("GO:0043202", {"lysosome_lumen": 1}, "specific", "lysosome lumen"),
        ("GO:0005765", {"lysosome_membrane": 1}, "specific", "lysosome membrane"),
        ("GO:0005764", {"lysosome_lumen": 0.5, "lysosome_membrane": 0.5}, "coarse", "lysosome"),
        ("GO:0031983", {"vesicle_lumen": 1}, "specific", "vesicle lumen"),
        ("GO:0012506", {"vesicle_membrane": 1}, "specific", "vesicle membrane"),
        ("GO:0031982", {"vesicle_lumen": 0.5, "vesicle_membrane": 0.5}, "coarse", "vesicle"),
        ("GO:0031410", {"vesicle_lumen": 0.5, "vesicle_membrane": 0.5}, "coarse", "cytoplasmic vesicle"),
        ("GO:0005759", {"mitochondrial_matrix": 1}, "specific", "mitochondrial matrix"),
        ("GO:0005743", {"mitochondrial_inner_membrane": 1}, "specific", "mitochondrial inner membrane"),
        ("GO:0005758", {"mitochondrial_intermembrane_space": 1}, "specific", "mitochondrial intermembrane space"),
        ("GO:0005741", {"mitochondrial_outer_membrane": 1}, "specific", "mitochondrial outer membrane"),
        (
            "GO:0005739",
            {
                "mitochondrial_matrix": 0.25,
                "mitochondrial_inner_membrane": 0.25,
                "mitochondrial_intermembrane_space": 0.25,
                "mitochondrial_outer_membrane": 0.25,
            },
            "coarse",
            "mitochondrion",
        ),
        ("GO:0005782", {"peroxisome_matrix": 1}, "specific", "peroxisomal matrix"),
        ("GO:0005778", {"peroxisome_membrane": 1}, "specific", "peroxisomal membrane"),
        ("GO:0005777", {"peroxisome_matrix": 0.5, "peroxisome_membrane": 0.5}, "coarse", "peroxisome"),
        ("GO:0060170", {"ciliary_membrane": 1}, "specific", "ciliary membrane"),
        ("GO:0097014", {"cilium_interior": 1}, "specific", "ciliary plasm"),
        ("GO:0005929", {"cilium_interior": 0.6, "ciliary_membrane": 0.4}, "coarse", "cilium"),
        ("GO:0005811", {"lipid_droplet_surface": 1}, "specific", "lipid droplet"),
        ("GO:0005737", {"cytoplasm_general": 1}, "generic", "cytoplasm"),
        ("GO:0016020", {"generic_membrane": 1}, "generic", "membrane"),
        ("GO:0005622", {"intracellular_unspecified": 1}, "generic", "intracellular anatomical structure"),
    ]
    return {
        go_id: Anchor(go_id, normalize(zones), specificity, description)
        for go_id, zones, specificity, description in raw
    }


def canonical_rules() -> tuple[np.ndarray, pd.DataFrame]:
    index = {zone: position for position, zone in enumerate(CANONICAL_ZONES)}
    matrix = np.zeros((len(CANONICAL_ZONES), len(CANONICAL_ZONES)), dtype=float)
    records: list[dict[str, object]] = []

    def add(
        a: str,
        b: str,
        weight: float,
        relation: str,
        rationale: str,
        source_url: str,
        review_status: str = "provisional_scientist_review_required",
    ) -> None:
        ia, ib = index[a], index[b]
        matrix[ia, ib] = weight
        matrix[ib, ia] = weight
        records.append(
            {
                "zone_a": a,
                "zone_b": b,
                "compatibility_weight": weight,
                "relation_class": relation,
                "rationale": rationale,
                "source_url": source_url,
                "review_status": review_status,
            }
        )

    for zone in CANONICAL_ZONES:
        add(
            zone,
            zone,
            1.0,
            "same_canonical_zone",
            "Same resolved physical zone.",
            TOPOLOGY_URL,
        )

    interfaces = [
        ("cytosol_bulk", "plasma_membrane", "Cytosolic face of plasma membrane.", HPA_URL),
        ("extracellular_space", "plasma_membrane", "Extracellular face of plasma membrane.", HPA_URL),
        ("plasma_membrane", "cell_surface_junction", "Junctions and adhesions are membrane-associated interfaces.", HPA_URL),
        ("cytoskeleton", "cell_surface_junction", "Cortical cytoskeleton anchors cell-surface junctions.", HPA_URL),
        ("cytosol_bulk", "nuclear_envelope", "Cytosolic face and nuclear pores.", NUCLEAR_URL),
        ("nucleus_interior", "nuclear_envelope", "Nucleoplasmic face and nuclear pores.", NUCLEAR_URL),
        ("cytosol_bulk", "er_membrane", "Cytosolic face of ER membrane.", ER_URL),
        ("er_lumen", "er_membrane", "Luminal face of ER membrane.", ER_URL),
        ("cytosol_bulk", "golgi_membrane", "Cytosolic face of Golgi membrane.", TOPOLOGY_URL),
        ("golgi_lumen", "golgi_membrane", "Luminal face of Golgi membrane.", TOPOLOGY_URL),
        ("cytosol_bulk", "endosome_membrane", "Cytosolic face of endosomal membrane.", TOPOLOGY_URL),
        ("endosome_lumen", "endosome_membrane", "Luminal face of endosomal membrane.", TOPOLOGY_URL),
        ("cytosol_bulk", "lysosome_membrane", "Cytosolic face of lysosomal membrane.", TOPOLOGY_URL),
        ("lysosome_lumen", "lysosome_membrane", "Luminal face of lysosomal membrane.", TOPOLOGY_URL),
        ("cytosol_bulk", "vesicle_membrane", "Cytosolic face of vesicle membrane.", TOPOLOGY_URL),
        ("vesicle_lumen", "vesicle_membrane", "Luminal face of vesicle membrane.", TOPOLOGY_URL),
        ("cytosol_bulk", "mitochondrial_outer_membrane", "Cytosolic face of mitochondrial outer membrane.", MITOCHONDRIA_URL),
        (
            "mitochondrial_outer_membrane",
            "mitochondrial_intermembrane_space",
            "Inner face of mitochondrial outer membrane.",
            MITOCHONDRIA_URL,
        ),
        (
            "mitochondrial_intermembrane_space",
            "mitochondrial_inner_membrane",
            "Outer face of mitochondrial inner membrane.",
            MITOCHONDRIA_URL,
        ),
        (
            "mitochondrial_inner_membrane",
            "mitochondrial_matrix",
            "Matrix face of mitochondrial inner membrane.",
            MITOCHONDRIA_URL,
        ),
        ("cytosol_bulk", "peroxisome_membrane", "Cytosolic face of peroxisomal membrane.", TOPOLOGY_URL),
        ("peroxisome_matrix", "peroxisome_membrane", "Matrix face of peroxisomal membrane.", TOPOLOGY_URL),
        ("cilium_interior", "ciliary_membrane", "Ciliary membrane bounds ciliary interior.", HPA_URL),
        ("cytosol_bulk", "lipid_droplet_surface", "Lipid droplets expose a cytosolic phospholipid surface.", TOPOLOGY_URL),
    ]
    for a, b, rationale, url in interfaces:
        add(a, b, 0.7, "direct_physical_interface", rationale, url)

    add(
        "plasma_membrane",
        "cell_surface_junction",
        0.8,
        "direct_specialized_interface",
        "Cell junctions are specialized plasma-membrane interfaces.",
        HPA_URL,
    )
    add(
        "cytosol_bulk",
        "cytoskeleton",
        0.6,
        "shared_bulk_neighborhood",
        "Cytoskeletal structures are embedded in the cytosolic phase.",
        HPA_URL,
    )
    add(
        "cytosol_bulk",
        "cell_surface_junction",
        0.6,
        "shared_bulk_neighborhood",
        "Cytosolic proteins can access intracellular junctional complexes.",
        HPA_URL,
    )
    add(
        "extracellular_space",
        "cell_surface_junction",
        0.5,
        "shared_interface_neighborhood",
        "Extracellular domains and matrix contact cell-surface junctions.",
        HPA_URL,
    )
    add(
        "cytoplasm_general",
        "cytosol_bulk",
        0.6,
        "coarse_annotation_overlap",
        "Cytoplasm is broader than cytosol and does not resolve organelles.",
        TOPOLOGY_URL,
    )
    for membrane in [
        "plasma_membrane",
        "nuclear_envelope",
        "er_membrane",
        "golgi_membrane",
        "endosome_membrane",
        "lysosome_membrane",
        "vesicle_membrane",
        "mitochondrial_outer_membrane",
        "mitochondrial_inner_membrane",
        "peroxisome_membrane",
        "ciliary_membrane",
    ]:
        add(
            "generic_membrane",
            membrane,
            0.5,
            "coarse_annotation_overlap",
            "Generic membrane annotation overlaps this resolved membrane class.",
            TOPOLOGY_URL,
        )
    add(
        "generic_membrane",
        "cytosol_bulk",
        0.25,
        "orientation_unresolved",
        "A generic membrane may expose a cytosolic face, but organelle and orientation are unknown.",
        TOPOLOGY_URL,
    )
    for zone in [
        "cytosol_bulk",
        "cytoskeleton",
        "nucleus_interior",
        "nuclear_envelope",
        "er_lumen",
        "er_membrane",
        "golgi_lumen",
        "golgi_membrane",
        "endosome_lumen",
        "endosome_membrane",
        "lysosome_lumen",
        "lysosome_membrane",
        "vesicle_lumen",
        "vesicle_membrane",
        "mitochondrial_matrix",
        "mitochondrial_inner_membrane",
        "mitochondrial_intermembrane_space",
        "mitochondrial_outer_membrane",
        "peroxisome_matrix",
        "peroxisome_membrane",
        "cilium_interior",
        "ciliary_membrane",
        "lipid_droplet_surface",
    ]:
        add(
            "intracellular_unspecified",
            zone,
            0.15,
            "weak_coarse_compatibility",
            "Broad intracellular annotation supplies only weak compatibility.",
            TOPOLOGY_URL,
        )

    transient = [
        ("nucleus_interior", "cytosol_bulk", 0.25, "Selective nucleocytoplasmic shuttling.", NUCLEAR_URL),
        ("nuclear_envelope", "er_membrane", 0.4, "Outer nuclear membrane is continuous with ER membrane.", NUCLEAR_URL),
        ("nuclear_envelope", "er_lumen", 0.2, "Perinuclear space is continuous with ER lumen; face is unresolved.", NUCLEAR_URL),
        ("er_membrane", "vesicle_membrane", 0.25, "Secretory-pathway vesicular continuity.", TOPOLOGY_URL),
        ("er_lumen", "vesicle_lumen", 0.25, "Secretory cargo can transit from ER lumen in vesicles.", TOPOLOGY_URL),
        ("golgi_membrane", "vesicle_membrane", 0.30, "Golgi exchanges membrane through transport vesicles.", TOPOLOGY_URL),
        ("golgi_lumen", "vesicle_lumen", 0.30, "Golgi luminal cargo transits in vesicles.", TOPOLOGY_URL),
        ("endosome_membrane", "vesicle_membrane", 0.30, "Endosomal membrane exchanges via vesicular traffic.", TOPOLOGY_URL),
        ("endosome_lumen", "vesicle_lumen", 0.30, "Endosomal cargo exchanges via vesicular traffic.", TOPOLOGY_URL),
        ("lysosome_membrane", "endosome_membrane", 0.30, "Endosome-to-lysosome maturation/traffic.", TOPOLOGY_URL),
        ("lysosome_lumen", "endosome_lumen", 0.30, "Endosome-to-lysosome cargo continuity.", TOPOLOGY_URL),
        ("lysosome_membrane", "vesicle_membrane", 0.25, "Lysosomal membrane participates in vesicular traffic.", TOPOLOGY_URL),
        ("lysosome_lumen", "vesicle_lumen", 0.25, "Lysosomal cargo participates in vesicular traffic.", TOPOLOGY_URL),
        ("plasma_membrane", "vesicle_membrane", 0.30, "Endocytosis and exocytosis exchange membrane.", TOPOLOGY_URL),
        ("extracellular_space", "vesicle_lumen", 0.20, "Exocytosis makes vesicle lumen topologically extracellular.", TOPOLOGY_URL),
        ("golgi_membrane", "er_membrane", 0.20, "ER-Golgi exchange is mediated by transport carriers.", TOPOLOGY_URL),
        ("golgi_lumen", "er_lumen", 0.20, "Luminal cargo moves between ER and Golgi.", TOPOLOGY_URL),
        ("cilium_interior", "cytosol_bulk", 0.25, "Ciliary base permits selective exchange with cytoplasm.", HPA_URL),
        ("ciliary_membrane", "plasma_membrane", 0.35, "Ciliary membrane is continuous with plasma membrane.", HPA_URL),
        ("lipid_droplet_surface", "er_membrane", 0.30, "Lipid droplets originate from and contact ER membrane.", TOPOLOGY_URL),
    ]
    for a, b, weight, rationale, url in transient:
        add(a, b, weight, "transient_or_topological_continuity", rationale, url)

    for zone in ["cytosol_bulk", "cytoskeleton", "cilium_interior"]:
        add(
            "sperm_specialized_structure",
            zone,
            0.25,
            "low_confidence_specialized_structure",
            "Sperm-specific HPA structure is insufficiently resolved in the general topology.",
            HPA_URL,
        )

    return matrix, pd.DataFrame(records)


def manual_hpa_mapping(location: str) -> tuple[dict[str, float], str, str]:
    one = {
        "Aggresome": "cytosol_bulk",
        "Cytoplasmic bodies": "cytosol_bulk",
        "Cytosol": "cytosol_bulk",
        "Rods & Rings": "cytosol_bulk",
        "Actin filaments": "cytoskeleton",
        "Centriolar satellite": "cytoskeleton",
        "Centrosome": "cytoskeleton",
        "Cytokinetic bridge": "cytoskeleton",
        "Flagellar centriole": "cytoskeleton",
        "Intermediate filaments": "cytoskeleton",
        "Microtubule ends": "cytoskeleton",
        "Microtubules": "cytoskeleton",
        "Midbody": "cytoskeleton",
        "Midbody ring": "cytoskeleton",
        "Mitotic spindle": "cytoskeleton",
        "Kinetochore": "nucleus_interior",
        "Mitotic chromosome": "nucleus_interior",
        "Nuclear bodies": "nucleus_interior",
        "Nuclear speckles": "nucleus_interior",
        "Nucleoli": "nucleus_interior",
        "Nucleoli fibrillar center": "nucleus_interior",
        "Nucleoli rim": "nucleus_interior",
        "Nucleoplasm": "nucleus_interior",
        "Nuclear membrane": "nuclear_envelope",
        "Cell Junctions": "cell_surface_junction",
        "Lipid droplets": "lipid_droplet_surface",
        "Plasma membrane": "plasma_membrane",
        "Primary cilium tip": "cilium_interior",
    }
    if location in one:
        return {one[location]: 1.0}, "manual_specific_mapping", "high"
    maps = {
        "Acrosome": ({"vesicle_lumen": 0.5, "vesicle_membrane": 0.5}, "manual_coarse_organelle_mapping", "moderate"),
        "Annulus": ({"sperm_specialized_structure": 1.0}, "manual_specialized_mapping", "low"),
        "Basal body": ({"cytoskeleton": 0.7, "cilium_interior": 0.3}, "manual_interface_mapping", "moderate"),
        "Calyx": ({"sperm_specialized_structure": 1.0}, "manual_specialized_mapping", "low"),
        "Cleavage furrow": ({"cytoskeleton": 0.5, "cell_surface_junction": 0.5}, "manual_interface_mapping", "moderate"),
        "Connecting piece": ({"sperm_specialized_structure": 0.5, "cytoskeleton": 0.5}, "manual_specialized_mapping", "low"),
        "End piece": ({"sperm_specialized_structure": 0.5, "cilium_interior": 0.5}, "manual_specialized_mapping", "low"),
        "Endoplasmic reticulum": ({"er_lumen": 0.5, "er_membrane": 0.5}, "manual_coarse_organelle_mapping", "moderate"),
        "Endosomes": ({"endosome_lumen": 0.5, "endosome_membrane": 0.5}, "manual_coarse_organelle_mapping", "moderate"),
        "Equatorial segment": ({"sperm_specialized_structure": 1.0}, "manual_specialized_mapping", "low"),
        "Focal adhesion sites": ({"cell_surface_junction": 0.6, "cytoskeleton": 0.4}, "manual_interface_mapping", "moderate"),
        "Golgi apparatus": ({"golgi_lumen": 0.5, "golgi_membrane": 0.5}, "manual_coarse_organelle_mapping", "moderate"),
        "Lysosomes": ({"lysosome_lumen": 0.5, "lysosome_membrane": 0.5}, "manual_coarse_organelle_mapping", "moderate"),
        "Mid piece": ({"sperm_specialized_structure": 0.5, "cilium_interior": 0.5}, "manual_specialized_mapping", "low"),
        "Mitochondria": (
            {
                "mitochondrial_matrix": 0.25,
                "mitochondrial_inner_membrane": 0.25,
                "mitochondrial_intermembrane_space": 0.25,
                "mitochondrial_outer_membrane": 0.25,
            },
            "manual_coarse_organelle_mapping",
            "moderate",
        ),
        "Perinuclear theca": ({"sperm_specialized_structure": 0.7, "nuclear_envelope": 0.3}, "manual_specialized_mapping", "low"),
        "Peroxisomes": ({"peroxisome_matrix": 0.5, "peroxisome_membrane": 0.5}, "manual_coarse_organelle_mapping", "moderate"),
        "Primary cilium": ({"cilium_interior": 0.6, "ciliary_membrane": 0.4}, "manual_coarse_organelle_mapping", "moderate"),
        "Primary cilium transition zone": ({"cilium_interior": 0.5, "ciliary_membrane": 0.5}, "manual_interface_mapping", "moderate"),
        "Principal piece": ({"sperm_specialized_structure": 0.5, "cilium_interior": 0.5}, "manual_specialized_mapping", "low"),
        "Vesicles": ({"vesicle_lumen": 0.5, "vesicle_membrane": 0.5}, "manual_coarse_organelle_mapping", "moderate"),
    }
    if location not in maps:
        return {}, "unmapped", "low"
    weights, method, confidence = maps[location]
    return normalize(weights), method, confidence


def shortest_anchor_mapping(
    go_id: str,
    terms: dict[str, GoTerm],
    alt_to_primary: dict[str, str],
    anchors: dict[str, Anchor],
) -> tuple[dict[str, float], str, str, str]:
    resolved = alt_to_primary.get(go_id, go_id)
    if resolved not in terms:
        return {}, "go_id_not_in_release", "low", ""
    queue: deque[tuple[str, int]] = deque([(resolved, 0)])
    visited: set[str] = set()
    matches: list[tuple[Anchor, int]] = []
    minimum_distance: int | None = None
    while queue:
        current, distance = queue.popleft()
        if current in visited or (
            minimum_distance is not None and distance > minimum_distance
        ):
            continue
        visited.add(current)
        if current in anchors:
            matches.append((anchors[current], distance))
            minimum_distance = distance
            continue
        for parent in terms[current].parents:
            primary = alt_to_primary.get(parent, parent)
            if primary in terms:
                queue.append((primary, distance + 1))
    if not matches:
        return {}, "no_anchor_ancestor", "low", ""
    accumulated: dict[str, float] = {}
    for anchor, _distance in matches:
        for zone, weight in anchor.zones.items():
            accumulated[zone] = accumulated.get(zone, 0.0) + weight
    zones = normalize(accumulated)
    specificities = {anchor.specificity for anchor, _ in matches}
    confidence = (
        "high"
        if specificities == {"specific"}
        else "moderate"
        if "generic" not in specificities
        else "low"
    )
    evidence = ";".join(
        f"{anchor.go_id}:{anchor.description}:distance={distance}"
        for anchor, distance in matches
    )
    return zones, "go_ancestor_anchor", confidence, evidence


def keyword_fallback(name: str) -> tuple[dict[str, float], str]:
    value = name.lower()
    patterns: list[tuple[tuple[str, ...], dict[str, float], str]] = [
        (("mitochondrial inner membrane",), {"mitochondrial_inner_membrane": 1}, "mitochondrial inner membrane"),
        (("mitochondrial outer membrane",), {"mitochondrial_outer_membrane": 1}, "mitochondrial outer membrane"),
        (("mitochondrial matrix",), {"mitochondrial_matrix": 1}, "mitochondrial matrix"),
        (("intermembrane space",), {"mitochondrial_intermembrane_space": 1}, "mitochondrial intermembrane space"),
        (("endoplasmic reticulum membrane", "er membrane"), {"er_membrane": 1}, "ER membrane"),
        (("endoplasmic reticulum lumen", "er lumen"), {"er_lumen": 1}, "ER lumen"),
        (("golgi membrane",), {"golgi_membrane": 1}, "Golgi membrane"),
        (("golgi lumen",), {"golgi_lumen": 1}, "Golgi lumen"),
        (("endosome membrane", "endosomal membrane"), {"endosome_membrane": 1}, "endosome membrane"),
        (("endosome lumen", "endosomal lumen"), {"endosome_lumen": 1}, "endosome lumen"),
        (("lysosome membrane", "lysosomal membrane"), {"lysosome_membrane": 1}, "lysosome membrane"),
        (("lysosome lumen", "lysosomal lumen"), {"lysosome_lumen": 1}, "lysosome lumen"),
        (("peroxisome membrane", "peroxisomal membrane"), {"peroxisome_membrane": 1}, "peroxisome membrane"),
        (("peroxisome matrix", "peroxisomal matrix"), {"peroxisome_matrix": 1}, "peroxisome matrix"),
        (("ciliary membrane",), {"ciliary_membrane": 1}, "ciliary membrane"),
        (("plasma membrane",), {"plasma_membrane": 1}, "plasma membrane"),
        (("nuclear envelope", "nuclear membrane"), {"nuclear_envelope": 1}, "nuclear envelope"),
        (("extracellular", "basement membrane"), {"extracellular_space": 1}, "extracellular"),
        (("cell junction", "focal adhesion"), {"cell_surface_junction": 1}, "cell-surface junction"),
        (("lipid droplet",), {"lipid_droplet_surface": 1}, "lipid droplet"),
        (
            ("chromosome", "chromatin", "nucleol", "nucleoplasm", "nuclear body", "nuclear speckle"),
            {"nucleus_interior": 1},
            "nucleus interior",
        ),
        (
            (
                "ruffle",
                "lamellipodium",
                "pseudopodium",
                "cell leading edge",
                "growth cone",
                "brush border",
                "uropod",
                "cell periphery",
                "cell surface",
                "apical part of cell",
                "basal part of cell",
                "subapical part of cell",
                "basolateral part of cell",
                "glial limiting end-foot",
            ),
            {
                "plasma_membrane": 0.4,
                "cytoskeleton": 0.35,
                "cytosol_bulk": 0.25,
            },
            "coarse cortical or cell-surface region",
        ),
        (
            (
                "axon",
                "dendrite",
                "neuron projection",
                "cell projection",
                "neuronal cell body",
                "cell body",
                "perikaryon",
                "astrocyte projection",
                "glial cell projection",
            ),
            {
                "cytoplasm_general": 0.5,
                "plasma_membrane": 0.25,
                "cytoskeleton": 0.25,
            },
            "coarse cellular projection or regional morphology",
        ),
        (
            (
                "myelin sheath",
                "internode region",
                "node of ranvier",
                "paranode region",
                "schmidt-lanterman incisure",
            ),
            {
                "generic_membrane": 0.7,
                "cytoskeleton": 0.15,
                "extracellular_space": 0.15,
            },
            "coarse myelin or axoglial membrane region",
        ),
        (
            (
                "intercellular bridge",
                "contractile ring",
                "flemming body",
                "basal body patch",
            ),
            {"cytoskeleton": 1},
            "cytoskeletal cell-division structure",
        ),
        (
            ("photoreceptor inner segment",),
            {
                "cytoplasm_general": 0.5,
                "generic_membrane": 0.25,
                "cilium_interior": 0.25,
            },
            "coarse photoreceptor cellular region",
        ),
        (
            ("sperm head-tail coupling apparatus", "lateral loop"),
            {"sperm_specialized_structure": 1},
            "sperm-specific structure",
        ),
        (
            ("organelle",),
            {"intracellular_unspecified": 1},
            "unresolved organelle",
        ),
        (
            ("phospholipid-translocating atpase complex",),
            {"generic_membrane": 1},
            "membrane-embedded ATPase complex",
        ),
        (
            ("cytoskeleton", "microtubule", "actin filament", "centrosome", "spindle", "midbody", "centriolar"),
            {"cytoskeleton": 1},
            "cytoskeleton",
        ),
        (("cytosol",), {"cytosol_bulk": 1}, "cytosol"),
        (("cytoplasmic",), {"cytoplasm_general": 1}, "cytoplasm"),
        (
            ("vesicle membrane", "autophagosome membrane", "phagosome membrane", "vacuolar membrane"),
            {"vesicle_membrane": 1},
            "generic vesicle membrane",
        ),
        (
            ("vesicle lumen", "vacuole lumen", "granule lumen"),
            {"vesicle_lumen": 1},
            "generic vesicle lumen",
        ),
        (("cilium", "ciliary"), {"cilium_interior": 0.6, "ciliary_membrane": 0.4}, "cilium"),
        (("membrane",), {"generic_membrane": 1}, "generic membrane"),
        (("intracellular",), {"intracellular_unspecified": 1}, "intracellular"),
    ]
    for needles, weights, explanation in patterns:
        if any(needle in value for needle in needles):
            return normalize(weights), explanation
    return {}, ""


def mapping_string(weights: dict[str, float]) -> str:
    return ";".join(f"{zone}:{weight:.6g}" for zone, weight in weights.items())


def project_native_matrix(
    labels: list[str],
    mappings: dict[str, dict[str, float]],
    canonical: np.ndarray,
) -> np.ndarray:
    zone_index = {zone: index for index, zone in enumerate(CANONICAL_ZONES)}
    result = np.zeros((len(labels), len(labels)), dtype=float)
    for left, label_a in enumerate(labels):
        result[left, left] = 1.0
        map_a = mappings[label_a]
        for right in range(left + 1, len(labels)):
            label_b = labels[right]
            map_b = mappings[label_b]
            score = 0.0
            same_resolved_zone = False
            for zone_a, weight_a in map_a.items():
                for zone_b, weight_b in map_b.items():
                    score += (
                        weight_a
                        * weight_b
                        * canonical[zone_index[zone_a], zone_index[zone_b]]
                    )
                    same_resolved_zone = same_resolved_zone or zone_a == zone_b
            if same_resolved_zone:
                score = min(score, 0.6)
            score = float(np.clip(score, 0.0, 1.0))
            result[left, right] = score
            result[right, left] = score
    return result


def write_matrix(path: Path, labels: list[str], matrix: np.ndarray) -> None:
    frame = pd.DataFrame(matrix, index=labels, columns=labels)
    frame.index.name = "native_location"
    frame.to_csv(path, sep="\t", float_format="%.6g")


def build_mpkccd_matrices() -> tuple[list[str], np.ndarray, np.ndarray, pd.DataFrame]:
    labels = ["1K", "4K", "17K", "200Kp", "200Ks"]
    primary = np.eye(len(labels), dtype=float)
    sensitivity = np.eye(len(labels), dtype=float)
    index = {label: i for i, label in enumerate(labels)}
    records: list[dict[str, object]] = []
    candidates = [
        ("1K", "4K", 0.35, "Neighboring sequential pellets; possible technical spillover."),
        ("4K", "17K", 0.35, "Neighboring sequential pellets; possible technical spillover."),
        ("17K", "200Kp", 0.35, "Neighboring sequential pellets; possible technical spillover."),
        ("1K", "17K", 0.15, "Second-neighbor pellets; weaker possible technical spillover."),
        ("4K", "200Kp", 0.15, "Second-neighbor pellets; weaker possible technical spillover."),
        ("17K", "200Ks", 0.10, "Pellet-to-final-supernatant boundary; highly provisional."),
        ("200Kp", "200Ks", 0.15, "Pellet and supernatant from the final separation; not physical adjacency."),
    ]
    for left, right, weight, rationale in candidates:
        i, j = index[left], index[right]
        sensitivity[i, j] = weight
        sensitivity[j, i] = weight
        records.append(
            {
                "fraction_a": left,
                "fraction_b": right,
                "primary_weight": 0.0,
                "sensitivity_weight": weight,
                "interpretation": "measurement_domain_proximity_only",
                "rationale": rationale,
                "review_status": "do_not_use_as_biological_adjacency_without_validation",
                "source_url": MPKCCD_URL,
            }
        )
    return labels, primary, sensitivity, pd.DataFrame(records)


def validate_matrix(name: str, matrix: np.ndarray) -> dict[str, object]:
    return {
        "matrix": name,
        "rows": int(matrix.shape[0]),
        "columns": int(matrix.shape[1]),
        "symmetric": bool(np.allclose(matrix, matrix.T, atol=1e-12)),
        "bounded_zero_to_one": bool(
            np.all((matrix >= -1e-12) & (matrix <= 1 + 1e-12))
        ),
        "diagonal_all_one": bool(np.allclose(np.diag(matrix), 1.0)),
        "nonzero_off_diagonal_cells": int(
            np.count_nonzero(matrix - np.diag(np.diag(matrix)))
        ),
        "minimum": float(matrix.min()),
        "maximum": float(matrix.max()),
    }


def main() -> int:
    args = parse_args()
    project = args.project_root.resolve()
    output_dir = project / "data/colocalization/compatibility_matrices"
    output_dir.mkdir(parents=True, exist_ok=True)
    obo_path = output_dir / f"go-basic-{GO_RELEASE}.obo"
    hpa_dictionary_path = (
        project
        / "data/edge_characterization/localization/hpa/v25.1/processed/"
        "hpa_location_dictionary.tsv"
    )
    compartments_dictionary_path = (
        project
        / "data/colocalization/compartments/processed/"
        "compartments_location_dictionary.tsv"
    )
    for path in [obo_path, hpa_dictionary_path, compartments_dictionary_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    terms, alt_to_primary = parse_obo(obo_path)
    anchors = make_anchors()
    canonical, rule_table = canonical_rules()

    canonical_path = output_dir / "C_canonical_physical_zones.tsv"
    write_matrix(canonical_path, CANONICAL_ZONES, canonical)
    rule_table.to_csv(
        output_dir / "canonical_compatibility_rules.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        [
            {
                "canonical_zone": zone,
                "definition": ZONE_DEFINITIONS[zone],
                "review_status": "provisional_scientist_review_required",
            }
            for zone in CANONICAL_ZONES
        ]
    ).to_csv(
        output_dir / "canonical_zone_dictionary.tsv", sep="\t", index=False
    )

    hpa_dictionary = pd.read_csv(hpa_dictionary_path, sep="\t", dtype=str).fillna("")
    hpa_rows: list[dict[str, object]] = []
    hpa_mappings: dict[str, dict[str, float]] = {}
    for row in hpa_dictionary.itertuples(index=False):
        weights, method, confidence = manual_hpa_mapping(row.location)
        hpa_mappings[row.location] = weights
        hpa_rows.append(
            {
                "native_location": row.location,
                "go_ids": row.go_ids,
                "canonical_zone_weights": mapping_string(weights),
                "mapping_method": method,
                "mapping_confidence": confidence,
                "scientist_review_priority": (
                    "routine" if confidence == "high" else "high"
                ),
                "review_status": "provisional_scientist_review_required",
                "rationale": (
                    "Specific standard location."
                    if confidence == "high"
                    else "Coarse or specialized HPA label distributed across plausible physical zones."
                ),
                "source_url": HPA_DATA_URL,
            }
        )
    hpa_mapping_frame = pd.DataFrame(hpa_rows)
    hpa_mapping_frame.to_csv(
        output_dir / "hpa_location_to_canonical_mapping.tsv",
        sep="\t",
        index=False,
    )
    hpa_labels = hpa_dictionary["location"].tolist()
    hpa_matrix = project_native_matrix(hpa_labels, hpa_mappings, canonical)
    hpa_matrix_path = output_dir / "C_hpa_locations.tsv"
    write_matrix(hpa_matrix_path, hpa_labels, hpa_matrix)

    compartments_dictionary = pd.read_csv(
        compartments_dictionary_path, sep="\t", dtype=str
    ).fillna("")
    compartments_rows: list[dict[str, object]] = []
    compartments_mappings: dict[str, dict[str, float]] = {}
    compartments_labels: list[str] = []
    for row in compartments_dictionary.itertuples(index=False):
        go_id = str(row.go_id)
        resolved = alt_to_primary.get(go_id, go_id)
        term = terms.get(resolved)
        native_name = (
            term.name
            if term is not None
            and (not row.location_name or row.location_name.startswith("GO:"))
            else str(row.location_name)
        )
        label = f"{go_id} | {native_name}"
        weights, method, confidence, evidence = shortest_anchor_mapping(
            go_id, terms, alt_to_primary, anchors
        )
        fallback = ""
        if not weights:
            weights, fallback = keyword_fallback(native_name)
            if weights:
                method = "keyword_fallback"
                confidence = "low"
                evidence = fallback
        compartments_labels.append(label)
        compartments_mappings[label] = weights
        if not weights:
            review_priority = "critical"
        elif confidence == "low":
            review_priority = "high"
        elif confidence == "moderate":
            review_priority = "medium"
        else:
            review_priority = "routine"
        compartments_rows.append(
            {
                "matrix_label": label,
                "go_id": go_id,
                "resolved_go_id": resolved,
                "location_name": native_name,
                "profiled_node_count": row.profiled_node_count,
                "maximum_score_in_universe": row.maximum_score_in_universe,
                "canonical_zone_weights": mapping_string(weights),
                "mapping_method": method,
                "mapping_confidence": confidence,
                "mapping_evidence": evidence,
                "scientist_review_priority": review_priority,
                "review_status": "provisional_scientist_review_required",
                "source_url": GO_URL,
            }
        )
    compartments_mapping_frame = pd.DataFrame(compartments_rows)
    compartments_mapping_frame.to_csv(
        output_dir / "compartments_go_to_canonical_mapping.tsv",
        sep="\t",
        index=False,
    )
    compartments_matrix = project_native_matrix(
        compartments_labels, compartments_mappings, canonical
    )
    compartments_matrix_path = output_dir / "C_compartments_go_locations.tsv"
    write_matrix(
        compartments_matrix_path, compartments_labels, compartments_matrix
    )

    (
        mpk_labels,
        mpk_primary,
        mpk_sensitivity,
        mpk_rules,
    ) = build_mpkccd_matrices()
    mpk_primary_path = output_dir / "C_mpkccd_fractions_primary.tsv"
    mpk_sensitivity_path = output_dir / "C_mpkccd_fractions_sensitivity.tsv"
    write_matrix(mpk_primary_path, mpk_labels, mpk_primary)
    write_matrix(mpk_sensitivity_path, mpk_labels, mpk_sensitivity)
    mpk_rules.to_csv(
        output_dir / "mpkccd_fraction_compatibility_rules.tsv",
        sep="\t",
        index=False,
    )

    validations = [
        validate_matrix("canonical_physical_zones", canonical),
        validate_matrix("hpa_locations", hpa_matrix),
        validate_matrix("compartments_go_locations", compartments_matrix),
        validate_matrix("mpkccd_fractions_primary", mpk_primary),
        validate_matrix("mpkccd_fractions_sensitivity", mpk_sensitivity),
    ]
    if not all(
        row["symmetric"]
        and row["bounded_zero_to_one"]
        and row["diagonal_all_one"]
        for row in validations
    ):
        raise RuntimeError(f"Matrix validation failed: {validations}")

    matrix_files = [
        canonical_path,
        hpa_matrix_path,
        compartments_matrix_path,
        mpk_primary_path,
        mpk_sensitivity_path,
    ]
    summary = {
        "purpose": "Scientist-reviewable compartment compatibility kernels for adjacency-aware colocalization.",
        "interpretation": "C[a,b] is a provisional compatibility weight for capacity to interact across native locations; it is not a probability and not direct evidence of PPI.",
        "native_exact_match_weight": 1.0,
        "different_native_labels_same_resolved_zone_cap": 0.6,
        "direct_interface_default_weight": 0.7,
        "model_status": "provisional_scientist_review_required",
        "canonical_zone_count": len(CANONICAL_ZONES),
        "hpa_location_count": len(hpa_labels),
        "compartments_go_location_count": len(compartments_labels),
        "compartments_mapping_method_counts": compartments_mapping_frame[
            "mapping_method"
        ].value_counts(dropna=False).to_dict(),
        "compartments_mapping_confidence_counts": compartments_mapping_frame[
            "mapping_confidence"
        ].value_counts(dropna=False).to_dict(),
        "compartments_unmapped_count": int(
            compartments_mapping_frame["canonical_zone_weights"].eq("").sum()
        ),
        "hpa_mapping_confidence_counts": hpa_mapping_frame[
            "mapping_confidence"
        ].value_counts(dropna=False).to_dict(),
        "mpkccd_primary_choice": "Identity matrix. Fractions are measurement bins, not named biological compartments.",
        "mpkccd_sensitivity_choice": "A separate technical-proximity matrix is supplied but must not be interpreted as physical compartment adjacency.",
        "go_release": GO_RELEASE,
        "sources": {
            "go_ontology": GO_URL,
            "hpa_subcellular": HPA_URL,
            "hpa_data": HPA_DATA_URL,
            "compartments": COMPARTMENTS_URL,
            "mpkccd_fractionation": MPKCCD_URL,
            "cell_topology": TOPOLOGY_URL,
            "nuclear_topology": NUCLEAR_URL,
            "mitochondrial_topology": MITOCHONDRIA_URL,
            "er_topology": ER_URL,
        },
        "validations": validations,
        "input_sha256": {
            path.name: sha256_file(path)
            for path in [
                obo_path,
                hpa_dictionary_path,
                compartments_dictionary_path,
            ]
        },
        "matrix_sha256": {path.name: sha256_file(path) for path in matrix_files},
    }
    (output_dir / "compatibility_matrix_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
