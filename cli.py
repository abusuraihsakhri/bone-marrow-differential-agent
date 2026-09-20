#!/usr/bin/env python3
"""CLI for the bone marrow differential calculator."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from bone_marrow_differential import (
    BoneMarrowCellCounts,
    BoneMarrowDifferentialAnalyzer,
    ClinicalCaseInput,
    DysplasiaFeatures,
    IronStoreGrade,
    format_clinical_report,
)

COUNT_ARGS = {
    "blasts": "blasts", "promyelocytes": "promyelocytes", "myelocytes": "myelocytes",
    "metamyelocytes": "metamyelocytes", "bands": "band_neutrophils", "segs": "segmented_neutrophils",
    "eosinophils": "eosinophils", "basophils": "basophils", "monocytes": "monocytes",
    "pronormo": "pronormoblasts", "baso_normo": "basophilic_normoblasts",
    "poly_normo": "polychromatophilic_normoblasts", "ortho_normo": "orthochromatophilic_normoblasts",
    "lymphocytes": "lymphocytes", "plasma_cells": "plasma_cells",
}
CSV_ALIASES = {
    "blasts": ["blasts", "blast_count"], "promyelocytes": ["promyelocytes", "promyelo"],
    "myelocytes": ["myelocytes", "myelo"], "metamyelocytes": ["metamyelocytes", "metamyelo"],
    "band_neutrophils": ["band_neutrophils", "bands"],
    "segmented_neutrophils": ["segmented_neutrophils", "segs", "neutrophils"],
    "eosinophils": ["eosinophils", "eos"], "basophils": ["basophils", "baso"],
    "monocytes": ["monocytes", "monos"], "pronormoblasts": ["pronormoblasts", "pronormo"],
    "basophilic_normoblasts": ["basophilic_normoblasts", "baso_normo"],
    "polychromatophilic_normoblasts": ["polychromatophilic_normoblasts", "poly_normo"],
    "orthochromatophilic_normoblasts": ["orthochromatophilic_normoblasts", "ortho_normo"],
    "lymphocytes": ["lymphocytes", "lymphs"], "plasma_cells": ["plasma_cells", "plasma"],
    "megakaryocytes": ["megakaryocytes", "megas"], "histiocytes": ["histiocytes"],
    "mast_cells": ["mast_cells"], "other_cells": ["other_cells", "other"],
}


def _demo_cases() -> Dict[str, ClinicalCaseInput]:
    normal = BoneMarrowCellCounts(
        blasts=5, promyelocytes=10, myelocytes=40, metamyelocytes=60, band_neutrophils=70,
        segmented_neutrophils=115, eosinophils=15, basophils=5, monocytes=10, pronormoblasts=5,
        basophilic_normoblasts=15, polychromatophilic_normoblasts=65,
        orthochromatophilic_normoblasts=35, lymphocytes=40, plasma_cells=5,
        megakaryocytes=2, histiocytes=2, mast_cells=1,
    )
    return {
        "normal": ClinicalCaseInput("DEMO-NORM-01", 45, normal, 55),
        "aml": ClinicalCaseInput(
            "DEMO-AML-01", 62,
            BoneMarrowCellCounts(
                blasts=185, promyelocytes=20, myelocytes=30, metamyelocytes=25,
                band_neutrophils=30, segmented_neutrophils=50, eosinophils=5, basophils=2,
                monocytes=15, pronormoblasts=8, basophilic_normoblasts=12,
                polychromatophilic_normoblasts=35, orthochromatophilic_normoblasts=23,
                lymphocytes=50, plasma_cells=10,
            ),
            90, 24, cytogenetics_or_mutations=["NPM1 mutated", "FLT3-ITD"],
        ),
        "mds_rs": ClinicalCaseInput(
            "DEMO-MDS-RS", 71,
            BoneMarrowCellCounts(
                blasts=12, promyelocytes=15, myelocytes=45, metamyelocytes=50,
                band_neutrophils=60, segmented_neutrophils=85, eosinophils=10, basophils=3,
                monocytes=10, pronormoblasts=10, basophilic_normoblasts=25,
                polychromatophilic_normoblasts=90, orthochromatophilic_normoblasts=45,
                lymphocytes=40, plasma_cells=10,
            ),
            60, dysplasia=DysplasiaFeatures(erythroid_dysplasia_pct=25, ring_sideroblasts_pct=28, sf3b1_mutation_detected=True),
            iron_store_grade=IronStoreGrade.GRADE_4, cytopenia_documented=True,
        ),
        "aplastic": ClinicalCaseInput(
            "DEMO-HYPOCELLULAR-01", 28,
            BoneMarrowCellCounts(
                blasts=2, promyelocytes=2, myelocytes=5, metamyelocytes=8, band_neutrophils=10,
                segmented_neutrophils=15, eosinophils=2, basophils=1, monocytes=3,
                pronormoblasts=1, basophilic_normoblasts=2, polychromatophilic_normoblasts=6,
                orthochromatophilic_normoblasts=8, lymphocytes=50, plasma_cells=5,
            ),
            8,
        ),
    }


def run_demo(name: str = "all") -> int:
    cases = _demo_cases()
    selected = cases.items() if name == "all" else [(name, cases[name])]
    for _, case in selected:
        print(format_clinical_report(BoneMarrowDifferentialAnalyzer.analyze(case)), "\n")
    return 0


def interactive_mode() -> int:
    print("Bone Marrow Differential — interactive entry")
    try:
        case_id = input("Case ID [BM-001]: ").strip() or "BM-001"
        age = int(input("Patient age [50]: ").strip() or "50")
        cellularity = float(input("Core cellularity % [50]: ").strip() or "50")
        values = {}
        for field, label in (
            ("blasts", "Blasts"), ("promyelocytes", "Promyelocytes"), ("myelocytes", "Myelocytes"),
            ("metamyelocytes", "Metamyelocytes"), ("band_neutrophils", "Bands"),
            ("segmented_neutrophils", "Segmented neutrophils"), ("eosinophils", "Eosinophils"),
            ("basophils", "Basophils"), ("monocytes", "Monocytes"), ("pronormoblasts", "Pronormoblasts"),
            ("basophilic_normoblasts", "Basophilic normoblasts"),
            ("polychromatophilic_normoblasts", "Polychromatophilic normoblasts"),
            ("orthochromatophilic_normoblasts", "Orthochromatophilic normoblasts"),
            ("lymphocytes", "Lymphocytes"), ("plasma_cells", "Plasma cells"),
        ):
            values[field] = int(input(f"  {label} [0]: ").strip() or "0")
        case = ClinicalCaseInput(case_id, age, BoneMarrowCellCounts(**values), cellularity)
        print("\n" + format_clinical_report(BoneMarrowDifferentialAnalyzer.analyze(case)))
        return 0
    except (ValueError, EOFError, KeyboardInterrupt) as exc:
        print(f"Input error: {exc}", file=sys.stderr)
        return 1


def _raw(row: Dict[str, str], keys: List[str]) -> Optional[str]:
    for key in keys:
        value = row.get(key)
        if value not in (None, ""):
            return value
    return None


def _float(row: Dict[str, str], keys: List[str], default: float) -> float:
    value = _raw(row, keys)
    return default if value is None else float(value)


def _int(row: Dict[str, str], keys: List[str], default: int = 0) -> int:
    value = _raw(row, keys)
    if value is None:
        return default
    number = float(value)
    if not number.is_integer():
        raise ValueError(f"{keys[0]} must be an integer count, got {value!r}")
    return int(number)


def _bool(row: Dict[str, str], keys: List[str]) -> Optional[bool]:
    value = _raw(row, keys)
    if value is None:
        return None
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes", "y", "positive"}:
        return True
    if normalized in {"0", "false", "no", "n", "negative"}:
        return False
    raise ValueError(f"{keys[0]} must be a boolean-like value, got {value!r}")


def process_batch_csv(input_path: str, output_path: Optional[str] = None) -> int:
    path = Path(input_path)
    if not path.is_file():
        print(f"Error: input file not found: {input_path}", file=sys.stderr)
        return 1
    results: List[Dict[str, Any]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            print("Error: CSV has no header row.", file=sys.stderr)
            return 1
        for index, row in enumerate(reader, 1):
            case_id = (row.get("case_id") or f"CASE-{index:03d}").strip()
            try:
                age_value = _float(row, ["patient_age", "age"], 50)
                if not age_value.is_integer():
                    raise ValueError(f"patient_age must be an integer number of years, got {age_value!r}")
                counts = BoneMarrowCellCounts(**{field: _int(row, aliases) for field, aliases in CSV_ALIASES.items()})
                dysplasia = DysplasiaFeatures(
                    erythroid_dysplasia_pct=_float(row, ["erythroid_dysplasia_pct", "erythroid_dysp"], 0),
                    granulocytic_dysplasia_pct=_float(row, ["granulocytic_dysplasia_pct", "granulocytic_dysp"], 0),
                    megakaryocytic_dysplasia_pct=_float(row, ["megakaryocytic_dysplasia_pct", "megakaryocytic_dysp"], 0),
                    ring_sideroblasts_pct=_float(row, ["ring_sideroblasts_pct", "ring_sideroblasts", "rs_pct"], 0),
                    sf3b1_mutation_detected=bool(_bool(row, ["sf3b1_mutated", "sf3b1"]) or False),
                    auer_rods_present=bool(_bool(row, ["auer_rods_present", "auer_rods"]) or False),
                )
                genetics = [x.strip() for x in (row.get("genetics") or "").split(";") if x.strip()]
                case = ClinicalCaseInput(
                    case_id, int(age_value), counts,
                    _float(row, ["core_cellularity_pct", "cellularity"], 50),
                    _float(row, ["peripheral_blood_blast_pct", "pb_blasts"], 0),
                    _float(row, ["peripheral_blood_monocyte_abs_k_ul", "pb_monos"], 0.5),
                    dysplasia, cytogenetics_or_mutations=genetics,
                    cytopenia_documented=_bool(row, ["cytopenia_documented"]),
                    persistent_pb_monocytosis_documented=_bool(row, ["persistent_pb_monocytosis_documented"]),
                    plasma_cell_clonality_documented=_bool(row, ["plasma_cell_clonality_documented"]),
                    myeloma_defining_event_documented=_bool(row, ["myeloma_defining_event_documented"]),
                    aplastic_anemia_pb_criteria_documented=_bool(row, ["aplastic_anemia_pb_criteria_documented"]),
                )
                report = BoneMarrowDifferentialAnalyzer.analyze(case)
                results.append({
                    "case_id": report.case_id, "patient_age": report.patient_age,
                    "total_cells_counted": report.total_cells_counted, "marrow_blast_pct": report.marrow_blast_pct,
                    "me_ratio": "inf" if report.me_ratio == float("inf") else report.me_ratio,
                    "cellularity_status": report.cellularity.status.value,
                    "dysplasia_degree": report.dysplasia_degree.value,
                    "interpretive_category": report.primary_diagnostic_category, "detail": report.subclassification,
                    "who5_interpretation": report.who5_interpretation, "icc2022_interpretation": report.icc2022_interpretation,
                    "review_flags": "; ".join(report.critical_alerts), "analysis_error": "",
                })
            except (TypeError, ValueError) as exc:
                results.append({
                    "case_id": case_id, "patient_age": "", "total_cells_counted": "", "marrow_blast_pct": "",
                    "me_ratio": "", "cellularity_status": "", "dysplasia_degree": "", "interpretive_category": "",
                    "detail": "", "who5_interpretation": "", "icc2022_interpretation": "", "review_flags": "",
                    "analysis_error": str(exc),
                })
    if not results:
        print("Warning: no data rows found.", file=sys.stderr)
        return 0
    fieldnames = list(results[0])
    if output_path:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        handle = output.open("w", encoding="utf-8", newline="")
    else:
        handle = sys.stdout
    try:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader(); writer.writerows(results)
    finally:
        if output_path:
            handle.close()
    if output_path:
        print(f"Batch processing complete: {len(results)} row(s) -> {output_path}")
    return 0


def _case_from_json(path: str) -> ClinicalCaseInput:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    counts = BoneMarrowCellCounts(**{k: v for k, v in data.get("counts", {}).items() if k in BoneMarrowCellCounts.__dataclass_fields__})
    dysplasia = DysplasiaFeatures(**{k: v for k, v in data.get("dysplasia", {}).items() if k in DysplasiaFeatures.__dataclass_fields__})
    return ClinicalCaseInput(
        data.get("case_id", "CASE-FILE"), int(data.get("patient_age", 50)), counts,
        float(data.get("core_cellularity_pct", 50)), float(data.get("peripheral_blood_blast_pct", 0)),
        float(data.get("peripheral_blood_monocyte_abs_k_ul", 0.5)), dysplasia,
        IronStoreGrade(data.get("iron_store_grade", 3)), list(data.get("cytogenetics_or_mutations", [])),
        clinical_history=str(data.get("clinical_history", "")),
        cytopenia_documented=data.get("cytopenia_documented"),
        persistent_pb_monocytosis_documented=data.get("persistent_pb_monocytosis_documented"),
        plasma_cell_clonality_documented=data.get("plasma_cell_clonality_documented"),
        myeloma_defining_event_documented=data.get("myeloma_defining_event_documented"),
        aplastic_anemia_pb_criteria_documented=data.get("aplastic_anemia_pb_criteria_documented"),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Bone marrow differential calculator and conservative interpretation aid")
    sub = parser.add_subparsers(dest="command")
    batch = sub.add_parser("batch", help="Batch-process CSV records")
    batch.add_argument("-i", "--input", required=True); batch.add_argument("-o", "--output")
    parser.add_argument("--interactive", action="store_true")
    parser.add_argument("--demo", choices=["normal", "aml", "mds_rs", "aplastic", "all"])
    parser.add_argument("--file", "-f"); parser.add_argument("--json", "-j", action="store_true")
    parser.add_argument("--case-id", default="CASE-001"); parser.add_argument("--age", type=int, default=50)
    parser.add_argument("--cellularity", type=float, default=50); parser.add_argument("--pb-blasts", type=float, default=0)
    parser.add_argument("--pb-monos", type=float, default=0.5)
    for arg in COUNT_ARGS:
        parser.add_argument("--" + arg.replace("_", "-"), dest=arg, type=int, default=None)
    parser.add_argument("--erythroid-dysp", type=float, default=0); parser.add_argument("--granulocytic-dysp", type=float, default=0)
    parser.add_argument("--megakaryocytic-dysp", type=float, default=0); parser.add_argument("--ring-sideroblasts", type=float, default=0)
    parser.add_argument("--sf3b1", action="store_true"); parser.add_argument("--auer-rods", action="store_true")
    parser.add_argument("--genetics", nargs="*", default=[]); parser.add_argument("--iron-grade", type=int, choices=range(7), default=3)
    parser.add_argument("--cytopenia-documented", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--persistent-pb-monocytosis", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--plasma-cell-clonality", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--myeloma-defining-event", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--aplastic-pb-criteria", action=argparse.BooleanOptionalAction, default=None)
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser(); args = parser.parse_args(argv)
    if args.command == "batch":
        return process_batch_csv(args.input, args.output)
    if args.interactive:
        return interactive_mode()
    if args.demo:
        return run_demo(args.demo)
    try:
        if args.file:
            case = _case_from_json(args.file)
        else:
            provided = {field: getattr(args, arg) for arg, field in COUNT_ARGS.items() if getattr(args, arg) is not None}
            if not provided:
                parser.error("No differential counts were provided. Use count arguments, --file, --interactive, or --demo.")
            case = ClinicalCaseInput(
                args.case_id, args.age, BoneMarrowCellCounts(**provided), args.cellularity, args.pb_blasts, args.pb_monos,
                DysplasiaFeatures(args.erythroid_dysp, args.granulocytic_dysp, args.megakaryocytic_dysp,
                                  args.auer_rods, args.ring_sideroblasts, args.sf3b1),
                IronStoreGrade(args.iron_grade), args.genetics,
                cytopenia_documented=args.cytopenia_documented,
                persistent_pb_monocytosis_documented=args.persistent_pb_monocytosis,
                plasma_cell_clonality_documented=args.plasma_cell_clonality,
                myeloma_defining_event_documented=args.myeloma_defining_event,
                aplastic_anemia_pb_criteria_documented=args.aplastic_pb_criteria,
            )
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        print(f"Error: {exc}", file=sys.stderr); return 1
    print(report.to_json() if args.json else format_clinical_report(report)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
