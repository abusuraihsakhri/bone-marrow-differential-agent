#!/usr/bin/env python3
"""
Command-Line Interface for Bone Marrow Differential & Hematopathology Agent
===========================================================================
Supports interactive case input, direct argument specification, batch JSON/CSV processing,
pre-configured demo scenarios, and structured JSON output.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from typing import List, Optional

from bone_marrow_differential import (
    BoneMarrowCellCounts,
    BoneMarrowDifferentialAnalyzer,
    ClinicalCaseInput,
    DysplasiaFeatures,
    IronStoreGrade,
    format_clinical_report,
)


def run_demo(scenario: str = "all") -> int:
    """Runs validated clinical benchmark scenarios."""
    scenarios = {
        "normal": ClinicalCaseInput(
            case_id="DEMO-NORM-01",
            patient_age=45,
            counts=BoneMarrowCellCounts(
                blasts=5, promyelocytes=10, myelocytes=40, metamyelocytes=60,
                band_neutrophils=70, segmented_neutrophils=115, eosinophils=15, basophils=5,
                monocytes=10, pronormoblasts=5, basophilic_normoblasts=15,
                polychromatophilic_normoblasts=65, orthochromatophilic_normoblasts=35,
                lymphocytes=40, plasma_cells=5, megakaryocytes=2, histiocytes=2, mast_cells=1
            ),
            core_cellularity_pct=55.0,
            peripheral_blood_blast_pct=0.0
        ),
        "aml": ClinicalCaseInput(
            case_id="DEMO-AML-01",
            patient_age=62,
            counts=BoneMarrowCellCounts(
                blasts=185, promyelocytes=20, myelocytes=30, metamyelocytes=25,
                band_neutrophils=30, segmented_neutrophils=50, eosinophils=5, basophils=2,
                monocytes=15, pronormoblasts=8, basophilic_normoblasts=12,
                polychromatophilic_normoblasts=35, orthochromatophilic_normoblasts=23,
                lymphocytes=50, plasma_cells=10
            ),
            core_cellularity_pct=90.0,
            peripheral_blood_blast_pct=24.0,
            cytogenetics_or_mutations=["NPM1 mutated", "FLT3-ITD"]
        ),
        "mds_rs": ClinicalCaseInput(
            case_id="DEMO-MDS-RS",
            patient_age=71,
            counts=BoneMarrowCellCounts(
                blasts=12, promyelocytes=15, myelocytes=45, metamyelocytes=50,
                band_neutrophils=60, segmented_neutrophils=85, eosinophils=10, basophils=3,
                monocytes=10, pronormoblasts=10, basophilic_normoblasts=25,
                polychromatophilic_normoblasts=90, orthochromatophilic_normoblasts=45,
                lymphocytes=40, plasma_cells=10
            ),
            core_cellularity_pct=60.0,
            dysplasia=DysplasiaFeatures(
                erythroid_dysplasia_pct=25.0,
                ring_sideroblasts_pct=28.0,
                sf3b1_mutation_detected=True
            ),
            iron_store_grade=IronStoreGrade.GRADE_4
        ),
        "aplastic": ClinicalCaseInput(
            case_id="DEMO-APLASTIC-01",
            patient_age=28,
            counts=BoneMarrowCellCounts(
                blasts=2, promyelocytes=2, myelocytes=5, metamyelocytes=8,
                band_neutrophils=10, segmented_neutrophils=15, eosinophils=2, basophils=1,
                monocytes=3, pronormoblasts=1, basophilic_normoblasts=2,
                polychromatophilic_normoblasts=6, orthochromatophilic_normoblasts=8,
                lymphocytes=50, plasma_cells=5
            ),
            core_cellularity_pct=8.0,
            peripheral_blood_blast_pct=0.0
        )
    }

    selected = scenarios.items() if scenario == "all" else [(scenario, scenarios[scenario])] if scenario in scenarios else []
    if not selected:
        print(f"Unknown scenario: {scenario}. Choose from: {list(scenarios.keys())} or 'all'")
        return 1

    for name, case in selected:
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        print(format_clinical_report(report))
        print("\n")
    return 0


def interactive_mode() -> int:
    """Guides the user through entering bone marrow cell differential data."""
    print("=" * 60)
    print(" Interactive Bone Marrow Differential Data Entry")
    print("=" * 60)
    try:
        case_id = input("Enter Case ID [BM-2026-001]: ").strip() or "BM-2026-001"
        age_str = input("Enter Patient Age in years [50]: ").strip() or "50"
        age = int(age_str)
        cell_str = input("Enter Biopsy Core Cellularity % [50]: ").strip() or "50"
        cellularity = float(cell_str)

        print("\nEnter Aspirate Differential Cell Counts (Recommended: 500 total cells):")
        blasts = int(input("  Blasts [0]: ").strip() or "0")
        promyelo = int(input("  Promyelocytes [0]: ").strip() or "0")
        myelo = int(input("  Myelocytes [0]: ").strip() or "0")
        metamyelo = int(input("  Metamyelocytes [0]: ").strip() or "0")
        bands = int(input("  Band Neutrophils [0]: ").strip() or "0")
        segs = int(input("  Segmented Neutrophils [0]: ").strip() or "0")
        eos = int(input("  Eosinophils [0]: ").strip() or "0")
        baso = int(input("  Basophils [0]: ").strip() or "0")
        monos = int(input("  Monocytes [0]: ").strip() or "0")
        pronormo = int(input("  Pronormoblasts [0]: ").strip() or "0")
        baso_normo = int(input("  Basophilic Normoblasts [0]: ").strip() or "0")
        poly_normo = int(input("  Polychromatophilic Normoblasts [0]: ").strip() or "0")
        ortho_normo = int(input("  Orthochromatophilic Normoblasts [0]: ").strip() or "0")
        lymphs = int(input("  Lymphocytes [0]: ").strip() or "0")
        plasma = int(input("  Plasma Cells [0]: ").strip() or "0")

        counts = BoneMarrowCellCounts(
            blasts=blasts, promyelocytes=promyelo, myelocytes=myelo,
            metamyelocytes=metamyelo, band_neutrophils=bands,
            segmented_neutrophils=segs, eosinophils=eos, basophils=baso,
            monocytes=monos, pronormoblasts=pronormo,
            basophilic_normoblasts=baso_normo,
            polychromatophilic_normoblasts=poly_normo,
            orthochromatophilic_normoblasts=ortho_normo,
            lymphocytes=lymphs, plasma_cells=plasma
        )

        pb_blast_str = input("\nPeripheral Blood Blast % [0.0]: ").strip() or "0.0"
        pb_blast = float(pb_blast_str)

        rs_str = input("Ring Sideroblasts % [0.0]: ").strip() or "0.0"
        rs_pct = float(rs_str)

        sf3b1_str = input("SF3B1 mutation detected? (y/n) [n]: ").strip().lower()
        sf3b1 = sf3b1_str in ("y", "yes", "true", "1")

        case = ClinicalCaseInput(
            case_id=case_id,
            patient_age=age,
            counts=counts,
            core_cellularity_pct=cellularity,
            peripheral_blood_blast_pct=pb_blast,
            dysplasia=DysplasiaFeatures(
                ring_sideroblasts_pct=rs_pct,
                sf3b1_mutation_detected=sf3b1
            )
        )

        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        print("\n" + format_clinical_report(report))
        return 0

    except Exception as e:
        print(f"Error during interactive processing: {e}", file=sys.stderr)
        return 1


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Bone Marrow Differential & Hematopathology Diagnostic Engine (WHO 2022 / ICC)"
    )
    parser.add_argument("--interactive", "-i", action="store_true", help="Launch interactive differential entry mode")
    parser.add_argument("--demo", choices=["normal", "aml", "mds_rs", "aplastic", "all"], help="Run benchmark demo scenario")
    parser.add_argument("--case-id", default="CASE-001", help="Clinical case identifier")
    parser.add_argument("--age", type=int, default=50, help="Patient age in years")
    parser.add_argument("--cellularity", type=float, default=50.0, help="Biopsy core cellularity percentage")
    parser.add_argument("--pb-blasts", type=float, default=0.0, help="Peripheral blood blast percentage")
    parser.add_argument("--pb-monos", type=float, default=0.5, help="Peripheral blood absolute monocyte count (k/uL)")

    # Differential counts
    parser.add_argument("--blasts", type=int, default=0, help="Blast count")
    parser.add_argument("--promyelocytes", type=int, default=0, help="Promyelocyte count")
    parser.add_argument("--myelocytes", type=int, default=0, help="Myelocyte count")
    parser.add_argument("--metamyelocytes", type=int, default=0, help="Metamyelocyte count")
    parser.add_argument("--bands", type=int, default=0, help="Band neutrophil count")
    parser.add_argument("--segs", type=int, default=0, help="Segmented neutrophil count")
    parser.add_argument("--eosinophils", type=int, default=0, help="Eosinophil count")
    parser.add_argument("--basophils", type=int, default=0, help="Basophil count")
    parser.add_argument("--monocytes", type=int, default=0, help="Monocyte count")
    parser.add_argument("--pronormo", type=int, default=0, help="Pronormoblast count")
    parser.add_argument("--baso-normo", type=int, default=0, help="Basophilic normoblast count")
    parser.add_argument("--poly-normo", type=int, default=0, help="Polychromatophilic normoblast count")
    parser.add_argument("--ortho-normo", type=int, default=0, help="Orthochromatophilic normoblast count")
    parser.add_argument("--lymphocytes", type=int, default=0, help="Lymphocyte count")
    parser.add_argument("--plasma-cells", type=int, default=0, help="Plasma cell count")

    # Dysplasia & Genetics
    parser.add_argument("--erythroid-dysp", type=float, default=0.0, help="Erythroid dysplasia percentage")
    parser.add_argument("--granulocytic-dysp", type=float, default=0.0, help="Granulocytic dysplasia percentage")
    parser.add_argument("--megakaryocytic-dysp", type=float, default=0.0, help="Megakaryocytic dysplasia percentage")
    parser.add_argument("--ring-sideroblasts", type=float, default=0.0, help="Ring sideroblasts % of erythroid cells")
    parser.add_argument("--sf3b1", action="store_true", help="SF3B1 somatic mutation detected")
    parser.add_argument("--auer-rods", action="store_true", help="Auer rods present on morphological review")
    parser.add_argument("--genetics", nargs="*", default=[], help="Cytogenetic findings / somatic mutations")
    parser.add_argument("--iron-grade", type=int, choices=[0, 1, 2, 3, 4, 5, 6], default=3, help="Prussian blue iron store grade (0-6)")

    # Output options
    parser.add_argument("--json", "-j", action="store_true", help="Output results as JSON")
    parser.add_argument("--file", "-f", help="Load case JSON file")

    args = parser.parse_args(argv)

    if args.interactive:
        return interactive_mode()

    if args.demo:
        return run_demo(args.demo)

    if args.file:
        with open(args.file, "r") as fp:
            data = json.load(fp)
        counts_data = data.get("counts", {})
        counts = BoneMarrowCellCounts(**{k: v for k, v in counts_data.items() if hasattr(BoneMarrowCellCounts, k)})
        dysp_data = data.get("dysplasia", {})
        dysp = DysplasiaFeatures(**{k: v for k, v in dysp_data.items() if hasattr(DysplasiaFeatures, k)})
        case = ClinicalCaseInput(
            case_id=data.get("case_id", "CASE-FILE"),
            patient_age=data.get("patient_age", 50),
            counts=counts,
            core_cellularity_pct=data.get("core_cellularity_pct", 50.0),
            peripheral_blood_blast_pct=data.get("peripheral_blood_blast_pct", 0.0),
            peripheral_blood_monocyte_abs_k_ul=data.get("peripheral_blood_monocyte_abs_k_ul", 0.5),
            dysplasia=dysp,
            iron_store_grade=IronStoreGrade(data.get("iron_store_grade", 3)),
            cytogenetics_or_mutations=data.get("cytogenetics_or_mutations", [])
        )
    else:
        # Default or command-line counts
        counts = BoneMarrowCellCounts(
            blasts=args.blasts,
            promyelocytes=args.promyelocytes,
            myelocytes=args.myelocytes,
            metamyelocytes=args.metamyelocytes,
            band_neutrophils=args.bands,
            segmented_neutrophils=args.segs,
            eosinophils=args.eosinophils,
            basophils=args.basophils,
            monocytes=args.monocytes,
            pronormoblasts=args.pronormo,
            basophilic_normoblasts=args.baso_normo,
            polychromatophilic_normoblasts=args.poly_normo,
            orthochromatophilic_normoblasts=args.ortho_normo,
            lymphocytes=args.lymphocytes,
            plasma_cells=args.plasma_cells
        )

        # If all counts are zero, populate with default normal distribution
        if counts.total_count() == 0:
            counts = BoneMarrowCellCounts(
                blasts=5, promyelocytes=10, myelocytes=40, metamyelocytes=60,
                band_neutrophils=70, segmented_neutrophils=115, eosinophils=15, basophils=5,
                monocytes=10, pronormoblasts=5, basophilic_normoblasts=15,
                polychromatophilic_normoblasts=65, orthochromatophilic_normoblasts=35,
                lymphocytes=40, plasma_cells=5, megakaryocytes=2, histiocytes=2
            )

        dysp = DysplasiaFeatures(
            erythroid_dysplasia_pct=args.erythroid_dysp,
            granulocytic_dysplasia_pct=args.granulocytic_dysp,
            megakaryocytic_dysplasia_pct=args.megakaryocytic_dysp,
            auer_rods_present=args.auer_rods,
            ring_sideroblasts_pct=args.ring_sideroblasts,
            sf3b1_mutation_detected=args.sf3b1
        )

        case = ClinicalCaseInput(
            case_id=args.case_id,
            patient_age=args.age,
            counts=counts,
            core_cellularity_pct=args.cellularity,
            peripheral_blood_blast_pct=args.pb_blasts,
            peripheral_blood_monocyte_abs_k_ul=args.pb_monos,
            dysplasia=dysp,
            iron_store_grade=IronStoreGrade(args.iron_grade),
            cytogenetics_or_mutations=args.genetics
        )

    report = BoneMarrowDifferentialAnalyzer.analyze(case)

    if args.json:
        print(report.to_json())
    else:
        print(format_clinical_report(report))

    return 0


if __name__ == "__main__":
    sys.exit(main())
