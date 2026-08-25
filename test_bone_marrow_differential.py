"""
Unit Test Suite for Bone Marrow Differential & Hematopathology Agent
====================================================================
Comprehensive test verification covering count math, M:E ratios, cellularity,
dysplasia, WHO 2022 diagnostic classification, edge cases, and CLI workflows.
"""

import json
import unittest
from bone_marrow_differential import (
    BoneMarrowCellCounts,
    CellularityAssessment,
    CellularityStatus,
    DysplasiaFeatures,
    DysplasiaDegree,
    IronStoreGrade,
    ClinicalCaseInput,
    BoneMarrowReport,
    BoneMarrowDifferentialAnalyzer,
    format_clinical_report,
)
import cli


class TestBoneMarrowCellCounts(unittest.TestCase):
    """Test manual cell count differential mathematical calculations."""

    def test_total_count_standard_500(self):
        counts = BoneMarrowCellCounts(
            blasts=10, promyelocytes=20, myelocytes=50, metamyelocytes=60,
            band_neutrophils=70, segmented_neutrophils=100, eosinophils=15, basophils=5,
            monocytes=10, pronormoblasts=5, basophilic_normoblasts=15,
            polychromatophilic_normoblasts=60, orthochromatophilic_normoblasts=40,
            lymphocytes=30, plasma_cells=5, megakaryocytes=2, histiocytes=2, mast_cells=1
        )
        self.assertEqual(counts.total_count(), 500)

    def test_granulocytic_and_erythroid_totals(self):
        counts = BoneMarrowCellCounts(
            blasts=10, promyelocytes=10, myelocytes=30, metamyelocytes=40,
            band_neutrophils=50, segmented_neutrophils=60, eosinophils=10, basophils=5,
            pronormoblasts=5, basophilic_normoblasts=15,
            polychromatophilic_normoblasts=40, orthochromatophilic_normoblasts=20
        )
        self.assertEqual(counts.total_granulocytic(), 215)
        self.assertEqual(counts.total_erythroid(), 80)

    def test_blast_percentage_calculation(self):
        counts = BoneMarrowCellCounts(blasts=25, segmented_neutrophils=475)
        self.assertEqual(counts.total_count(), 500)
        self.assertAlmostEqual(counts.blast_percentage(), 5.0, places=2)

    def test_non_erythroid_blast_percentage(self):
        # 20 blasts, 80 segs (total myeloid 100), 100 erythroid precursors -> total 200 cells
        # Non-erythroid total = 200 - 100 = 100. Non-erythroid blast % = 20/100 = 20%
        counts = BoneMarrowCellCounts(
            blasts=20, segmented_neutrophils=80,
            polychromatophilic_normoblasts=100
        )
        self.assertEqual(counts.total_count(), 200)
        self.assertAlmostEqual(counts.blast_percentage(), 10.0, places=2)
        self.assertAlmostEqual(counts.non_erythroid_blast_percentage(), 20.0, places=2)

    def test_me_ratio_normal(self):
        counts = BoneMarrowCellCounts(
            segmented_neutrophils=300,
            polychromatophilic_normoblasts=100
        )
        self.assertAlmostEqual(counts.myeloid_to_erythroid_ratio(), 3.0, places=2)

    def test_me_ratio_zero_erythroid(self):
        counts = BoneMarrowCellCounts(segmented_neutrophils=100)
        self.assertEqual(counts.myeloid_to_erythroid_ratio(), float("inf"))

    def test_validation_negative_counts(self):
        counts = BoneMarrowCellCounts(blasts=-5, segmented_neutrophils=100)
        issues = counts.validate()
        self.assertTrue(any("Negative count" in i for i in issues))

    def test_validation_insufficient_cells(self):
        counts = BoneMarrowCellCounts(segmented_neutrophils=50)
        issues = counts.validate()
        self.assertTrue(any("below standard minimal" in i for i in issues))


class TestCellularityAssessment(unittest.TestCase):
    """Test age-adjusted cellularity logic and reference intervals."""

    def test_expected_cellularity_adult(self):
        cell_50 = CellularityAssessment(patient_age=50, observed_cellularity_pct=50.0)
        self.assertEqual(cell_50.expected_cellularity_pct, 50.0)
        self.assertEqual(cell_50.status, CellularityStatus.NORMOCELLULAR)

    def test_expected_cellularity_elderly_normocellular(self):
        cell_80 = CellularityAssessment(patient_age=80, observed_cellularity_pct=25.0)
        self.assertEqual(cell_80.expected_cellularity_pct, 20.0)
        self.assertEqual(cell_80.status, CellularityStatus.NORMOCELLULAR)

    def test_severe_hypocellularity(self):
        cell = CellularityAssessment(patient_age=30, observed_cellularity_pct=5.0)
        self.assertEqual(cell.status, CellularityStatus.SEVERELY_HYPOCELLULAR)

    def test_marked_hypercellularity(self):
        cell = CellularityAssessment(patient_age=70, observed_cellularity_pct=95.0)
        self.assertEqual(cell.status, CellularityStatus.SEVERELY_HYPERCELLULAR)

    def test_invalid_age_and_pct(self):
        with self.assertRaises(ValueError):
            CellularityAssessment(patient_age=-5, observed_cellularity_pct=50.0)
        with self.assertRaises(ValueError):
            CellularityAssessment(patient_age=50, observed_cellularity_pct=150.0)


class TestDysplasiaAssessment(unittest.TestCase):
    """Test morphological dysplasia scoring and WHO lineage criteria."""

    def test_no_dysplasia(self):
        d = DysplasiaFeatures(erythroid_dysplasia_pct=5.0, granulocytic_dysplasia_pct=2.0)
        self.assertEqual(d.dysplastic_lineages_count(), 0)
        self.assertEqual(d.get_dysplasia_degree(), DysplasiaDegree.NONE)

    def test_single_lineage_dysplasia(self):
        d = DysplasiaFeatures(erythroid_dysplasia_pct=15.0, granulocytic_dysplasia_pct=5.0)
        self.assertEqual(d.dysplastic_lineages_count(), 1)
        self.assertEqual(d.get_dysplasia_degree(), DysplasiaDegree.SINGLE_LINEAGE)

    def test_multilineage_dysplasia(self):
        d = DysplasiaFeatures(
            erythroid_dysplasia_pct=20.0,
            granulocytic_dysplasia_pct=18.0,
            megakaryocytic_dysplasia_pct=12.0
        )
        self.assertEqual(d.dysplastic_lineages_count(), 3)
        self.assertEqual(d.get_dysplasia_degree(), DysplasiaDegree.MULTILINEAGE)

    def test_ring_sideroblasts_threshold(self):
        # 15% ring sideroblasts without SF3B1 qualifies as erythroid dysplasia
        d1 = DysplasiaFeatures(ring_sideroblasts_pct=16.0, sf3b1_mutation_detected=False)
        self.assertEqual(d1.dysplastic_lineages_count(), 1)

        # 8% ring sideroblasts WITH SF3B1 qualifies as erythroid dysplasia (>=5% threshold)
        d2 = DysplasiaFeatures(ring_sideroblasts_pct=8.0, sf3b1_mutation_detected=True)
        self.assertEqual(d2.dysplastic_lineages_count(), 1)


class TestDiagnosticClassification(unittest.TestCase):
    """Test WHO 2022 / ICC diagnostic classifications."""

    def test_aml_by_marrow_blast_percentage(self):
        counts = BoneMarrowCellCounts(blasts=110, segmented_neutrophils=390)
        case = ClinicalCaseInput(
            case_id="TEST-AML-01",
            patient_age=60,
            counts=counts,
            core_cellularity_pct=80.0
        )
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        self.assertEqual(report.primary_diagnostic_category, "Acute Myeloid Leukemia (AML)")
        self.assertTrue(any("Blasts >= 20%" in c for c in report.who_2022_criteria_matched))

    def test_aml_by_defining_genetics_low_blasts(self):
        counts = BoneMarrowCellCounts(blasts=55, segmented_neutrophils=445)  # 11% blasts
        case = ClinicalCaseInput(
            case_id="TEST-APL",
            patient_age=40,
            counts=counts,
            cytogenetics_or_mutations=["PML::RARA fusion detected by FISH"]
        )
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        self.assertIn("AML", report.primary_diagnostic_category)
        self.assertIn("PML::RARA", report.subclassification)

    def test_mds_ib1(self):
        counts = BoneMarrowCellCounts(blasts=35, segmented_neutrophils=465)  # 7.0% blasts
        case = ClinicalCaseInput(
            case_id="TEST-MDS-IB1",
            patient_age=72,
            counts=counts
        )
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        self.assertIn("Increased Blasts", report.primary_diagnostic_category)
        self.assertIn("MDS-IB1", report.subclassification)
        self.assertEqual(report.ipss_r_blast_score_category, "5% to 10% (IPSS-R Score: 2)")

    def test_mds_ib2_due_to_auer_rods(self):
        counts = BoneMarrowCellCounts(blasts=20, segmented_neutrophils=480)  # 4% blasts
        case = ClinicalCaseInput(
            case_id="TEST-MDS-AUER",
            patient_age=68,
            counts=counts,
            dysplasia=DysplasiaFeatures(auer_rods_present=True)
        )
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        self.assertIn("MDS-IB2", report.subclassification)

    def test_mds_with_del5q(self):
        counts = BoneMarrowCellCounts(blasts=10, segmented_neutrophils=490)  # 2% blasts
        case = ClinicalCaseInput(
            case_id="TEST-DEL5Q",
            patient_age=65,
            counts=counts,
            dysplasia=DysplasiaFeatures(erythroid_dysplasia_pct=15.0),
            cytogenetics_or_mutations=["del(5q)(q31q33) isolated"]
        )
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        self.assertEqual(report.primary_diagnostic_category, "Myelodysplastic Neoplasm (MDS)")
        self.assertIn("del(5q)", report.subclassification)

    def test_mds_sf3b1_ring_sideroblasts(self):
        counts = BoneMarrowCellCounts(blasts=10, segmented_neutrophils=490)
        case = ClinicalCaseInput(
            case_id="TEST-MDS-RS",
            patient_age=75,
            counts=counts,
            dysplasia=DysplasiaFeatures(ring_sideroblasts_pct=22.0, sf3b1_mutation_detected=True)
        )
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        self.assertIn("Ring Sideroblasts", report.subclassification)

    def test_multiple_myeloma_high_plasma_cells(self):
        counts = BoneMarrowCellCounts(plasma_cells=320, segmented_neutrophils=180)  # 64% plasma cells
        case = ClinicalCaseInput(
            case_id="TEST-MYELOMA",
            patient_age=67,
            counts=counts
        )
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        self.assertEqual(report.primary_diagnostic_category, "Plasma Cell Neoplasm")
        self.assertIn("Multiple Myeloma", report.subclassification)

    def test_aplastic_anemia_severe_hypocellularity(self):
        counts = BoneMarrowCellCounts(blasts=2, segmented_neutrophils=198)
        case = ClinicalCaseInput(
            case_id="TEST-APLASTIC",
            patient_age=25,
            counts=counts,
            core_cellularity_pct=5.0
        )
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        self.assertIn("Aplastic Anemia", report.subclassification)

    def test_granulocytic_hyperplasia(self):
        counts = BoneMarrowCellCounts(
            myelocytes=100, metamyelocytes=150, band_neutrophils=100, segmented_neutrophils=100,
            polychromatophilic_normoblasts=50
        )  # M:E = 450 / 50 = 9:1
        case = ClinicalCaseInput(
            case_id="TEST-GRAN-HYPER",
            patient_age=50,
            counts=counts
        )
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        self.assertEqual(report.primary_diagnostic_category, "Myeloid Hyperplasia / Shift")

    def test_erythroid_hyperplasia(self):
        counts = BoneMarrowCellCounts(
            segmented_neutrophils=100,
            pronormoblasts=20, basophilic_normoblasts=50,
            polychromatophilic_normoblasts=150, orthochromatophilic_normoblasts=180
        )  # M:E = 100 / 400 = 0.25:1
        case = ClinicalCaseInput(
            case_id="TEST-ERYTH-HYPER",
            patient_age=45,
            counts=counts
        )
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        self.assertEqual(report.primary_diagnostic_category, "Erythroid Hyperplasia")


class TestReportFormattingAndCLI(unittest.TestCase):
    """Test text rendering, JSON export, and CLI commands."""

    def test_json_and_dict_serialization(self):
        counts = BoneMarrowCellCounts(blasts=5, segmented_neutrophils=495)
        case = ClinicalCaseInput(case_id="TEST-JSON", patient_age=50, counts=counts)
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        d = report.to_dict()
        self.assertIn("marrow_blast_pct", d)
        self.assertEqual(d["case_id"], "TEST-JSON")

        json_str = report.to_json()
        parsed = json.loads(json_str)
        self.assertEqual(parsed["case_id"], "TEST-JSON")

    def test_text_report_formatting(self):
        counts = BoneMarrowCellCounts(blasts=5, segmented_neutrophils=495)
        case = ClinicalCaseInput(case_id="TEST-TEXT", patient_age=50, counts=counts)
        report = BoneMarrowDifferentialAnalyzer.analyze(case)
        txt = format_clinical_report(report)
        self.assertIn("BONE MARROW DIFFERENTIAL & HISTOPATHOLOGY REPORT", txt)
        self.assertIn("TEST-TEXT", txt)

    def test_cli_demo_execution(self):
        self.assertEqual(cli.main(["--demo", "normal"]), 0)
        self.assertEqual(cli.main(["--demo", "aml"]), 0)
        self.assertEqual(cli.main(["--demo", "mds_rs"]), 0)
        self.assertEqual(cli.main(["--demo", "aplastic"]), 0)

    def test_cli_direct_args_json(self):
        ret = cli.main([
            "--case-id", "CLI-TEST",
            "--age", "55",
            "--blasts", "150",
            "--segs", "350",
            "--json"
        ])
        self.assertEqual(ret, 0)


if __name__ == "__main__":
    unittest.main()
