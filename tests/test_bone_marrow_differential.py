import csv
import json
import tempfile
import unittest
from pathlib import Path

import cli
from bone_marrow_differential import (
    BoneMarrowCellCounts, BoneMarrowDifferentialAnalyzer, CellularityAssessment,
    CellularityStatus, ClinicalCaseInput, DysplasiaDegree, DysplasiaFeatures,
    format_clinical_report,
)


def case(counts, **kwargs):
    return ClinicalCaseInput("TEST", 50, counts, **kwargs)


class CountTests(unittest.TestCase):
    def test_math(self):
        c = BoneMarrowCellCounts(blasts=25, segmented_neutrophils=375, polychromatophilic_normoblasts=100)
        self.assertEqual(c.total_count(), 500)
        self.assertEqual(c.blast_percentage(), 5.0)
        self.assertEqual(c.non_erythroid_blast_percentage(), 6.25)
        self.assertEqual(c.percentages()["blasts"], 5.0)
        self.assertEqual(BoneMarrowDifferentialAnalyzer.calculate_me_ratio(300, 100), 3.0)

    def test_infinite_me_is_strict_json(self):
        report = BoneMarrowDifferentialAnalyzer.analyze(case(BoneMarrowCellCounts(segmented_neutrophils=500)))
        self.assertEqual(report.me_ratio, float("inf"))
        self.assertIsNone(json.loads(report.to_json())["me_ratio"])
        self.assertIn("inf:1", format_clinical_report(report))

    def test_invalid_values_fail_analysis(self):
        with self.assertRaises(ValueError):
            BoneMarrowDifferentialAnalyzer.analyze(case(BoneMarrowCellCounts(blasts=-1, segmented_neutrophils=501)))
        with self.assertRaises(ValueError):
            BoneMarrowDifferentialAnalyzer.analyze(case(
                BoneMarrowCellCounts(segmented_neutrophils=500),
                dysplasia=DysplasiaFeatures(erythroid_dysplasia_pct=101),
            ))


class MorphologyTests(unittest.TestCase):
    def test_cellularity_and_dysplasia(self):
        self.assertEqual(CellularityAssessment(50, 50).status, CellularityStatus.NORMOCELLULAR)
        with self.assertRaises(ValueError):
            CellularityAssessment(-1, 50)
        d = DysplasiaFeatures(erythroid_dysplasia_pct=12, granulocytic_dysplasia_pct=11)
        self.assertEqual(d.dysplastic_lineages_count(), 2)
        self.assertEqual(d.get_dysplasia_degree(), DysplasiaDegree.MULTILINEAGE)
        self.assertEqual(DysplasiaFeatures(ring_sideroblasts_pct=25).dysplastic_lineages_count(), 0)


class InterpretationTests(unittest.TestCase):
    def test_blast_threshold_is_flag_not_full_subtype(self):
        report = BoneMarrowDifferentialAnalyzer.analyze(case(BoneMarrowCellCounts(blasts=100, segmented_neutrophils=400)))
        self.assertEqual(report.primary_diagnostic_category, "Acute leukemia blast threshold met")
        self.assertIn("lineage assignment", report.subclassification.lower())

    def test_dysplasia_alone_is_not_mds(self):
        report = BoneMarrowDifferentialAnalyzer.analyze(case(
            BoneMarrowCellCounts(segmented_neutrophils=400, polychromatophilic_normoblasts=100),
            dysplasia=DysplasiaFeatures(erythroid_dysplasia_pct=20),
        ))
        self.assertEqual(report.primary_diagnostic_category, "Significant marrow dysplasia pattern")
        self.assertTrue(any("not sufficient" in x for x in report.advisory_recommendations))

    def test_plasmacytosis_is_not_myeloma_without_context(self):
        report = BoneMarrowDifferentialAnalyzer.analyze(case(BoneMarrowCellCounts(plasma_cells=320, segmented_neutrophils=180)))
        self.assertEqual(report.primary_diagnostic_category, "Marrow plasmacytosis")
        self.assertNotIn("multiple myeloma", report.subclassification.lower())

    def test_hypocellularity_is_not_aplastic_anemia_by_itself(self):
        report = BoneMarrowDifferentialAnalyzer.analyze(case(
            BoneMarrowCellCounts(segmented_neutrophils=450, lymphocytes=50), core_cellularity_pct=5,
        ))
        self.assertEqual(report.primary_diagnostic_category, "Severe marrow hypocellularity")
        self.assertIn("cannot be assigned", report.subclassification)

    def test_monocytosis_is_not_cmml_by_itself(self):
        report = BoneMarrowDifferentialAnalyzer.analyze(case(
            BoneMarrowCellCounts(monocytes=100, segmented_neutrophils=400),
            peripheral_blood_monocyte_abs_k_ul=1.2, persistent_pb_monocytosis_documented=True,
        ))
        self.assertNotIn("CMML", report.primary_diagnostic_category)
        self.assertTrue(any("CMML assessment" in x for x in report.advisory_recommendations))

    def test_negative_genetic_text_is_not_detected(self):
        report = BoneMarrowDifferentialAnalyzer.analyze(case(
            BoneMarrowCellCounts(blasts=20, segmented_neutrophils=480),
            cytogenetics_or_mutations=["NPM1 not detected"],
        ))
        self.assertNotEqual(report.primary_diagnostic_category, "AML-defining genetic abnormality flag")

    def test_report_states_limitations(self):
        text = format_clinical_report(BoneMarrowDifferentialAnalyzer.analyze(case(
            BoneMarrowCellCounts(segmented_neutrophils=350, polychromatophilic_normoblasts=150)
        )))
        self.assertIn("not a standalone diagnosis", text)
        self.assertIn("WHO-HAEM5", text)
        self.assertIn("ICC 2022", text)


class CLITests(unittest.TestCase):
    def test_demos(self):
        for name in ("normal", "aml", "mds_rs", "aplastic"):
            self.assertEqual(cli.main(["--demo", name]), 0)

    def test_direct_json_and_missing_counts(self):
        self.assertEqual(cli.main(["--case-id", "CLI", "--blasts", "100", "--segs", "400", "--json"]), 0)
        with self.assertRaises(SystemExit) as caught:
            cli.main([])
        self.assertEqual(caught.exception.code, 2)

    def test_batch_does_not_silently_default_bad_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp) / "in.csv", Path(tmp) / "out.csv"
            src.write_text("case_id,patient_age,blasts,segs\nBAD,nope,5,495\nFRACTION,50.5,5,495\nGOOD,50,5,495\n", encoding="utf-8")
            self.assertEqual(cli.process_batch_csv(str(src), str(dst)), 0)
            with dst.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertTrue(rows[0]["analysis_error"])
            self.assertIn("integer number of years", rows[1]["analysis_error"])
            self.assertEqual(rows[2]["analysis_error"], "")
            self.assertEqual(rows[2]["total_cells_counted"], "500")


if __name__ == "__main__":
    unittest.main()
