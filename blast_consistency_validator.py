#!/usr/bin/env python3
"""
Blast-count consistency validator for bone marrow specimens.

Cross-validates three independent blast measurements that must agree within
laboratory tolerance before a diagnosis is signed out:

  1. Aspirate smear differential (500-cell count)
  2. Core-biopsy/clot section estimate (area-based)
  3. Flow cytometry blast percentage (CD45-dim SSC-low gate)

Discrepancies beyond tolerance trigger specific laboratory recommendations.
Stdlib only.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


ASPIRATE_VS_BIOPSY_TOLERANCE_PCT = 5.0   # conventional sign-out tolerance
MORPHOLOGY_VS_FLOW_TOLERANCE_PCT = 10.0  # flow lags morphology in fibrotic marrows


@dataclass
class BlastMeasurements:
    aspirate_blast_pct: float
    biopsy_blast_pct: Optional[float] = None
    flow_blast_pct: Optional[float] = None
    aspirate_hemodilute: bool = False
    flow_markers: List[str] = field(default_factory=list)


@dataclass
class ConsistencyVerdict:
    pairwise: List[Dict[str, Any]]
    max_discrepancy: float
    consistent: bool
    recommendation: str
    adjudicated_blast_pct: float


def _pair(name: str, a: float, b: float, tol: float) -> Dict[str, Any]:
    disc = abs(a - b)
    return {"comparison": name, "values": [round(a, 1), round(b, 1)],
            "discrepancy": round(disc, 1), "tolerance": tol,
            "consistent": disc <= tol}


def validate_blast_consistency(m: BlastMeasurements) -> ConsistencyVerdict:
    pairwise: List[Dict[str, Any]] = []
    recs: List[str] = []

    ref = m.aspirate_blast_pct  # aspirate differential is the primary measurement

    if m.biopsy_blast_pct is not None:
        p = _pair("aspirate vs biopsy", ref, m.biopsy_blast_pct,
                  ASPIRATE_VS_BIOPSY_TOLERANCE_PCT)
        pairwise.append(p)
        if not p["consistent"]:
            if m.aspirate_blast_pct > m.biopsy_blast_pct:
                recs.append("Aspirate exceeds biopsy: suspect focal blast aggregates or "
                            "biopsy decalcification artifact; recount additional levels")
            else:
                recs.append("Biopsy exceeds aspirate: suspect aspirate hemodilution or "
                            "fibrosis-related dry tap; use immunohistochemistry (CD34) on biopsy")

    if m.flow_blast_pct is not None:
        p = _pair("aspirate vs flow", ref, m.flow_blast_pct,
                  MORPHOLOGY_VS_FLOW_TOLERANCE_PCT)
        pairwise.append(p)
        if not p["consistent"]:
            if m.aspirate_blast_pct > m.flow_blast_pct:
                recs.append("Morphology exceeds flow: blasts may be hypogranular/lysis-fragile; "
                            "review CD34/CD117 IHC and repeat flow with viability rescue")
            else:
                recs.append("Flow exceeds morphology: possible abnormal immature myeloid "
                            "population miscounted as mature cells on smear; re-review slides")

    if m.aspirate_hemodilute:
        recs.append("Aspirate flagged hemodilute: peripheral-blood dilution underestimates "
                    "marrow blasts; rely on biopsy + flow")

    max_disc = max((p["discrepancy"] for p in pairwise), default=0.0)
    consistent = all(p["consistent"] for p in pairwise)

    candidates = [v for v in (m.aspirate_blast_pct, m.biopsy_blast_pct,
                              m.flow_blast_pct) if v is not None]
    adjudicated = m.aspirate_blast_pct
    if not consistent:
        # Prefer a 2-of-3 consensus within aspirate-vs-flow tolerance;
        # fall back to the lowest measurement when no pair agrees.
        consensus = None
        for i in range(len(candidates)):
            for j in range(i + 1, len(candidates)):
                if abs(candidates[i] - candidates[j]) <= MORPHOLOGY_VS_FLOW_TOLERANCE_PCT:
                    consensus = (candidates[i] + candidates[j]) / 2.0
        adjudicated = consensus if consensus is not None else min(candidates)
        recs.insert(0, f"DISCREPANT counts (max {max_disc:.1f}%); "
                       f"adjudicated estimate {adjudicated:.1f}% pending re-review")

    recommendation = "; ".join(recs) if recs else \
        "All blast measurements concordant; proceed with WHO classification"

    return ConsistencyVerdict(pairwise=pairwise, max_discrepancy=round(max_disc, 1),
                              consistent=consistent, recommendation=recommendation,
                              adjudicated_blast_pct=round(adjudicated, 1))


def check_aml_threshold(verdict: ConsistencyVerdict,
                        mutations: Optional[List[str]] = None) -> Dict[str, Any]:
    """Apply the WHO 20% AML threshold to the adjudicated blast percentage."""
    mutations = mutations or []
    pct = verdict.adjudicated_blast_pct
    aml_by_blast = pct >= 20.0
    genetic_aml = bool({"PML::RARA", "RUNX1::RUNX1T1", "CBFB::MYH11"} & set(mutations))
    return {
        "adjudicated_blast_pct": pct,
        "meets_aml_blast_threshold": aml_by_blast,
        "genetic_aml_defining_rearrangement": genetic_aml,
        "classification": "AML" if (aml_by_blast or genetic_aml) else
                          "MDS/CMML range - apply lineage-specific criteria",
    }


if __name__ == "__main__":
    good = BlastMeasurements(aspirate_blast_pct=12.0, biopsy_blast_pct=14.0,
                             flow_blast_pct=13.5)
    bad = BlastMeasurements(aspirate_blast_pct=25.0, biopsy_blast_pct=10.0,
                            flow_blast_pct=24.0, aspirate_hemodilute=False)
    for label, m in (("Concordant case", good), ("Discrepant case", bad)):
        v = validate_blast_consistency(m)
        print(f"\n=== {label} ===")
        for p in v.pairwise:
            flag = "OK " if p["consistent"] else "FAIL"
            print(f"[{flag}] {p['comparison']}: {p['values']} diff={p['discrepancy']}% "
                  f"(tol {p['tolerance']}%)")
        print(f"verdict: consistent={v.consistent}  adjudicated={v.adjudicated_blast_pct}%")
        print(f"action : {v.recommendation}")
        print(check_aml_threshold(v))
