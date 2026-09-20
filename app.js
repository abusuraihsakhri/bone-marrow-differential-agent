"use strict";

const PYODIDE_BASE = "https://cdn.jsdelivr.net/pyodide/v314.0.7/full/";
const countFields = [
  "blasts", "promyelocytes", "myelocytes", "metamyelocytes", "band_neutrophils",
  "segmented_neutrophils", "eosinophils", "basophils", "monocytes", "pronormoblasts",
  "basophilic_normoblasts", "polychromatophilic_normoblasts", "orthochromatophilic_normoblasts",
  "lymphocytes", "plasma_cells", "megakaryocytes", "histiocytes", "mast_cells", "other_cells"
];
const dysplasiaFields = ["erythroid_dysplasia_pct", "granulocytic_dysplasia_pct", "megakaryocytic_dysplasia_pct", "ring_sideroblasts_pct"];
const triStateFields = [
  "cytopenia_documented", "persistent_pb_monocytosis_documented", "plasma_cell_clonality_documented",
  "myeloma_defining_event_documented", "aplastic_anemia_pb_criteria_documented"
];

let pyodide = null;
let lastReport = null;

const $ = (id) => document.getElementById(id);

function setRuntime(text, state = "") {
  const el = $("runtimeStatus");
  el.classList.remove("ready", "error");
  if (state) el.classList.add(state);
  el.querySelector("span:last-child").textContent = text;
}

function parseNumber(id, fallback = 0) {
  const value = $(id).value.trim();
  if (value === "") return fallback;
  const parsed = Number(value);
  if (!Number.isFinite(parsed)) throw new Error(`${id.replaceAll("_", " ")} must be a finite number.`);
  return parsed;
}

function parseInteger(id, fallback = 0) {
  const value = parseNumber(id, fallback);
  if (!Number.isInteger(value)) throw new Error(`${id.replaceAll("_", " ")} must be an integer.`);
  return value;
}

function parseTriState(id) {
  const raw = $(id).value;
  if (raw === "true") return true;
  if (raw === "false") return false;
  return null;
}

function updateTotal() {
  let total = 0;
  for (const field of countFields) {
    const value = Number($(field).value || 0);
    if (Number.isFinite(value) && value > 0) total += value;
  }
  $("cellTotal").textContent = `${total} cells`;
  $("cellTotal").classList.toggle("warning", total > 0 && total < 500);
}

function collectCase() {
  const counts = {};
  for (const field of countFields) counts[field] = parseInteger(field);
  if (Object.values(counts).some((value) => value < 0)) throw new Error("Cell counts cannot be negative.");
  if (Object.values(counts).reduce((a, b) => a + b, 0) === 0) throw new Error("Enter at least one marrow differential count.");

  const dysplasia = {};
  for (const field of dysplasiaFields) dysplasia[field] = parseNumber(field);
  dysplasia.sf3b1_mutation_detected = $("sf3b1_mutation_detected").checked;
  dysplasia.auer_rods_present = $("auer_rods_present").checked;

  const genetics = $("genetics").value.split(";").map((item) => item.trim()).filter(Boolean);
  const payload = {
    case_id: $("case_id").value.trim() || "CASE-001",
    patient_age: parseInteger("patient_age"),
    core_cellularity_pct: parseNumber("core_cellularity_pct"),
    peripheral_blood_blast_pct: parseNumber("peripheral_blood_blast_pct"),
    peripheral_blood_monocyte_abs_k_ul: parseNumber("peripheral_blood_monocyte_abs_k_ul", 0.5),
    iron_store_grade: parseInteger("iron_store_grade", 3),
    counts,
    dysplasia,
    cytogenetics_or_mutations: genetics
  };
  for (const field of triStateFields) payload[field] = parseTriState(field);
  return payload;
}

function clearList(listId) {
  const list = $(listId);
  while (list.firstChild) list.removeChild(list.firstChild);
}

function fillList(listId, sectionId, items) {
  clearList(listId);
  const section = $(sectionId);
  if (!items || items.length === 0) {
    section.hidden = true;
    return;
  }
  for (const text of items) {
    const li = document.createElement("li");
    li.textContent = text;
    $(listId).appendChild(li);
  }
  section.hidden = false;
}

function renderReport(report) {
  lastReport = report;
  $("emptyState").hidden = true;
  $("errorBox").hidden = true;
  $("resultContent").hidden = false;
  $("downloadButton").disabled = false;
  $("metricTotal").textContent = String(report.total_cells_counted);
  $("metricBlasts").textContent = `${Number(report.marrow_blast_pct).toFixed(2)}%`;
  $("metricME").textContent = report.me_ratio === null ? "∞" : `${Number(report.me_ratio).toFixed(2)}:1`;
  $("metricCellularity").textContent = report.cellularity.status;
  $("category").textContent = report.primary_diagnostic_category;
  $("detail").textContent = report.subclassification;
  $("whoResult").textContent = report.who5_interpretation;
  $("iccResult").textContent = report.icc2022_interpretation;
  fillList("flagsList", "flagsSection", report.critical_alerts);
  fillList("notesList", "notesSection", report.advisory_recommendations);

  const limitations = $("limitations");
  limitations.replaceChildren();
  const ul = document.createElement("ul");
  for (const item of report.interpretation_limitations || []) {
    const li = document.createElement("li");
    li.textContent = item;
    ul.appendChild(li);
  }
  limitations.appendChild(ul);

  const pctContainer = $("percentages");
  pctContainer.replaceChildren();
  for (const [key, value] of Object.entries(report.differential_percentages || {})) {
    if (Number(value) <= 0) continue;
    const row = document.createElement("div");
    row.className = "percent-row";
    const label = document.createElement("span");
    label.textContent = key.replaceAll("_", " ");
    const val = document.createElement("strong");
    val.textContent = `${Number(value).toFixed(2)}%`;
    row.append(label, val);
    pctContainer.appendChild(row);
  }
}

function showError(message) {
  const box = $("errorBox");
  box.textContent = message;
  box.hidden = false;
}

async function analyzeCase(payload) {
  pyodide.globals.set("case_json", JSON.stringify(payload));
  const raw = pyodide.runPython(`
import json
from bone_marrow_differential import BoneMarrowCellCounts, DysplasiaFeatures, ClinicalCaseInput, IronStoreGrade, BoneMarrowDifferentialAnalyzer
_data = json.loads(case_json)
_counts = BoneMarrowCellCounts(**_data["counts"])
_dysplasia = DysplasiaFeatures(**_data["dysplasia"])
_case = ClinicalCaseInput(
    case_id=_data["case_id"],
    patient_age=_data["patient_age"],
    counts=_counts,
    core_cellularity_pct=_data["core_cellularity_pct"],
    peripheral_blood_blast_pct=_data["peripheral_blood_blast_pct"],
    peripheral_blood_monocyte_abs_k_ul=_data["peripheral_blood_monocyte_abs_k_ul"],
    dysplasia=_dysplasia,
    iron_store_grade=IronStoreGrade(_data["iron_store_grade"]),
    cytogenetics_or_mutations=_data["cytogenetics_or_mutations"],
    cytopenia_documented=_data["cytopenia_documented"],
    persistent_pb_monocytosis_documented=_data["persistent_pb_monocytosis_documented"],
    plasma_cell_clonality_documented=_data["plasma_cell_clonality_documented"],
    myeloma_defining_event_documented=_data["myeloma_defining_event_documented"],
    aplastic_anemia_pb_criteria_documented=_data["aplastic_anemia_pb_criteria_documented"],
)
BoneMarrowDifferentialAnalyzer.analyze(_case).to_json()
  `);
  return JSON.parse(raw);
}

function applyDemo() {
  const demo = {
    blasts: 5, promyelocytes: 10, myelocytes: 40, metamyelocytes: 60, band_neutrophils: 70,
    segmented_neutrophils: 115, eosinophils: 15, basophils: 5, monocytes: 10, pronormoblasts: 5,
    basophilic_normoblasts: 15, polychromatophilic_normoblasts: 65, orthochromatophilic_normoblasts: 35,
    lymphocytes: 40, plasma_cells: 5, megakaryocytes: 2, histiocytes: 2, mast_cells: 1, other_cells: 0
  };
  $("case_id").value = "DEMO-NORM-01";
  $("patient_age").value = "45";
  $("core_cellularity_pct").value = "55";
  $("peripheral_blood_blast_pct").value = "0";
  for (const [key, value] of Object.entries(demo)) $(key).value = String(value);
  updateTotal();
}

function clearForm() {
  $("caseForm").reset();
  $("case_id").value = "CASE-001";
  $("patient_age").value = "50";
  $("core_cellularity_pct").value = "50";
  $("peripheral_blood_blast_pct").value = "0";
  $("peripheral_blood_monocyte_abs_k_ul").value = "0.5";
  $("iron_store_grade").value = "3";
  for (const field of countFields) $(field).value = "0";
  updateTotal();
  lastReport = null;
  $("resultContent").hidden = true;
  $("emptyState").hidden = false;
  $("errorBox").hidden = true;
  $("downloadButton").disabled = true;
}

function downloadJson() {
  if (!lastReport) return;
  const blob = new Blob([JSON.stringify(lastReport, null, 2)], {type: "application/json"});
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = `${(lastReport.case_id || "case").replace(/[^a-z0-9_-]+/gi, "_")}-bone-marrow-analysis.json`;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

function setTheme(theme) {
  document.documentElement.dataset.theme = theme;
  $("themeIcon").textContent = theme === "dark" ? "☀" : "☾";
  $("themeToggle").setAttribute("aria-label", theme === "dark" ? "Switch to light theme" : "Switch to dark theme");
  localStorage.setItem("bm-theme", theme);
}

function initTabs() {
  for (const button of document.querySelectorAll("[data-input-tab]")) {
    button.addEventListener("click", () => {
      const target = button.dataset.inputTab;
      for (const other of document.querySelectorAll("[data-input-tab]")) {
        const active = other.dataset.inputTab === target;
        other.classList.toggle("active", active);
        other.setAttribute("aria-selected", String(active));
      }
      for (const pane of document.querySelectorAll("[data-input-pane]")) {
        const active = pane.dataset.inputPane === target;
        pane.hidden = !active;
        pane.classList.toggle("active", active);
      }
    });
  }
}

async function initPython() {
  try {
    setRuntime("Loading Python…");
    pyodide = await loadPyodide({indexURL: PYODIDE_BASE});
    const source = await fetch("./bone_marrow_differential.py", {cache: "no-store"});
    if (!source.ok) throw new Error(`Could not load Python engine (${source.status}).`);
    pyodide.FS.writeFile("bone_marrow_differential.py", await source.text());
    pyodide.runPython("import bone_marrow_differential");
    setRuntime("Python ready", "ready");
    $("analyzeButton").disabled = false;
    if (new URLSearchParams(window.location.search).get("smoke") === "1") {
      applyDemo();
      const report = await analyzeCase(collectCase());
      renderReport(report);
      document.body.dataset.smoke = "pass";
    }
  } catch (error) {
    console.error(error);
    setRuntime("Runtime unavailable", "error");
    showError("Python runtime failed to load. Check your network connection and reload the page.");
  }
}

function init() {
  initTabs();
  const storedTheme = localStorage.getItem("bm-theme");
  const initialTheme = storedTheme || "light";
  setTheme(initialTheme);
  $("themeToggle").addEventListener("click", () => setTheme(document.documentElement.dataset.theme === "dark" ? "light" : "dark"));
  $("demoButton").addEventListener("click", applyDemo);
  $("clearButton").addEventListener("click", clearForm);
  $("downloadButton").addEventListener("click", downloadJson);
  $("countsGrid").addEventListener("input", updateTotal);
  $("caseForm").addEventListener("submit", async (event) => {
    event.preventDefault();
    if (!pyodide) return;
    const button = $("analyzeButton");
    button.disabled = true;
    button.querySelector("span:first-child").textContent = "Analyzing…";
    try {
      const report = await analyzeCase(collectCase());
      renderReport(report);
    } catch (error) {
      console.error(error);
      showError(error.message || String(error));
    } finally {
      button.disabled = false;
      button.querySelector("span:first-child").textContent = "Analyze case";
    }
  });
  updateTotal();
  initPython();
}

document.addEventListener("DOMContentLoaded", init);
