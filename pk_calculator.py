#!/usr/bin/env python3
"""
PK Calculator — human PK prediction from preclinical animal data.

Implements 8 methods (see references/formulas.md):
    IVIVE (CL), Allometry CL, Allometry Vss, FCIM (CL),
    Wajima CL, Wajima Vss, Obach Vss, Oie-Tozer Vss.

Usage (CLI):
    python pk_calculator.py <method> --input data.json [--output result.json]
    python pk_calculator.py list-methods

Or import as a module:
    from pk_calculator import allometry_cl, ivive, ...
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, asdict
from typing import Any


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SPECIES_BW = {  # kg
    "mouse": 0.03,
    "rat": 0.25,
    "dog": 8.6,
    "monkey": 5.9,
    "human": 60.0,
}

IVIVE_DEFAULTS = {
    #                    MPPGL, HPGL(M/g), Qh(mL/min/kg), Liver(g/kg)
    "human":  {"MPPGL": 45, "HPGL": 139, "Qh": 20.7, "Liver": 25.7},
    "mouse":  {"MPPGL": 45, "HPGL": 135, "Qh": 90.0, "Liver": 87.5},
    "rat":    {"MPPGL": 45, "HPGL": 117, "Qh": 55.2, "Liver": 40.0},
    "dog":    {"MPPGL": 45, "HPGL": 215, "Qh": 30.9, "Liver": 32.0},
    "monkey": {"MPPGL": 45, "HPGL": 120, "Qh": 44.0, "Liver": 30.0},
}

BRAIN_SCALE = {  # for BrW = BW * scale / 100
    "mouse": 1.65,
    "rat": 0.57,
    "dog": 0.78,
    "monkey": 1.65,
    "human": 2.0,
}

OIE_TOZER_VOL = {
    #             Vp,     Ve,    Vr,    R_e/i
    "rat":    {"Vp": 0.0313, "Ve": 0.265, "Vr": 0.364, "Rei": 1.4},
    "monkey": {"Vp": 0.0448, "Ve": 0.208, "Vr": 0.485, "Rei": 1.4},
    "dog":    {"Vp": 0.0515, "Ve": 0.216, "Vr": 0.450, "Rei": 1.4},
    "human":  {"Vp": 0.0436, "Ve": 0.151, "Vr": 0.380, "Rei": 1.4},
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def brw(species: str, bw: float | None = None) -> float:
    """Brain weight in kg."""
    if bw is None:
        bw = SPECIES_BW[species]
    return bw * BRAIN_SCALE[species] / 100.0


def mlp_hours(species: str, bw: float | None = None) -> float:
    """Maximum lifespan potential in hours."""
    if bw is None:
        bw = SPECIES_BW[species]
    b = brw(species, bw)
    return 185.4 * (b ** 0.636) * (bw ** -0.225) * 8760.0


def loglog_regression(bws: list[float], ys: list[float]) -> dict[str, float]:
    """Linear regression in log10(bw) vs log10(y). Returns a, b, log10a, R2."""
    n = len(bws)
    if n < 2:
        raise ValueError("Need at least 2 data points for regression.")
    xs = [math.log10(x) for x in bws]
    lys = [math.log10(y) for y in ys]
    sx, sy = sum(xs), sum(lys)
    sxy = sum(x * y for x, y in zip(xs, lys))
    sxx = sum(x * x for x in xs)
    denom = n * sxx - sx * sx
    if denom == 0:
        raise ValueError("Regression denominator is zero (identical BW?).")
    b = (n * sxy - sx * sy) / denom
    log10_a = (sy - b * sx) / n
    a = 10 ** log10_a
    y_mean = sy / n
    ss_tot = sum((y - y_mean) ** 2 for y in lys)
    ss_res = sum((y - (log10_a + b * x)) ** 2 for x, y in zip(xs, lys))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return {"a": a, "b": b, "log10_a": log10_a, "R2": r2}


def fu_inc_from_clogp(clogp: float) -> float:
    """Estimate fu,inc from cLogP/logD."""
    return 1.0 / (1.0 + 10 ** (0.56 * clogp - 1.41))


# ---------------------------------------------------------------------------
# 1. IVIVE
# ---------------------------------------------------------------------------

def ivive(
    species: str = "human",
    CLint_mxs: float | None = None,     # μL/min/mg
    CLint_hep: float | None = None,     # μL/min/million cells
    fu_p: float = 1.0,
    fu_inc: float | None = None,
    cLogP: float | None = None,
    CL_obs: float | None = None,        # mL/min/kg, optional
    MPPGL: float | None = None,
    HPGL: float | None = None,
    Qh: float | None = None,
    Liver: float | None = None,
) -> dict[str, Any]:
    """Well-stirred-model IVIVE. Returns predictions for microsomes and hepatocytes."""
    d = IVIVE_DEFAULTS[species]
    MPPGL = MPPGL if MPPGL is not None else d["MPPGL"]
    HPGL = HPGL if HPGL is not None else d["HPGL"]
    Qh = Qh if Qh is not None else d["Qh"]
    Liver = Liver if Liver is not None else d["Liver"]

    if fu_inc is None:
        if cLogP is None:
            fu_inc = 1.0
        else:
            fu_inc = fu_inc_from_clogp(cLogP)

    out: dict[str, Any] = {
        "species": species,
        "inputs": {"fu_p": fu_p, "fu_inc": fu_inc, "cLogP": cLogP,
                   "MPPGL": MPPGL, "HPGL": HPGL, "Qh": Qh, "Liver": Liver},
    }

    def _pred(clint: float, per_g_factor: float) -> dict[str, float]:
        clint_u = clint / fu_inc
        sf = per_g_factor * Liver * fu_p * 0.001
        cl_int_scaled = clint * sf
        clh = (Qh * cl_int_scaled) / (Qh + cl_int_scaled)
        cl_int_scaled_u = clint_u * sf
        clh_u = (Qh * cl_int_scaled_u) / (Qh + cl_int_scaled_u)
        r = {"CLint": clint, "CLint_u": clint_u,
             "CLh_pred": clh, "CLh_pred_u": clh_u}
        if CL_obs is not None:
            r["IVIVC"] = CL_obs / clh if clh else None
            r["IVIVC_u"] = CL_obs / clh_u if clh_u else None
        return r

    if CLint_mxs is not None:
        out["microsomes"] = _pred(CLint_mxs, MPPGL)
    if CLint_hep is not None:
        out["hepatocytes"] = _pred(CLint_hep, HPGL)
    return out


# ---------------------------------------------------------------------------
# 2. Allometry CL
# ---------------------------------------------------------------------------

def allometry_cl(
    species_data: list[dict],       # [{species, bw, CL, fu?}]
    human_weight: float = 60.0,
    fu_human: float | None = None,
    fu_correction: bool = False,
) -> dict[str, Any]:
    """Allometric scaling with auto model selection (simple / MLP / BrW)."""
    prepared = []
    for row in species_data:
        sp = row["species"].lower()
        bw = row.get("bw") or SPECIES_BW.get(sp)
        if bw is None:
            raise ValueError(f"Missing body weight for {sp}")
        cl_mlmin_kg = row["CL"]
        cl_lh = cl_mlmin_kg * 0.06 * bw   # total L/h
        cl_used = cl_lh
        if fu_correction:
            fu = row.get("fu")
            if fu is None:
                raise ValueError(f"fu required for {sp} when fu_correction is True")
            cl_used = cl_lh / fu
        b_ = brw(sp, bw)
        m_ = mlp_hours(sp, bw)
        prepared.append({
            "species": sp, "bw": bw, "CL_input": cl_mlmin_kg,
            "CL_L_h": cl_lh, "CL_used": cl_used,
            "BrW": b_, "MLP": m_,
            "CL_MLP": cl_used * m_, "CL_BrW": cl_used * b_,
        })

    bws = [p["bw"] for p in prepared]
    cls = [p["CL_used"] for p in prepared]

    initial = loglog_regression(bws, cls)
    b0 = initial["b"]

    if 0.71 < b0 <= 1.0:
        model = "MLP"
        ys = [p["CL_MLP"] for p in prepared]
    elif b0 > 1.0:
        model = "BrW"
        ys = [p["CL_BrW"] for p in prepared]
    else:
        model = "simple"
        ys = cls

    fit = loglog_regression(bws, ys)
    a, b = fit["a"], fit["b"]

    log_yh = math.log10(a) + b * math.log10(human_weight)
    y_h = 10 ** log_yh

    h_brw = brw("human", human_weight)
    h_mlp = mlp_hours("human", human_weight)

    if model == "simple":
        cl_h = y_h
    elif model == "MLP":
        cl_h = y_h / h_mlp
    else:
        cl_h = y_h / h_brw

    per_species = []
    for p, y in zip(prepared, ys):
        yhat = 10 ** (math.log10(a) + b * math.log10(p["bw"]))
        per_species.append({
            **{k: p[k] for k in ("species", "bw", "CL_input", "CL_L_h", "CL_used")},
            "Y_observed": y, "Y_predicted": yhat,
            "residual_log10": math.log10(y) - math.log10(yhat),
        })

    out = {
        "method": "Allometry CL",
        "fu_correction": fu_correction,
        "human_weight": human_weight,
        "initial_fit": initial,
        "selected_model": model,
        "final_fit": fit,
        "per_species": per_species,
        "human_BrW": h_brw,
        "human_MLP_hours": h_mlp,
        "human_CL_L_per_h": cl_h,
    }
    if fu_correction and fu_human is not None:
        out["human_CL_bound_L_per_h"] = cl_h * fu_human
    return out


# ---------------------------------------------------------------------------
# 3. Allometry Vss
# ---------------------------------------------------------------------------

def allometry_vss(
    species_data: list[dict],   # [{species, bw, Vss, fu?}]
    human_weight: float = 60.0,
    fu_correction: bool = False,
) -> dict[str, Any]:
    prepared = []
    for row in species_data:
        sp = row["species"].lower()
        bw = row.get("bw") or SPECIES_BW.get(sp)
        vss_per_kg = row["Vss"]
        if fu_correction:
            fu = row.get("fu")
            if fu is None:
                raise ValueError(f"fu required for {sp} when fu_correction is True")
            vss_per_kg = vss_per_kg / fu
        vss_total = vss_per_kg * bw
        prepared.append({"species": sp, "bw": bw,
                         "Vss_input": row["Vss"], "Vss_used": vss_total})

    bws = [p["bw"] for p in prepared]
    ys = [p["Vss_used"] for p in prepared]
    fit = loglog_regression(bws, ys)
    a, b = fit["a"], fit["b"]
    vss_h = 10 ** (math.log10(a) + b * math.log10(human_weight))

    per_species = []
    for p in prepared:
        yhat = 10 ** (math.log10(a) + b * math.log10(p["bw"]))
        per_species.append({**p, "Vss_predicted": yhat,
                            "residual_log10": math.log10(p["Vss_used"]) - math.log10(yhat)})

    return {
        "method": "Allometry Vss",
        "fu_correction": fu_correction,
        "human_weight": human_weight,
        "fit": fit,
        "per_species": per_species,
        "human_Vss_L": vss_h,
        "human_Vss_L_per_kg": vss_h / human_weight,
    }


# ---------------------------------------------------------------------------
# 4. FCIM
# ---------------------------------------------------------------------------

def fcim(
    species_data: list[dict],  # [{species, bw, CL}]  (CL in mL/min/kg)
    fu_ref: float,
    fu_human: float,
    ref_species: str = "rat",
) -> dict[str, Any]:
    prepared = []
    for row in species_data:
        sp = row["species"].lower()
        bw = row.get("bw") or SPECIES_BW.get(sp)
        cl_total = row["CL"] * bw   # mL/min
        prepared.append({"species": sp, "bw": bw,
                         "CL_input": row["CL"], "CL_total_mL_min": cl_total})

    bws = [p["bw"] for p in prepared]
    ys = [p["CL_total_mL_min"] for p in prepared]
    fit = loglog_regression(bws, ys)
    a, b = fit["a"], fit["b"]

    Rfu = fu_ref / fu_human
    cl_human_mlmin = 33.35 * ((a / Rfu) ** 0.77)
    cl_human_lh = cl_human_mlmin * 0.06

    return {
        "method": "FCIM",
        "ref_species": ref_species,
        "fu_ref": fu_ref, "fu_human": fu_human, "Rfu": Rfu,
        "fit": fit,
        "per_species": prepared,
        "human_CL_mL_per_min": cl_human_mlmin,
        "human_CL_L_per_h": cl_human_lh,
    }


# ---------------------------------------------------------------------------
# 5. Wajima CL
# ---------------------------------------------------------------------------

def wajima_cl(CL_rat: float, CL_dog: float, MW: float, Ha: float) -> dict[str, Any]:
    log_r = math.log10(CL_rat)
    log_d = math.log10(CL_dog)
    log_cl = (0.433 * log_r
              + 1.0 * log_d
              + (-0.00627) * MW
              + 0.189 * Ha
              + (-0.00111) * log_d * MW
              + 0.0000144 * MW * MW
              + (-0.0004) * MW * Ha
              + (-0.707))
    cl_mlminkg = 10 ** log_cl
    return {
        "method": "Wajima CL",
        "inputs": {"CL_rat": CL_rat, "CL_dog": CL_dog, "MW": MW, "Ha": Ha,
                   "log10_CL_rat": log_r, "log10_CL_dog": log_d},
        "log10_CL_human": log_cl,
        "human_CL_mL_min_kg": cl_mlminkg,
        "human_CL_L_h_kg": cl_mlminkg * 0.06,
    }


# ---------------------------------------------------------------------------
# 6. Wajima Vss
# ---------------------------------------------------------------------------

def wajima_vss(Vss_rat: float, Vss_dog: float, human_weight: float = 60.0) -> dict[str, Any]:
    log_r = math.log10(Vss_rat * 1000.0)  # mL/kg
    log_d = math.log10(Vss_dog * 1000.0)
    log_vh = 0.07714 * log_r * log_d + 0.5147 * log_d + 0.586
    vss_h_L = (10 ** log_vh) * human_weight / 1000.0
    return {
        "method": "Wajima Vss",
        "inputs": {"Vss_rat_L_kg": Vss_rat, "Vss_dog_L_kg": Vss_dog,
                   "log10_Vss_rat_mL_kg": log_r, "log10_Vss_dog_mL_kg": log_d},
        "log10_Vss_human_mL_kg": log_vh,
        "human_Vss_L": vss_h_L,
        "human_Vss_L_per_kg": vss_h_L / human_weight,
    }


# ---------------------------------------------------------------------------
# 7. Obach Vss
# ---------------------------------------------------------------------------

def obach_vss(Vss_dog: float, fu_dog: float, fu_human: float,
              human_weight: float = 60.0) -> dict[str, Any]:
    vss_kg = Vss_dog * (fu_human / fu_dog)
    return {
        "method": "Obach Vss",
        "inputs": {"Vss_dog_L_kg": Vss_dog, "fu_dog": fu_dog,
                   "fu_human": fu_human, "human_weight": human_weight},
        "human_Vss_L_per_kg": vss_kg,
        "human_Vss_L": vss_kg * human_weight,
    }


# ---------------------------------------------------------------------------
# 8. Oie-Tozer Vss
# ---------------------------------------------------------------------------

def _fut(vdss: float, fu: float, Vp: float, Ve: float, Vr: float, Rei: float) -> float:
    denom = vdss - Vp - fu * Ve - (1 - fu) * Rei * Vp
    return (Vr * fu) / denom


def oie_tozer_vss(
    species_data: dict,   # {"rat":{"Vss":..,"fu":..}, "monkey":{...}, "dog":{...}}
    fu_human: float,
    human_weight: float = 60.0,
) -> dict[str, Any]:
    futs = {}
    for sp in ("rat", "monkey", "dog"):
        row = species_data.get(sp)
        if row is None:
            raise ValueError(f"Oie-Tozer requires data for {sp}")
        v = OIE_TOZER_VOL[sp]
        futs[sp] = _fut(row["Vss"], row["fu"], v["Vp"], v["Ve"], v["Vr"], v["Rei"])
    fut_avg = sum(futs.values()) / 3.0

    h = OIE_TOZER_VOL["human"]
    vdss_per_kg = (
        h["Vp"] * (1 + h["Rei"])
        + fu_human * h["Vp"] * (h["Ve"] / h["Vp"] - h["Rei"])
        + h["Vr"] * (fu_human / fut_avg)
    )
    return {
        "method": "Oie-Tozer Vss",
        "fut_per_species": futs,
        "fut_preclinical_avg": fut_avg,
        "fu_human": fu_human,
        "human_Vss_L_per_kg": vdss_per_kg,
        "human_Vss_L": vdss_per_kg * human_weight,
    }


# ---------------------------------------------------------------------------
# Half-life
# ---------------------------------------------------------------------------

def half_life(CL_L_per_h: float, Vss_L: float) -> float:
    """Terminal half-life in hours: t1/2 = 0.693 * Vss / CL."""
    return 0.693 * Vss_L / CL_L_per_h


# ---------------------------------------------------------------------------
# Dispatch / CLI
# ---------------------------------------------------------------------------

METHODS = {
    "ivive": ivive,
    "allometry-cl": allometry_cl,
    "allometry-vss": allometry_vss,
    "fcim": fcim,
    "wajima-cl": wajima_cl,
    "wajima-vss": wajima_vss,
    "obach-vss": obach_vss,
    "oie-tozer-vss": oie_tozer_vss,
    "half-life": lambda **kw: {"t_half_hours": half_life(**kw)},
}


def run(method: str, payload: dict) -> dict:
    if method not in METHODS:
        raise SystemExit(f"Unknown method '{method}'. Options: {', '.join(METHODS)}")
    return METHODS[method](**payload)


def _cli() -> None:
    p = argparse.ArgumentParser(description="PK Calculator — human PK prediction.")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("list-methods", help="List available methods.")

    run_p = sub.add_parser("run", help="Run a method with JSON input.")
    run_p.add_argument("method", choices=list(METHODS))
    run_p.add_argument("--input", "-i", required=True,
                       help="JSON file with keyword args for the method. Use '-' for stdin.")
    run_p.add_argument("--output", "-o", default="-",
                       help="Output JSON file (default: stdout).")

    args = p.parse_args()

    if args.cmd == "list-methods":
        for m in METHODS:
            print(m)
        return

    if args.input == "-":
        payload = json.load(sys.stdin)
    else:
        with open(args.input) as fh:
            payload = json.load(fh)

    result = run(args.method, payload)
    text = json.dumps(result, indent=2, default=float)
    if args.output == "-":
        print(text)
    else:
        with open(args.output, "w") as fh:
            fh.write(text)


if __name__ == "__main__":
    _cli()
