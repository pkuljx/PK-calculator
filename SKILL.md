---
name: pk-calculator
description: >
  Predict human pharmacokinetic (PK) parameters from preclinical animal data using established formula-based methods.
  Supports 8 methods: IVIVE (CL), Allometry (CL/Vss), FCIM (CL), Wajima (CL/Vss), Obach (Vss), and Oie-Tozer (Vss).
  Use this skill whenever the user asks to predict, estimate, or calculate human PK parameters (clearance, volume of distribution,
  half-life) from animal data — even if they don't name a specific method. Also use when the user mentions terms like
  "allometric scaling", "IVIVE", "FCIM", "interspecies extrapolation", "PK prediction", "human CL prediction",
  "Vss prediction", "fu correction", "well-stirred model", "hepatic clearance", or provides preclinical PK data
  (mouse/rat/dog/monkey CL or Vss values) and wants human estimates.
---

# PK Parameter Calculator

Predict human CL (clearance) and Vss (volume of distribution at steady state) from preclinical species data using validated pharmacokinetic formulas.

## Workflow

1. **Identify what the user wants to predict** (CL, Vss, or both) and what data they have
2. **Select the appropriate method(s)** — suggest the best method based on available data, or use the method the user requests
3. **Collect required inputs** — ask for any missing parameters
4. **Calculate** — follow the formulas in `references/formulas.md` exactly
5. **Present results** in a clear markdown table with units
6. **Optionally visualize** — if the user wants charts (log-log regression plots, bar comparisons), generate them

## Method Selection Guide

| Method | Predicts | Minimum Data Required |
|--------|----------|----------------------|
| IVIVE | CL | In vitro CLint (microsomes or hepatocytes), fu,p, cLogP or fu,inc |
| Allometry CL | CL | CL + body weight from ≥2 species |
| Allometry Vss | Vss | Vss + body weight from ≥2 species |
| FCIM | CL | CL + body weight from ≥2 species, fu for ref species + human |
| Wajima CL | CL | Rat CL, Dog CL, MW, number of H-bond acceptors |
| Wajima Vss | Vss | Rat Vss, Dog Vss |
| Obach Vss | Vss | Dog Vss, fu for dog + human |
| Oie-Tozer Vss | Vss | Rat/Monkey/Dog Vss + fu, human fu |

When the user has data from multiple species (mouse, rat, dog, monkey), suggest **Allometry** or **FCIM** for CL. When they only have rat + dog, **Wajima** is a good option. For Vss with just dog data + fu values, suggest **Obach**. For a more physiologically-based Vss prediction, suggest **Oie-Tozer**.

If the user provides enough data for multiple methods, run all applicable methods and present a comparison table.

## Input Conventions

- **CL** input is typically in **mL/min/kg** (per-kg basis). Convert as needed:
  - mL/min/kg → L/h: multiply by 0.06
  - mL/min/kg → L/h (total): multiply by 0.06 × body weight (kg)
- **Vss** input is typically in **L/kg**
- **fu** (fraction unbound) is dimensionless, range 0–1
- **Body weights** in kg. Common defaults: Mouse 0.03, Rat 0.25, Dog 8.6, Monkey 5.9, Human 60–70

## Output Format

Present results as markdown tables. Always show:
- The method name and formula used
- Key intermediate values
- The final predicted human parameter with units
- When multiple methods are used, a summary comparison table

## Formulas

All formulas, constants, and step-by-step calculation procedures are in `references/formulas.md`. Read that file before performing any calculation.

## Half-life Estimation

Once both CL and Vss are predicted, estimate human half-life:
```
t1/2 = (0.693 × Vss) / CL
```
Make sure units are consistent (e.g., Vss in L, CL in L/h → t1/2 in hours).

## Important Notes

- When doing allometric scaling, the regression is always in **log10 space** (log10(BW) vs log10(parameter))
- For Allometry CL, after getting the initial b value from logCL regression, select the correction model:
  - b ≤ 0.71: use simple allometry (no correction)
  - 0.71 < b ≤ 1.0: use MLP (maximum lifespan potential) correction
  - b > 1.0: use BrW (brain weight) correction
- FCIM does NOT apply fu correction to input CL values — it uses total CL directly
- Wajima formulas use log10 of the input values
- Oie-Tozer uses species-specific physiological volume constants — use the exact values in the formulas reference
