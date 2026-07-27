# PK Calculator — Complete Formula Reference

## Table of Contents
1. [Constants and Defaults](#constants-and-defaults)
2. [IVIVE (CL)](#1-ivive-cl)
3. [Allometry CL](#2-allometry-cl)
4. [Allometry Vss](#3-allometry-vss)
5. [FCIM (CL)](#4-fcim-cl)
6. [Wajima CL](#5-wajima-cl)
7. [Wajima Vss](#6-wajima-vss)
8. [Obach Vss](#7-obach-vss)
9. [Oie-Tozer Vss](#8-oie-tozer-vss)

---

## Constants and Defaults

### Species Body Weights (default, kg)
| Species | Weight |
|---------|--------|
| Mouse   | 0.03   |
| Rat     | 0.25   |
| Dog     | 8.6    |
| Monkey  | 5.9    |
| Human   | 60     |

### IVIVE Species Defaults
| Parameter | Human | Mouse | Rat  | Dog  | Cynomolgus Monkey |
|-----------|-------|-------|------|------|-------------------|
| MPPGL (mg/g) | 45 | 45 | 45 | 45 | 45 |
| HPGL (million/g) | 139 | 135 | 117 | 215 | 120 |
| Qh (mL/min/kg) | 20.7 | 90 | 55.2 | 30.9 | 44 |
| Liver weight (g/kg) | 25.7 | 87.5 | 40 | 32 | 30 |

### Brain Weight Scale Factors (for Allometry CL correction)
| Species | Scale |
|---------|-------|
| Mouse   | 1.65  |
| Rat     | 0.57  |
| Dog     | 0.78  |
| Monkey  | 1.65  |
| Human   | 2.0   |

Brain weight: `BrW = body_weight × scale / 100` (in kg)

### Maximum Lifespan Potential (MLP)
```
MLP = 185.4 × BrW^0.636 × body_weight^(-0.225) × 8760   (unit: hours)
```

### Oie-Tozer Physiological Volumes (L/kg)
| Species | Vp     | Ve    | Vr    | R_e/i |
|---------|--------|-------|-------|-------|
| Rat     | 0.0313 | 0.265 | 0.364 | 1.4   |
| Monkey  | 0.0448 | 0.208 | 0.485 | 1.4   |
| Dog     | 0.0515 | 0.216 | 0.450 | 1.4   |
| Human   | 0.0436 | 0.151 | 0.380 | 1.4   |

---

## 1. IVIVE (CL)

Predicts hepatic clearance from in vitro intrinsic clearance using the well-stirred model.

### Inputs
- `CLint_mxs`: intrinsic clearance from microsomes (μL/min/mg protein)
- `CLint_hep`: intrinsic clearance from hepatocytes (μL/min/million cells)
- `fu,p`: fraction unbound in plasma
- `fu,inc`: fraction unbound in incubation (from experiment or calculated from cLogP/logD)
- `MPPGL`, `HPGL`, `Qh`, `Liver`: species-specific constants (see table above)
- `CL,obs`: observed in vivo CL (optional, for IVIVC calculation)

### fu,inc from cLogP/logD
```
fu,inc = 1 / (1 + 10^(0.56 × cLogP - 1.41))
```

### Unbound intrinsic clearance
```
CLint,u = CLint / fu,inc
```

### Hepatic clearance prediction — Microsomes pathway
```
scaling_factor = MPPGL × Liver × fu,p × 0.001
CL_int_scaled = CLint × scaling_factor
CLh_pred = (Qh × CL_int_scaled) / (Qh + CL_int_scaled)
```
For the unbound pathway, use `CLint,u` instead of `CLint`.

### Hepatic clearance prediction — Hepatocytes pathway
```
scaling_factor = HPGL × Liver × fu,p × 0.001
CL_int_scaled = CLint × scaling_factor
CLh_pred = (Qh × CL_int_scaled) / (Qh + CL_int_scaled)
```

### IVIVC (in vitro-in vivo correlation)
```
IVIVC = CL,obs / CLh_pred
```

### Output
Report for both microsomes and hepatocytes:
- CLint,u
- CLh_pred (from CLint)
- CLh_pred,u (from CLint,u)
- IVIVC and IVIVC,u (if CL,obs provided)

Units: CLh in mL/min/kg.

---

## 2. Allometry CL

Predicts human CL from multi-species data using allometric scaling with automatic model selection.

### Inputs
- Per species: body weight (kg), CL (mL/min/kg), optionally fu
- Human body weight (default 60 kg)
- Whether to use fu correction (CL(u) = CL/fu mode)

### Step 1: Unit conversion
For each species:
```
CL_L_per_h = CL_mL_min_kg × 0.06 × body_weight    (total CL in L/h)
```
If fu correction is enabled:
```
CL_for_calc = CL_L_per_h / fu
```
Otherwise use `CL_L_per_h` directly.

### Step 2: Compute correction metrics
For each species:
```
BrW = body_weight × scale / 100
MLP = 185.4 × BrW^0.636 × body_weight^(-0.225) × 8760
CL_MLP = CL × MLP
CL_BrW = CL × BrW
```

### Step 3: Log-log linear regression
Perform linear regression in log10 space:
```
log10(Y) = log10(a) + b × log10(body_weight)
```

First fit `log10(CL)` vs `log10(BW)` to get initial b value.

### Step 4: Model selection based on b
- If 0.55 < b ≤ 0.71: use **simple allometry** (Y = CL)
- If 0.71 < b ≤ 1.0: use **MLP correction** (Y = CL × MLP)
- If b > 1.0: use **BrW correction** (Y = CL × BrW)
- If b ≤ 0.55: use **simple allometry** (default fallback)

Re-fit the selected model's Y variable against body weight.

### Step 5: Human prediction
```
log10(Y_human) = log10(a) + b × log10(human_weight)
Y_human = 10^(log10(Y_human))
```

Then extract CL:
- Simple allometry: `human_CL = Y_human`
- MLP correction: `human_CL = Y_human / human_MLP`
- BrW correction: `human_CL = Y_human / human_BrW`

Where human BrW and MLP use human scale = 2.0.

### Step 6: R² calculation
```
R² = 1 - (SS_residual / SS_total)
```

### Output
- Model parameters: a, b, R² for both initial and selected model
- Per-species: observed vs predicted values, residuals
- Human predicted CL in L/h
- If fu correction enabled: also report CL = CL(u) × fu_human

---

## 3. Allometry Vss

Predicts human Vss using simple allometric scaling.

### Inputs
- Per species: body weight (kg), Vss (L/kg), optionally fu
- Human body weight (default 60 kg)

### Calculation
For each species:
```
Vss_total = Vss_per_kg × body_weight    (L, total)
```
If fu correction: `Vss_total = (Vss_per_kg / fu) × body_weight`

Log-log regression:
```
log10(Vss_total) = log10(a) + b × log10(body_weight)
```

Human prediction:
```
log10(Vss_human) = log10(a) + b × log10(human_weight)
Vss_human = 10^(log10(Vss_human))    (in L)
```

### Output
- a, b, R²
- Per-species observed vs predicted
- Human Vss in L

---

## 4. FCIM (CL)

Fixed Cutoff Intercept Method — uses allometric intercept and fu ratio.

### Inputs
- Per species: body weight (kg), CL (mL/min/kg)
- Reference species fu and human fu
- Reference species selection (default: rat)

### Step 1: Convert CL to total (mL/min)
```
CL_total = CL_mL_min_kg × body_weight    (mL/min)
```
**FCIM does NOT apply fu correction to CL.**

### Step 2: Log-log regression
```
log10(CL_total) = log10(a) + b × log10(body_weight)
```

### Step 3: Human CL prediction
```
Rfu = fu_ref_species / fu_human
human_CL = 33.35 × (a / Rfu)^0.77    (mL/min)
```

Convert: `human_CL_L_per_h = human_CL × 0.06`

### Output
- a, b, R²
- Rfu value
- Human CL in mL/min and L/h

---

## 5. Wajima CL

Multi-parameter regression for human CL prediction.

### Inputs
- `CL_rat`: rat clearance (mL/min/kg)
- `CL_dog`: dog clearance (mL/min/kg)
- `MW`: molecular weight
- `Ha`: number of hydrogen bond acceptors

### Formula
```
log10(CL_human) = 0.433 × log10(CL_rat)
                + 1.0 × log10(CL_dog)
                + (-0.00627) × MW
                + 0.189 × Ha
                + (-0.00111) × log10(CL_dog) × MW
                + 0.0000144 × MW²
                + (-0.0004) × MW × Ha
                + (-0.707)

CL_human = 10^(log10(CL_human))    (mL/min/kg)
```

Convert: `CL_L_h_kg = CL_human × 0.06` (L/h/kg)

### Output
- Input parameters and their log10 values
- log10(CL_pred) and CL_pred in mL/min/kg
- CL_pred in L/h/kg

---

## 6. Wajima Vss

Regression for human Vss using rat and dog data.

### Inputs
- `Vss_rat`: rat Vss (L/kg)
- `Vss_dog`: dog Vss (L/kg)
- `human_weight`: human body weight (kg, default 60)

### Formula
Note: Vss values are converted to mL/kg (×1000) before taking log10.
```
log10(Vss_rat_mL) = log10(Vss_rat × 1000)
log10(Vss_dog_mL) = log10(Vss_dog × 1000)

log10(Vss_human) = 0.07714 × log10(Vss_rat_mL) × log10(Vss_dog_mL)
                 + 0.5147 × log10(Vss_dog_mL)
                 + 0.586

Vss_human = 10^(log10(Vss_human)) × human_weight / 1000    (L)
```

### Output
- Input Vss values and their log10 (in mL/kg)
- log10(Vss_human)
- Vss_human in L

---

## 7. Obach Vss

Simple fu-ratio scaling from dog to human.

### Inputs
- `Vss_dog`: dog Vss (L/kg)
- `fu_dog`: dog fraction unbound
- `fu_human`: human fraction unbound
- `human_weight`: human body weight (kg, default 60)

### Formula
```
Vss_human = Vss_dog × (fu_human / fu_dog) × human_weight    (L)
```

### Output
- Input parameters
- Vss_human in L

---

## 8. Oie-Tozer Vss

Physiologically-based Vss prediction using tissue binding (fut) from preclinical species.

### Inputs
- Rat, Monkey, Dog: Vss (L/kg) and fu
- Human: fu
- Human body weight (kg, default 60)

### Step 1: Calculate fut for each preclinical species
```
fut = (Vr × fu) / (Vdss - Vp - (fu × Ve) - ((1 - fu) × R_e/i × Vp))
```
Where Vp, Ve, Vr, R_e/i are species-specific physiological volumes (see constants table).

Use these species-specific values:
- **Rat**: Vp=0.0313, Ve=0.265, Vr=0.364, R_e/i=1.4
- **Monkey**: Vp=0.0448, Ve=0.208, Vr=0.485, R_e/i=1.4
- **Dog**: Vp=0.0515, Ve=0.216, Vr=0.450, R_e/i=1.4

### Step 2: Average preclinical fut
```
fut_preclinical = (fut_rat + fut_monkey + fut_dog) / 3
```

### Step 3: Predict human Vdss
Using human physiological volumes: Vp=0.0436, Ve=0.151, Vr=0.380, R_e/i=1.4
```
Vdss_human = Vp × (1 + R_e/i)
           + fu_human × Vp × (Ve/Vp - R_e/i)
           + Vr × (fu_human / fut_preclinical)
```

Final: `Vdss_human_total = Vdss_human × human_weight` (L)

### Output
- fut for each species
- fut_preclinical (average)
- Vdss_human per kg (L/kg)
- Vdss_human total (L)

---

## Linear Regression in Log10 Space

Used by Allometry CL, Allometry Vss, and FCIM. Here's the procedure:

Given n data points (x_i, y_i) where x = log10(BW) and y = log10(parameter):

```
b = (n × Σ(x_i × y_i) - Σx_i × Σy_i) / (n × Σ(x_i²) - (Σx_i)²)
log10(a) = (Σy_i - b × Σx_i) / n
a = 10^(log10(a))

R² = 1 - SS_res / SS_tot
SS_tot = Σ(y_i - ȳ)²
SS_res = Σ(y_i - ŷ_i)²
ŷ_i = log10(a) + b × x_i
```
