"""
Extensión: incorporación de Janssen et al. 2019 (PLOS ONE 14(6): e0218152, CC-BY 4.0).

Estudio antes-después, unicéntrico (Amphia, Países Bajos), >=70 años, cirugía electiva por
cáncer colorrectal (CRC) o aneurisma de aorta abdominal (AAA). 267 prehabilitados vs 360 controles.

Como la población NO es GI alta, no se mezcla con el análisis principal. Se hace:
  1. Re-cálculo de OR crudos desde los conteos (validación contra lo reportado)
  2. Meta-análisis EXPLORATORIO "cirugía abdominal mayor" = 12 estudios GI alto + Janssen 2019,
     para los outcomes que se solapan, usando el OR AJUSTADO de Janssen cuando existe.

Uso: python scripts/extension_janssen2019.py
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from meta_analysis import fixed_effect, random_effects, forest, LABELS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA, FIG, RES = ROOT / "data", ROOT / "figures", ROOT / "results"


def crude_or(a, n1, c, n0):
    b, d = n1 - a, n0 - c
    lor = np.log((a * d) / (b * c))
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return lor, se


def se_from_ci(lo, hi):
    return (np.log(hi) - np.log(lo)) / (2 * 1.96)


def main():
    j = pd.read_csv(DATA / "janssen2019_outcomes.csv")
    rows = []
    for _, r in j.iterrows():
        out = dict(outcome=r.outcome)
        if pd.notna(r.events_prehab):
            lor, se = crude_or(r.events_prehab, r.n_prehab, r.events_control, r.n_control)
            out.update(pct_prehab=100 * r.events_prehab / r.n_prehab,
                       pct_control=100 * r.events_control / r.n_control,
                       OR_crude_calc=np.exp(lor), lo_crude=np.exp(lor - 1.96 * se), hi_crude=np.exp(lor + 1.96 * se))
        out.update(OR_unadj_reported=r.OR_unadj_reported, OR_adj_reported=r.OR_adj_reported,
                   lo_adj=r.lo_adj, hi_adj=r.hi_adj)
        rows.append(out)
    J = pd.DataFrame(rows)
    J.to_csv(RES / "janssen2019_recalculated.csv", index=False, float_format="%.3f")

    # ---- entradas de Janssen para el pool exploratorio
    jr = j.set_index("outcome")
    add = []
    a = jr.loc["any_complication_other_than_delirium"]
    add.append(("overall_complications", np.log(a.OR_adj_reported), se_from_ci(a.lo_adj, a.hi_adj), "ajustado"))
    a = jr.loc["readmission_30d"]
    add.append(("hospital_readmission", np.log(a.OR_adj_reported), se_from_ci(a.lo_adj, a.hi_adj), "ajustado"))
    a = jr.loc["anastomotic_leakage_CRC_only"]
    lor, se = crude_or(a.events_prehab, a.n_prehab, a.events_control, a.n_control)
    add.append(("anastomotic_leakage", lor, se, "crudo (solo CRC)"))

    eff = pd.read_csv(DATA / "study_effects.csv")
    chars = pd.read_csv(DATA / "study_characteristics.csv")
    df = eff.merge(chars, on="study", how="left")
    jan = pd.DataFrame([dict(outcome=o, study="Janssen 2019", logOR=l, SE=s, year=2019,
                             design_group="Observational", estimate=e) for o, l, s, e in add])
    df["estimate"] = "crudo (del paper GI alto)"
    pooled = pd.concat([df, jan], ignore_index=True)

    res = []
    for outc, _, _, est in add:
        base = pooled[(pooled.outcome == outc) & (pooled.study != "Janssen 2019")]
        ext = pooled[pooled.outcome == outc]
        for name, sub in [("Solo GI alto", base), ("GI alto + Janssen 2019", ext)]:
            fe = fixed_effect(sub.logOR.values, sub.SE.values)
            re = random_effects(sub.logOR.values, sub.SE.values)
            res.append(dict(outcome=outc, analysis=name, k=fe["k"],
                            OR_FE=np.exp(fe["mu"]), lo_FE=np.exp(fe["lo"]), hi_FE=np.exp(fe["hi"]),
                            OR_RE=np.exp(re["mu"]), lo_RE=np.exp(re["lo"]), hi_RE=np.exp(re["hi"]), p_RE=re["p"],
                            PI_lo=np.exp(re["pi_lo"]), PI_hi=np.exp(re["pi_hi"]), I2=fe["I2"], tau2=re["tau2"],
                            janssen_estimate=est if "Janssen" in name else ""))
            if "Janssen" in name:
                forest(sub, fe, re, f"{LABELS[outc]} — GI alto + Janssen 2019 (exploratorio)",
                       FIG / f"extension_forest_{outc}.png", fe_label="Efectos fijos (IV)")
    R = pd.DataFrame(res)
    R.to_csv(RES / "extension_pooled_with_janssen2019.csv", index=False, float_format="%.3f")

    pd.set_option("display.width", 220)
    print(J.round(2).to_string(index=False))
    print()
    print(R.drop(columns=["janssen_estimate"]).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
