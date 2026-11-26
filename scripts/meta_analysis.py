"""
Re-análisis reproducible: prehabilitación vs. no prehabilitación en cirugía de cáncer GI alto.

Fuente de datos: estimaciones por estudio (log OR y SE) publicadas en abierto (CC-BY 4.0) en
PLOS ONE 20(3): e0315734 (2025), Figs 2-4 y Tabla 1.

Qué agrega este re-análisis respecto del paper original (que usó solo efectos fijos):
  1. Efectos aleatorios (REML) con ajuste Hartung-Knapp-Sidik-Jonkman (HKSJ)
  2. Intervalos de predicción (PI 95%)
  3. Sensibilidad: solo ECA, solo observacionales, leave-one-out
  4. Subgrupos (diseño, tipo de cáncer, tratamiento, modalidad) con test Q entre grupos
  5. Test de Egger (solo k >= 8, interpretación cautelosa)

Uso:  python scripts/meta_analysis.py   (desde la raíz del repo)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats, optimize
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory

ROOT = Path(__file__).resolve().parents[1]
DATA, FIG, RES = ROOT / "data", ROOT / "figures", ROOT / "results"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

LABELS = {
    "overall_complications": "Complicaciones postoperatorias globales",
    "pulmonary_complications": "Complicaciones pulmonares",
    "anastomotic_leakage": "Filtración anastomótica",
    "cardiac_complications": "Complicaciones cardíacas",
    "wound_infection": "Infección de herida",
    "vocal_cord_palsy": "Parálisis de cuerda vocal",
    "chylothorax": "Quilotórax",
    "postoperative_bleeding": "Sangrado postoperatorio",
    "hospital_readmission": "Reingreso hospitalario",
    "icu_readmission": "Reingreso a UCI",
}


# ----------------------------------------------------------------------------- core
def fixed_effect(y, se):
    w = 1 / se**2
    mu = np.sum(w * y) / np.sum(w)
    s = np.sqrt(1 / np.sum(w))
    q = np.sum(w * (y - mu) ** 2)
    k = len(y)
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if q > 0 else 0.0
    p = 2 * stats.norm.sf(abs(mu / s))
    return dict(mu=mu, se=s, lo=mu - 1.96 * s, hi=mu + 1.96 * s, p=p, Q=q, I2=i2, k=k)


def tau2_reml(y, se):
    v = se**2
    if len(y) < 2:
        return 0.0

    def negll(t2):
        w = 1 / (v + t2)
        mu = np.sum(w * y) / np.sum(w)
        return 0.5 * (np.sum(np.log(v + t2)) + np.log(np.sum(w)) + np.sum(w * (y - mu) ** 2))

    res = optimize.minimize_scalar(negll, bounds=(0, 50), method="bounded")
    return float(res.x) if res.x > 1e-8 else 0.0


def random_effects(y, se):
    k = len(y)
    t2 = tau2_reml(y, se)
    w = 1 / (se**2 + t2)
    mu = np.sum(w * y) / np.sum(w)
    se_z = np.sqrt(1 / np.sum(w))
    out = dict(mu=mu, tau2=t2, k=k)
    if k >= 2:
        # HKSJ con piso (max(1, q*)) -> nunca más estrecho que el IC estándar
        qstar = np.sum(w * (y - mu) ** 2) / (k - 1)
        se_h = se_z * np.sqrt(max(1.0, qstar))
        tcrit = stats.t.ppf(0.975, k - 1)
        out.update(se=se_h, lo=mu - tcrit * se_h, hi=mu + tcrit * se_h,
                   p=2 * stats.t.sf(abs(mu / se_h), k - 1))
    else:
        out.update(se=se_z, lo=mu - 1.96 * se_z, hi=mu + 1.96 * se_z,
                   p=2 * stats.norm.sf(abs(mu / se_z)))
    if k >= 3:
        tp = stats.t.ppf(0.975, k - 2)
        half = tp * np.sqrt(t2 + se_z**2)
        out.update(pi_lo=mu - half, pi_hi=mu + half)
    else:
        out.update(pi_lo=np.nan, pi_hi=np.nan)
    return out


def egger(y, se):
    z = y / se
    prec = 1 / se
    X = np.column_stack([np.ones_like(prec), prec])
    beta, *_ = np.linalg.lstsq(X, z, rcond=None)
    resid = z - X @ beta
    df = len(y) - 2
    s2 = resid @ resid / df
    cov = s2 * np.linalg.inv(X.T @ X)
    t = beta[0] / np.sqrt(cov[0, 0])
    return beta[0], 2 * stats.t.sf(abs(t), df)


def q_between(df, group_col):
    """Test de diferencias entre subgrupos (efectos fijos dentro de grupo)."""
    rows, qw = [], 0.0
    fe_all = fixed_effect(df.logOR.values, df.SE.values)
    for g, sub in df.groupby(group_col):
        r = fixed_effect(sub.logOR.values, sub.SE.values)
        rr = random_effects(sub.logOR.values, sub.SE.values)
        qw += r["Q"]
        rows.append(dict(group=g, k=r["k"], OR_FE=np.exp(r["mu"]), lo_FE=np.exp(r["lo"]),
                         hi_FE=np.exp(r["hi"]), OR_RE=np.exp(rr["mu"]), lo_RE=np.exp(rr["lo"]),
                         hi_RE=np.exp(rr["hi"])))
    qb = fe_all["Q"] - qw
    dfb = df[group_col].nunique() - 1
    p = stats.chi2.sf(qb, dfb) if dfb > 0 else np.nan
    return pd.DataFrame(rows), qb, dfb, p


# ----------------------------------------------------------------------------- plots
def forest(sub, fe, re, title, path):
    sub = sub.sort_values("year")
    k = len(sub)
    fig, ax = plt.subplots(figsize=(10, 0.42 * (k + 5) + 1.2))
    ypos = np.arange(k, 0, -1) + 2
    tl = blended_transform_factory(ax.transAxes, ax.transData)
    for yy, (_, r) in zip(ypos, sub.iterrows()):
        orr, lo, hi = np.exp(r.logOR), np.exp(r.logOR - 1.96 * r.SE), np.exp(r.logOR + 1.96 * r.SE)
        col = "#5B2C83" if r.design_group == "RCT" else "#555555"
        ax.plot([lo, hi], [yy, yy], color=col, lw=1.2)
        ax.plot(orr, yy, "s", color=col, ms=4 + 8 * (1 / r.SE**2) / (1 / sub.SE**2).max())
        lab = f"{r.study}{' (ECA)' if r.design_group == 'RCT' else ''}"
        ax.text(-0.03, yy, lab, va="center", ha="right", fontsize=9, transform=tl)
        ax.text(1.03, yy, f"{orr:.2f} [{lo:.2f}, {hi:.2f}]", va="center", fontsize=8.5, transform=tl)

    def diamond(res, yy, color, label):
        m, l, h = np.exp(res["mu"]), np.exp(res["lo"]), np.exp(res["hi"])
        ax.fill([l, m, h, m], [yy, yy + 0.3, yy, yy - 0.3], color=color)
        ax.text(-0.03, yy, label, va="center", ha="right", fontsize=9, weight="bold", transform=tl)
        ax.text(1.03, yy, f"{m:.2f} [{l:.2f}, {h:.2f}]", va="center", fontsize=8.5, weight="bold", transform=tl)

    diamond(fe, 1.4, "#999999", "Efectos fijos (como el paper)")
    diamond(re, 0.5, "#5B2C83", "Efectos aleatorios REML-HKSJ")
    if not np.isnan(re["pi_lo"]):
        ax.plot([np.exp(re["pi_lo"]), np.exp(re["pi_hi"])], [-0.2, -0.2], color="#5B2C83", lw=3, alpha=0.35)
        ax.text(-0.03, -0.2, "Intervalo de predicción 95%", va="center", ha="right", fontsize=8.5, style="italic", transform=tl)
    ax.axvline(1, color="black", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(0.005, 200)
    ax.set_xticks([0.05, 0.2, 1, 5, 20])
    ax.set_xticklabels(["0.05", "0.2", "1", "5", "20"])
    ax.set_ylim(-0.8, k + 3)
    ax.set_yticks([])
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.set_xlabel("Odds ratio (← favorece prehabilitación | favorece control →)")
    fig.suptitle(f"{title}   (I² = {fe['I2']:.0f}%, τ² = {re['tau2']:.3f}, k = {k}; morado = ECA)", fontsize=11, x=0.02, ha="left")
    fig.subplots_adjust(left=0.33, right=0.76, top=0.9, bottom=0.12)
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    eff = pd.read_csv(DATA / "study_effects.csv")
    chars = pd.read_csv(DATA / "study_characteristics.csv")
    df = eff.merge(chars, on="study", how="left", validate="many_to_one")
    assert df.design.notna().all(), "estudio sin características"
    df["neoadjuvant"] = np.where(df.treatment == "S", "Solo cirugía", "Neoadyuvancia + cirugía")
    df["cancer_grp"] = df.cancer.replace({"ESCC": "EC", "GOJ": "EC"})

    summary, loo_rows, sub_rows = [], [], []
    for outc, sub in df.groupby("outcome", sort=False):
        y, se = sub.logOR.values, sub.SE.values
        fe, re = fixed_effect(y, se), random_effects(y, se)
        row = dict(outcome=outc, label=LABELS[outc], k=fe["k"],
                   k_RCT=int((sub.design_group == "RCT").sum()),
                   OR_FE=np.exp(fe["mu"]), lo_FE=np.exp(fe["lo"]), hi_FE=np.exp(fe["hi"]), p_FE=fe["p"],
                   OR_RE=np.exp(re["mu"]), lo_RE=np.exp(re["lo"]), hi_RE=np.exp(re["hi"]), p_RE=re["p"],
                   PI_lo=np.exp(re["pi_lo"]), PI_hi=np.exp(re["pi_hi"]),
                   I2=fe["I2"], tau2=re["tau2"])
        obs = sub[sub.design_group == "Observational"]
        if len(obs) >= 2:
            ro = random_effects(obs.logOR.values, obs.SE.values)
            row.update(OR_obs_RE=np.exp(ro["mu"]), lo_obs=np.exp(ro["lo"]), hi_obs=np.exp(ro["hi"]))
        rct = sub[sub.design_group == "RCT"]
        if len(rct) >= 1:
            rr = fixed_effect(rct.logOR.values, rct.SE.values)
            row.update(OR_RCT=np.exp(rr["mu"]), lo_RCT=np.exp(rr["lo"]), hi_RCT=np.exp(rr["hi"]))
        if fe["k"] >= 8:
            b0, pe = egger(y, se)
            row.update(egger_intercept=b0, egger_p=pe)
        summary.append(row)

        forest(sub, fe, re, LABELS[outc], FIG / f"forest_{outc}.png")

        if fe["k"] >= 4:
            for st in sub.study:
                s2 = sub[sub.study != st]
                r2 = random_effects(s2.logOR.values, s2.SE.values)
                loo_rows.append(dict(outcome=outc, omitted=st, OR_RE=np.exp(r2["mu"]),
                                     lo=np.exp(r2["lo"]), hi=np.exp(r2["hi"]), p=r2["p"]))
        if outc in ("overall_complications", "pulmonary_complications"):
            for gc, gname in [("design_group", "Diseño"), ("cancer_grp", "Tipo de cáncer"),
                              ("neoadjuvant", "Tratamiento"), ("modality", "Modalidad")]:
                tab, qb, dfb, p = q_between(sub, gc)
                tab.insert(0, "subgroup_var", gname)
                tab.insert(0, "outcome", outc)
                tab["Q_between"], tab["df"], tab["p_between"] = qb, dfb, p
                sub_rows.append(tab)

    S = pd.DataFrame(summary)
    S.to_csv(RES / "summary_by_outcome.csv", index=False, float_format="%.3f")
    pd.DataFrame(loo_rows).to_csv(RES / "leave_one_out.csv", index=False, float_format="%.3f")
    pd.concat(sub_rows).to_csv(RES / "subgroups.csv", index=False, float_format="%.3f")

    # Figura resumen
    S2 = S.iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for i, r in S2.iterrows():
        ax.plot([r.lo_RE, r.hi_RE], [i, i], color="#5B2C83", lw=2)
        ax.plot(r.OR_RE, i, "o", color="#5B2C83")
        ax.plot([r.lo_FE, r.hi_FE], [i - 0.25, i - 0.25], color="#AAAAAA", lw=1.5)
        ax.plot(r.OR_FE, i - 0.25, "s", color="#AAAAAA", ms=4)
    ax.set_yticks(range(len(S2)))
    ax.set_yticklabels([f"{l} (k={k})" for l, k in zip(S2.label, S2.k)])
    ax.axvline(1, color="black", lw=0.8)
    ax.set_xscale("log")
    ax.set_xticks([0.05, 0.2, 0.5, 1, 2, 5, 20])
    ax.set_xticklabels(["0.05", "0.2", "0.5", "1", "2", "5", "20"])
    ax.set_xlabel("OR (← favorece prehabilitación)")
    ax.set_title("Prehabilitación vs. control · morado: REML-HKSJ · gris: efectos fijos", loc="left", fontsize=10)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "summary_all_outcomes.png", dpi=150)
    plt.close(fig)

    pd.set_option("display.width", 220)
    print(S[["label", "k", "k_RCT", "OR_FE", "lo_FE", "hi_FE", "OR_RE", "lo_RE", "hi_RE", "p_RE",
             "PI_lo", "PI_hi", "I2"]].round(3).to_string(index=False))
    print(S[["label", "OR_obs_RE", "lo_obs", "hi_obs", "OR_RCT", "lo_RCT", "hi_RCT",
             "egger_intercept", "egger_p"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
