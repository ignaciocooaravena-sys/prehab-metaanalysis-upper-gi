# Prehabilitación vs. no prehabilitación en cirugía oncológica GI alta — re-análisis reproducible

Re-análisis abierto y reproducible de los datos por estudio publicados en:

> Boosting recovery before surgery: The impact of prehabilitation on upper gastrointestinal cancer patients – A quantitative comparative analysis. *PLOS ONE* 20(3): e0315734 (2025). https://doi.org/10.1371/journal.pone.0315734 — licencia CC-BY 4.0.

12 estudios (2 ECA + 10 observacionales), 832 pacientes con cáncer de esófago, unión gastroesofágica o gástrico sometidos a cirugía, comparando **prehabilitación** (ejercicio ± nutrición ± apoyo psicológico) vs. **cuidado habitual**.

## ¿Por qué datos agregados y no por paciente?

Se buscaron datasets abiertos de pacientes individuales (IPD) con brazo prehabilitado y brazo control en Dryad, Zenodo, Figshare, Harvard Dataverse, Mendeley Data y material suplementario de PLOS ONE (noviembre 2026). **No se encontró ninguno de acceso libre**: los ensayos grandes comparten IPD solo "a solicitud" o vía Vivli, y el único dataset abierto de prehabilitación en Dryad (Hôpital Foch, 2021) es de un solo brazo (todos prehabilitados, comparando por edad). Por eso se usaron las estimaciones por estudio (log OR y SE) publicadas en las Figs 2–4 del paper, que sí son abiertas.

## Qué agrega este re-análisis

El paper original usó **solo efectos fijos**. Aquí se agregan:

1. **Efectos aleatorios** REML con ajuste Hartung-Knapp-Sidik-Jonkman (más apropiado con pocos estudios heterogéneos en diseño).
2. **Intervalos de predicción 95%** (rango esperable del efecto en un nuevo centro).
3. **Sensibilidad por diseño** (ECA vs. observacional) y **leave-one-out**.
4. **Subgrupos** (diseño, tipo de cáncer, neoadyuvancia, modalidad) con test Q entre grupos.
5. **Test de Egger** para outcomes con k ≥ 8.

Validación: el modelo de efectos fijos reproduce exactamente los OR del paper (p. ej. complicaciones globales 0,59 [0,39–0,88]; pulmonares 0,54 [0,36–0,80]).

## Resultados principales

| Outcome | k (ECA) | OR efectos fijos (paper) | OR REML-HKSJ | Intervalo de predicción 95% | I² |
|---|---|---|---|---|---|
| Complicaciones globales | 8 (1) | 0,59 [0,39–0,88] | **0,59 [0,36–0,96]**, p = 0,04 | 0,36–0,98 | 0% |
| Complicaciones pulmonares | 9 (1) | 0,54 [0,36–0,80] | **0,56 [0,29–1,09]**, p = 0,08 | 0,21–1,52 | 40% |
| Filtración anastomótica | 7 (0) | 0,94 [0,46–1,92] | 0,92 [0,35–2,43] | 0,24–3,54 | 0% |
| Complicaciones cardíacas | 4 (0) | 0,79 [0,35–1,75] | 0,79 [0,21–2,89] | — | 0% |
| Reingreso hospitalario | 4 (1) | 1,06 [0,52–2,17] | 1,06 [0,31–3,62] | — | 0% |
| Reingreso a UCI | 2 (1) | 0,23 [0,06–0,89] | no estimable con 2 estudios | — | 0% |

Herida, parálisis de cuerda vocal, quilotórax y sangrado (k = 3–4): sin diferencias, IC muy amplios. Tabla completa en [`results/summary_by_outcome.csv`](results/summary_by_outcome.csv).

![Resumen](figures/summary_all_outcomes.png)

### Sensibilidad

- **Solo ECA**: complicaciones globales OR 0,54 [0,16–1,80] (1 ECA, n = 51); pulmonares OR 0,40 [0,12–1,33] (1 ECA, n = 48). La dirección coincide, pero ningún ECA tiene potencia por sí solo.
- **Leave-one-out (REML-HKSJ)**: en complicaciones globales, al retirar 6 de 8 estudios el IC cruza 1 (OR entre 0,52 y 0,64). En pulmonares, 7 de 9 omisiones dejan el resultado no significativo.
- **Subgrupo pulmonar por tratamiento**: neoadyuvancia + cirugía OR 0,47 [0,25–0,87] vs. solo cirugía OR 1,79 (FE) — Q entre grupos p = 0,045. El subgrupo "solo cirugía" son 2 estudios japoneses de gástrico, así que tratamiento, tipo de cáncer y país están confundidos.
- **Egger**: sin asimetría detectable (p = 0,39 y 0,30), pero con 8–9 estudios el test tiene muy poca potencia.

## Extensión: Janssen et al. 2019 (colorrectal + AAA, ≥70 años)

> Janssen TL, Steyerberg EW, et al. Multimodal prehabilitation to reduce the incidence of delirium and other adverse events in elderly patients undergoing elective major abdominal surgery: An uncontrolled before-and-after study. *PLOS ONE* 14(6): e0218152 (2019). https://doi.org/10.1371/journal.pone.0218152 — CC-BY 4.0.

Estudio antes-después unicéntrico (Amphia, Países Bajos): 267 prehabilitados (2015–2018) vs. 360 controles históricos (2013–2015), cáncer colorrectal (n = 484) o aneurisma de aorta abdominal (n = 143). Prehabilitación domiciliaria multimodal (ejercicio, nutrición, fragilidad, hierro EV si anemia), duración media 39 días, adherencia 74%.

**Como la población no es GI alta, no se mezcla con el análisis principal.** Se usa de dos formas ([`scripts/extension_janssen2019.py`](scripts/extension_janssen2019.py)):

**1. Resultados del estudio, recalculados desde los conteos** ([`results/janssen2019_recalculated.csv`](results/janssen2019_recalculated.csv))

| Outcome | Prehab | Control | OR crudo (recalculado) | OR ajustado (reportado) |
|---|---|---|---|---|
| Delirium (primario) | 8,2% | 11,7% | 0,68 [0,40–1,17] | **0,56 [0,32–0,98]** |
| Cualquier complicación | 40,8% | 36,9% | 1,18 [0,85–1,63] | 1,12 [0,80–1,57] |
| Clavien-Dindo III–V | 14,2% | 8,9% | **1,70 [1,03–2,80]** | — |
| Ingreso no planificado a UCI | 11,2% | 7,5% | 1,56 [0,90–2,70] | 1,54 [0,86–2,75] |
| Reingreso 30 días | 8,4% | 6,1% | 1,41 [0,76–2,60] | 1,42 [0,75–2,68] |
| Filtración anastomótica (solo CRC) | 7,1% | 2,8% | **2,67 [1,10–6,49]** | — |
| Mortalidad hospitalaria/30 d | — | — | — | 1,50 [0,61–3,72] |

Nota: el OR crudo de "cualquier complicación" calculado desde los conteos de la Tabla 2 (1,18) no coincide con el "no ajustado" de la Tabla 3 (1,09); probablemente la regresión excluyó casos con datos faltantes.

**2. Pool exploratorio "cirugía abdominal mayor" (GI alto + Janssen 2019)** ([`results/extension_pooled_with_janssen2019.csv`](results/extension_pooled_with_janssen2019.csv))

| Outcome | Solo GI alto (REML-HKSJ) | + Janssen 2019 (REML-HKSJ) | Efectos fijos con Janssen |
|---|---|---|---|
| Complicaciones globales | 0,59 [0,36–0,96], p = 0,04 | **0,71 [0,44–1,14], p = 0,13**, I² 37% | 0,86 [0,66–1,12] |
| Filtración anastomótica | 0,92 [0,35–2,43] | 1,22 [0,51–2,91] | 1,42 [0,81–2,47] |
| Reingreso hospitalario | 1,06 [0,31–3,62] | 1,25 [0,64–2,45] | 1,25 [0,78–2,01] |

Para Janssen se usó el OR ajustado cuando existía (complicaciones, reingreso) y el crudo para filtración (solo CRC).

![Extensión](figures/extension_forest_overall_complications.png)

**Lectura**: Janssen 2019 tiene más pacientes (627) que cualquier estudio del meta-análisis GI alto (máximo 111), y apunta en dirección contraria para complicaciones. Al sumarlo, el beneficio en complicaciones globales **deja de ser significativo** y aparece heterogeneidad (I² 0% → 37%). Su único hallazgo favorable es el delirium (OR ajustado 0,56), un outcome que los estudios GI alto no reportan. Su grupo prehabilitado era más enfermo al inicio (CCI ≥ 7: 41% vs. 31%; ASA ≥ 3: 56% vs. 41%) y fue operado en un período posterior, por lo que el diseño antes-después no permite separar el efecto de la prehabilitación del de los cambios en el tiempo.

## Conclusiones

1. **La dirección del efecto es consistente**: la prehabilitación se asocia a ~40% menos odds de complicaciones postoperatorias globales y ~45% menos odds de complicaciones pulmonares.
2. **La certeza es menor de lo que sugiere el paper**. Con un modelo de efectos aleatorios, el efecto sobre complicaciones pulmonares deja de ser significativo (IC hasta 1,09) y el intervalo de predicción incluye daño (hasta 1,52). El efecto global sobrevive (p = 0,04) pero es frágil al leave-one-out.
3. **La evidencia es casi toda observacional** (10/12 estudios). La prehabilitación no se asigna al azar en esos estudios: los pacientes más motivados o en mejor forma pueden tener más probabilidad de recibirla, lo que sesga a favor de la intervención.
4. **No hay señal de efecto** en filtración anastomótica, complicaciones cardíacas, reingreso ni otras complicaciones específicas.
5. **El beneficio pulmonar parece concentrarse en esófago con neoadyuvancia**, lo cual es biológicamente plausible (la neoadyuvancia deteriora capacidad funcional y la prehabilitación puede contrarrestarlo), pero es una hipótesis, no una conclusión.
6. **Al agregar un estudio grande de cirugía colorrectal/AAA en adultos mayores (Janssen 2019)**, el beneficio sobre complicaciones globales desaparece (OR 0,71 [0,44–1,14]). La prehabilitación no muestra un beneficio robusto sobre complicaciones en cirugía abdominal mayor en general; su efecto podría depender del tipo de cirugía, de la población y de la intensidad o supervisión del programa (Janssen fue domiciliario y no supervisado). El único beneficio en esa población fue menor delirium.
7. **Implicancia práctica**: el resultado apoya programas multimodales de prehabilitación en esofagectomía post-neoadyuvancia como intervención de bajo riesgo con beneficio probable, y refuerza la necesidad de ECA con potencia y de que se compartan datos individuales.

## Limitaciones

- Datos extraídos de las figuras del paper (log OR y SE redondeados a 2 decimales); no se re-extrajeron de los estudios primarios.
- Los OR de estudios observacionales son mayormente crudos (no ajustados por confusores).
- Posible solapamiento de pacientes entre cohortes del mismo grupo o centro (p. ej. pares de estudios con el mismo autor en años distintos): no verificable con los datos disponibles.
- Subgrupos con 1–3 estudios: exploratorios.

## Estructura

```
data/
  study_effects.csv           # log OR y SE por estudio y outcome (Figs 2-4 del paper)
  study_characteristics.csv   # diseño, país, cáncer, tratamiento, modalidad, n (Tabla 1)
  janssen2019_outcomes.csv    # conteos y OR reportados de Janssen 2019 (Tablas 2-3)
scripts/
  meta_analysis.py            # análisis principal (FE, REML-HKSJ, PI, LOO, subgrupos, Egger)
  extension_janssen2019.py    # recálculo de Janssen 2019 + pool exploratorio
results/
  summary_by_outcome.csv
  subgroups.csv
  leave_one_out.csv
  janssen2019_recalculated.csv
  extension_pooled_with_janssen2019.csv
figures/
  forest_*.png                # forest plot por outcome (GI alto)
  extension_forest_*.png      # forest plots con Janssen 2019
  summary_all_outcomes.png
```

## Reproducir

```bash
pip install -r requirements.txt
python scripts/meta_analysis.py
python scripts/extension_janssen2019.py
```

## Licencia

Código: MIT. Datos derivados de los artículos originales bajo CC-BY 4.0 — citar las fuentes originales.
