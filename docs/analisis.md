# Análisis

Resumen de los análisis del diseño (`docs/experimental_design.md`). Los resultados esperados son hipótesis: todavía no hay resultados con las redes reales.

| Análisis | Pregunta | Valores de entrada | Valores de salida | Resultado esperado | Redes nuevas | Figura |
|---|---|---|---|---|---|---|
| 1. AD vs CT con el mismo n (`matched`) | ¿Cuánto rewiring hay entre AD y CT en cada región? | Por región: la red `CT_all` y las 10 redes `AD_rep01` … `AD_rep10`; J_ref(n) del análisis 3 | Por réplica: J(AD, CT), ρ de Spearman del grado, NMI de las comunidades de Leiden, grado medio de AD y de CT e **índice de rewiring** = J(AD, CT) / J_ref(n) | Índice de rewiring < 1 en las 5 regiones (hay rewiring), con un valor distinto en cada región y poca dispersión entre las 10 réplicas | 55 | 1c, 1d, 1e |
| 2. Mismo n en todas las regiones (`equal_n`) | ¿El orden de las regiones por rewiring depende del n? | Por región de ROSMAP: 10 pares de redes `AD_repk` y `CT_repk` con n = 90; J_ref(90) del análisis 3 | Por réplica: J(AD, CT), índice de rewiring y orden de las 3 regiones. En total: % de réplicas con el mismo orden que el análisis 1 y **τ de Kendall** contra ese orden | El mismo orden que en el análisis 1 en la mayoría de las 10 réplicas (τ de Kendall cercana a 1) | 60 | S3b |
| 3. J de referencia (`curve`) | ¿Cuánto se parecen dos redes sin ninguna diferencia biológica? | Por grupo y n: 5 pares de redes de dos mitades sin personas en común (`repk_A` y `repk_B`) | Por réplica: J(A, B). Del modelo: los coeficientes a y b de **J_ref(n) = a + b · log(n)**, y J_ref en el n de cada comparación de los análisis 1 y 2 | J_ref sube con n, siguiendo a + b · log(n); los 3 grupos de validación caen sobre la curva | 60 | 1b |
| 4. Rewiring y DE | ¿Los genes recableados son DE, o la red detecta genes que el DE no ve? ¿Cuáles son los genes con más rewiring? | Las redes del análisis 1; el DE de Sage (`logFC`, `t`, `AveExpr`, `adj.P.Val`) | Por gen: cambio de vecinos, grado en AD y en CT, **rewiring ajustado**, recableado (sí/no), réplicas en el top 10% (recableado robusto) y grupo (recableado y DE, solo recableado, solo DE, sin cambio). Por región: **OR** recableados × DE, genes en común observados y esperados, p y q; OR con el top 5% y el 20%; OR de cada réplica; términos GO con su q | Solapamiento parcial: OR > 1, pero con muchos genes recableados no DE, con términos GO propios | 0 (usa las del 1) | 2a, 2b, S4 |
| 5. Hubs y DE | ¿Los hubs que se ganan o se pierden en AD son DE? | Las redes del análisis 1; el DE de Sage | Por gen: **PageRank** en CT, número de réplicas AD en que es hub y tipo de hub (compartido, ganado, perdido, no se clasifica). Por región y tipo de hub: **OR** tipo × DE, % DE observado y esperado, p y q | Hubs ganados y perdidos con más genes DE de lo esperado (OR > 1); hubs compartidos cerca de lo esperado | 0 (usa las del 1) | 2c |
| 6. Vulnerabilidad regional (predicción registrada el 2026-11-06) | ¿La región que se afecta antes (TCX) tiene más rewiring que la preservada (CBE)? | El índice de rewiring y la J de las 10 réplicas de TCX y de CBE (análisis 1 y 3) | Con el índice y con la J directa: promedio de cada región, estadístico **W y p de Wilcoxon**, y si la predicción se cumple (sí/no) | Índice de TCX menor que el de CBE, con Wilcoxon p < 0.05 entre las 10 réplicas de cada una | 0 (usa las del 1 y el 3) | 1c, S3a |
| 7. Genes de riesgo | ¿Los genes de riesgo de AD (Bellenguez et al., 2022) están entre los que cambian en la red? | Los 55 genes de riesgo; los genes recableados (análisis 4) y los hubs ganados y perdidos (análisis 5) | Por región y conjunto (recableados, hubs ganados, hubs perdidos): **OR**, p y q. Por gen de riesgo y región: percentil de rewiring ajustado, recableado (sí/no), tipo de hub y si es DE | Poco poder (55 genes): se espera una tendencia, no necesariamente significativa; la tabla por gen es el resultado principal | 0 (usa las del 1) | Tabla T4 |
| **Total** | | | | | **175** | |

Las 175 redes son de ARACNe-AP, con 100 bootstraps cada una. Los análisis 4–7 no necesitan redes nuevas.

**Abreviaturas:**
- **J:** índice de Jaccard de los enlaces (enlaces en ambas redes / enlaces en al menos una).
- **NMI:** información mutua normalizada.
- **OR:** odds ratio.
- **p:** en los análisis 4, 5 y 7, valor p por permutación (1000 permutaciones dentro de estratos de expresión y, cuando el conjunto son genes recableados, también de grado). En el análisis 6, valor p de la prueba de Wilcoxon.
- **q:** p corregido por Benjamini–Hochberg.
