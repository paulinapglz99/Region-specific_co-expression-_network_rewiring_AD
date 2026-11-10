# Diseño del paper

> Propuesta. Los resultados todavía no existen: todo lo que dice "esperamos" es una hipótesis, y la sección "Si los resultados no salen como esperamos" dice qué cambia en cada caso. Los análisis se refieren a `docs/experimental_design.md`.

## Aportación biológica

1. **El cambio en la coordinación génica de AD depende de la región.** Comparar redes entre regiones está sesgado por el tamaño de muestra: con más personas se detectan más enlaces. Por eso, hasta ahora no se ha podido decir con rigor qué región se recablea más. Este trabajo compara el rewiring AD vs CT entre cinco regiones con el n igualado y contra una referencia de cuánto se parecen dos redes sin diferencia biológica (J_ref).
2. **El rewiring sigue la vulnerabilidad regional de AD**, con una predicción registrada antes de inferir las redes: TCX tiene más rewiring que CBE.
3. **La red detecta genes que la expresión diferencial no ve:** genes que no cambian de nivel en AD, pero sí cambian de vecinos en la red. Además, se describe qué hubs se ganan o se pierden en cada región y si son genes DE.

**Aportación metodológica de apoyo:** un marco para comparar redes de co-expresión entre grupos con distinto n:
- n igualado;
- curva de J de referencia validada en grupos independientes;
- índice de rewiring;
- una sola prueba para cruzar conjuntos de genes, que controla la expresión (y el grado) y no supone genes independientes;
- validación con redes simuladas con respuesta conocida (prueba interna: no va como figura; se agrega si un revisor la pide).

Se presenta como lo que hace posible el resultado biológico, no como el resultado principal.

## Resultado principal

**Hipótesis:** la magnitud del rewiring AD vs CT difiere entre regiones del cerebro después de controlar el tamaño de muestra, la proporción de sexos y qué donantes se eligen. Además, la región que se afecta antes en AD (TCX) tiene más rewiring que la región preservada (CBE).

En una frase, para el abstract: *"El rewiring de las redes de co-expresión en AD es específico de región, sigue la vulnerabilidad regional y afecta a genes cuyo nivel de expresión no cambia."*

Esta frase solo se puede escribir si los análisis 1, 4 y 6 salen en la dirección esperada. Si no, ver "Si los resultados no salen como esperamos".

## Historia: secuencia de resultados

Cada paso responde una pregunta que deja abierto el anterior.

| # | Pregunta | Qué mostramos | Análisis | Figura |
|---|---|---|---|---|
| 1 | ¿Se pueden comparar redes de regiones con distinto n? | No directamente: J sube con n aunque no haya diferencia biológica. Medimos esa dependencia (J_ref) y la validamos en grupos independientes | 3 | 1a |
| 2 | ¿Cuánto rewiring hay en cada región? | Índice de rewiring por región, con el n igualado | 1, 3 | 1b |
| 3 | ¿Esto es biología de AD o ruido? | Se cumple la predicción registrada: TCX tiene más rewiring que CBE | 6 | 1b |
| 4 | ¿Qué cambia en la estructura de la red, además de los enlaces? | Grado, comunidades y densidad | 1 | 1c, 1d |
| 5 | ¿Los genes que se recablean son los genes DE? | En parte sí y en parte no. Hay genes recableados que no son DE | 4 | 2a |
| 6 | ¿Qué hacen los genes que solo detecta la red? | Sus funciones (enriquecimiento GO) | 4 | 2b |
| 7 | ¿Los hubs que cambian son DE? | Hubs compartidos, ganados y perdidos por región, contra el DE | 5 | 2c |

Pasos 1–4 → **Figura 1** · pasos 5–7 → **Figura 2**. La robustez (análisis 2, sensibilidad al umbral, estabilidad entre réplicas) y los genes de riesgo (análisis 7) van al suplemento.

## Figura 1. El rewiring de las redes de co-expresión en AD es específico de región y sigue la vulnerabilidad regional

**Mensaje:** con una comparación justa entre regiones, AD recablea las redes de forma distinta en cada región, y más en la región que se afecta antes (TCX) que en la preservada (CBE).

| Panel | Qué muestra | Análisis | De dónde sale | Estado |
|---|---|---|---|---|
| a | Curva J_ref(n): J entre dos mitades del mismo grupo, calibración (DLPFC AD) y validación (otros 3 grupos) | 3 | `reference_curve.tsv` | existe |
| b | **Índice de rewiring por región** (10 réplicas, promedio ± DE), línea en 1 = sin rewiring, con la predicción registrada TCX > CBE y su p de Wilcoxon. **Panel principal del paper** | 1, 3, 6 | `rewiring_index.tsv`, `prediction.tsv` | existe |
| c | Estructura: ρ de grado y NMI de comunidades, AD vs CT | 1 | `rewiring_index.tsv` | existe |
| d | Grado medio de las redes AD y CT | 1 | `rewiring_index.tsv` | existe |

Se genera con la regla `plot_network_metrics`, en `results/figures/figure1/` (`figure1.pdf` y un archivo por panel: `1a.pdf` … `1d.pdf`).

## Figura 2. La red detecta cambios que la expresión diferencial no ve

**Mensaje:** el rewiring y la expresión diferencial se solapan solo en parte. Los genes recableados que no son DE tienen funciones propias, y los hubs que cambian se comparan con el DE en cada región.

| Panel | Qué muestra | Análisis | De dónde sale | Estado |
|---|---|---|---|---|
| a | Genes recableados (top 10%) de cada región: fracción DE ("recableado y DE") y no DE ("solo recableado"), con una línea para la fracción DE esperada por azar, y el OR con su q encima de cada barra | 4 | `rewiring_tests.tsv` | existe |
| b | Términos GO de los genes "solo recableados": los 5 más significativos de cada región | 4 | `rewiring_go.tsv` | existe |
| c | Hubs compartidos, ganados y perdidos: fracción de genes DE observada contra la esperada por azar, por región | 5 | `hubs_tests.tsv` | existe |

Se genera con la regla `plot_networks_de`, en `results/figures/figure2/` (`figure2.pdf`, `2a.pdf`, `2b.pdf` y `2c.pdf`).

**Cómo se lee el panel a:** "recableado" es el top 10% por definición, así que todas las barras suman el 100% de los genes recableados. Lo que importa es dónde queda la parte DE respecto de la línea del esperado: por encima, los recableados son DE más de lo esperado (OR > 1). La parte "solo recableado" son los genes que la red detecta y el DE no. La línea del esperado es necesaria porque el número de genes DE cambia mucho entre regiones (TCX tiene muchos más).

## Figuras y tablas suplementarias

| # | Contenido | Estado |
|---|---|---|
| S1 | Demografía de los grupos (sexo, edad, APOE, RIN, donantes compartidos entre regiones) | existe (`results/figures/figureS1/`) |
| S2 | Diagrama completo del diseño (análisis 1–7) | existe (`docs/figures/experimental_design`) |
| S3 | J directa AD vs CT por región (con la predicción sobre la J directa), e índice con n = 90 en las 3 regiones de ROSMAP (análisis 2) | existe (`results/figures/figureS3/`) |
| S4 | Robustez del análisis 4: OR con el top 5% y el top 20%, y rango del OR entre réplicas | existe (`results/figures/figureS4/`) |
| T1 | Genes recableados robustos por región (top 10% en ≥ 8 de 10 réplicas), ordenados por rewiring ajustado: los genes con más rewiring | tabla (`rewiring_robust.tsv`) |
| T2 | Términos GO completos de los genes "solo recableados" | tabla (`rewiring_go.tsv`) |
| T3 | Lista de hubs ganados y perdidos por región | tabla (`hubs_genes.tsv`) |
| T4 | Genes de riesgo de AD (análisis 7): pruebas contra recableados y hubs, y tabla por gen y región | tablas (`risk_genes_tests.tsv`, `risk_genes_table.tsv`) |

**Validación con redes simuladas** (`planted`: control positivo; `null`: control negativo, con AD = CT). Es una prueba interna del código y no va como figura. Sus resultados están en `results/simulated/{planted,null}/` y se agregan al suplemento solo si un revisor los pide.

## Si los resultados no salen como esperamos

| Resultado | Qué cambia en el paper |
|---|---|
| El índice es ≈ 1 en todas las regiones (no hay rewiring) | El mensaje pasa a ser que AD cambia el nivel de expresión pero no la coordinación entre genes. La Figura 1 sigue igual. La Figura 2 se reduce a mostrar que los genes recableados no se distinguen de los demás |
| Hay rewiring, pero igual en todas las regiones | El mensaje pasa a ser que el rewiring de AD es común a todo el cerebro |
| No se cumple TCX > CBE | Se reporta igual (la predicción está registrada). La discusión trata por qué el rewiring no sigue la vulnerabilidad, por ejemplo por la composición celular |
| Rewiring y DE se solapan casi por completo | La red no añade información sobre qué genes cambian. La aportación se queda en la comparación entre regiones (Figura 1) |
| No hay términos GO significativos en los "solo recableados" | El panel 2b se reemplaza por la tabla de esos genes en el suplemento |
| No hay enriquecimiento de genes de riesgo | Se reporta la tabla descriptiva (T4) y el bajo poder de la prueba |

## Decisiones pendientes

1. **Enfoque:** biología primero, con el marco metodológico como apoyo (propuesta), o metodología primero.
2. **Revista objetivo:** define el número de figuras, el tamaño de los paneles y cuánto va al suplemento.
3. **Panel de esquema del diseño en la Figura 1:** agregarlo, o dejar el diseño solo en S2.

## Decisiones tomadas

- **Enriquecimiento funcional:** GO Biological Process (`config/config.yaml`, `go`).
- **Genes de riesgo:** solo los 55 genes de los 42 loci nuevos de Bellenguez et al. (2022); los 33 loci ya conocidos no se incluyen.
- **Red consenso:** se quitó del diseño, porque ningún análisis la usa.
- **Validación con redes simuladas:** prueba interna, sin figura en el paper.
