# Diseño del paper

> Propuesta. Los resultados todavía no existen: todo lo que dice "esperamos" es una hipótesis, y la sección "Si los resultados no salen como esperamos" dice qué cambia en cada caso. Los análisis se refieren a `docs/experimental_design.md`.

## Aportación biológica

1. **El cambio en la coordinación génica de AD depende de la región.** Comparar redes entre regiones está sesgado por el tamaño de muestra: con más personas se detectan más enlaces. Por eso, hasta ahora no se ha podido decir con rigor qué región se recablea más. Este trabajo compara el rewiring AD vs CT entre cinco regiones con el n igualado y contra una referencia de cuánto se parecen dos redes sin diferencia biológica (J_ref).
2. **La red detecta cambios que la expresión diferencial no ve.** Genes que no cambian de nivel en AD, pero sí cambian de vecinos en la red. Además, se describe qué pasa alrededor de los genes con el cambio de expresión más fuerte (sus redes ego) y qué genes ganan o pierden centralidad (hubs) en cada región.
3. **El rewiring concuerda con la biología conocida de AD:**
   - sigue la vulnerabilidad regional, con una predicción registrada antes de inferir las redes (TCX > CBE);
   - toca a los genes de riesgo genético de AD.

**Aportación metodológica de apoyo:** un marco para comparar redes de co-expresión entre grupos con distinto n:
- n igualado;
- curva de J de referencia validada en grupos independientes;
- índice de rewiring;
- validación con redes simuladas con respuesta conocida.

Se presenta como lo que hace posible el resultado biológico, no como el resultado principal.

## Resultado principal

**Hipótesis:** la magnitud del rewiring AD vs CT difiere entre regiones del cerebro después de controlar el tamaño de muestra, la proporción de sexos y qué donantes se eligen. Además, la región que se afecta antes en AD (TCX) tiene más rewiring que la región preservada (CBE).

En una frase, para el abstract: *"El rewiring de las redes de co-expresión en AD es específico de región, sigue la vulnerabilidad regional y afecta a genes cuyo nivel de expresión no cambia."*

Esta frase solo se puede escribir si los análisis 1, 4 y 6 salen en la dirección esperada. Si no, ver "Si los resultados no salen como esperamos".

## Historia: secuencia de resultados

Cada paso responde una pregunta que deja abierto el anterior.

| # | Pregunta | Qué mostramos | Análisis | Figura |
|---|---|---|---|---|
| 1 | ¿Se pueden comparar redes de regiones con distinto n? | No directamente: J sube con n aunque no haya diferencia biológica. Medimos esa dependencia (J_ref) y la validamos en grupos independientes | 3 | 1b |
| 2 | ¿Cuánto rewiring hay en cada región? | Índice de rewiring por región, con el n igualado | 1, 3 | 1c |
| 3 | ¿Ese resultado depende del n o de los donantes? | El orden de las regiones se mantiene con n = 90 en todas, y la dispersión entre réplicas es baja | 2 | 1d |
| 4 | ¿Qué cambia en la estructura de la red, además de los enlaces? | Grado, comunidades y densidad | 1 | 1e |
| 5 | ¿Los genes que se recablean son los genes DE? | En parte sí y en parte no. Hay genes recableados que no son DE | 4 | 2a, 2b |
| 6 | ¿Qué pasa alrededor de los genes más DE? | Sus vecindarios cambian más que los de genes comparables, y se rodean de otros genes DE | 4 (redes ego) | 2c |
| 7 | ¿Qué genes ganan o pierden centralidad, y en qué regiones? | Hubs ganados y perdidos por región; hubs exclusivos de una región, más de lo esperado por azar | 5 | 2d, 2e |
| 8 | ¿Esto es biología de AD o ruido? | Se cumple la predicción registrada: TCX tiene más rewiring que CBE | 6 | 3a |
| 9 | ¿Toca a los genes de riesgo de AD? | Los genes de riesgo están entre los recableados o los hubs que cambian | 7 | 3b, 3c |

Pasos 1–4 → **Figura 1** · pasos 5–7 → **Figura 2** · pasos 8–9 → **Figura 3**.

## Figura 1. El rewiring de las redes de co-expresión en AD es específico de región

**Mensaje:** con una comparación justa entre regiones, AD recablea las redes de forma distinta en cada región.

| Panel | Qué muestra | Análisis | De dónde sale | Estado |
|---|---|---|---|---|
| a | Esquema del diseño: 5 regiones de 2 cohortes, n de AD y CT, submuestreo con el mismo n, ARACNe-AP y comparación con J | — | versión reducida de `docs/figures/experimental_design` | por adaptar |
| b | Curva J_ref(n): J entre dos mitades del mismo grupo, calibración (DLPFC AD) y validación (otros 3 grupos) | 3 | `reference_curve.tsv` | existe (`plot_network_metrics`) |
| c | **Índice de rewiring por región** (10 réplicas, promedio ± DE), línea en 1 = sin rewiring. **Panel principal del paper** | 1, 3 | `rewiring_index.tsv` | existe |
| d | Robustez: índice con n = 90 en las 3 regiones de ROSMAP, y % de réplicas con el mismo orden que en c | 2 | `rewiring_index.tsv`, `region_order.tsv` | existe |
| e | Estructura: ρ de grado y NMI de comunidades AD vs CT, y grado medio de AD y CT | 1 | `rewiring_index.tsv` | existe |

**Al suplemento:**
- la J directa por región, que hoy es el panel a de `plot_network_metrics`;
- las métricas réplica por réplica.

## Figura 2. La red detecta cambios que la expresión diferencial no ve

**Mensaje:** el rewiring y la expresión diferencial se solapan solo en parte. Los genes más DE cambian de vecindario, y los hubs que cambian son en parte exclusivos de cada región.

| Panel | Qué muestra | Análisis | Estado |
|---|---|---|---|
| a | Efecto de ser DE sobre el rewiring ajustado, por región: estimador e intervalo de confianza, uno por región (forest plot) | 4 | por hacer |
| b | Cuatro grupos de genes por región: recableado y DE, solo recableado, solo DE, sin cambio (barras apiladas). Resalta "solo recableado" | 4 | por hacer |
| c | Redes ego de los genes más DE: J del ego y fracción de vecinos DE, contra genes no DE de la misma expresión y el mismo grado. Arriba, las redes ego CT y AD de uno de los 50 genes, como ejemplo | 4 (redes ego) | por hacer |
| d | Hubs por región: compartidos, ganados y perdidos en AD, y su fracción de genes DE contra el nulo | 5 | por hacer |
| e | Hubs exclusivos: número de genes que son hub ganado (o perdido) en 1, 2… regiones, observado contra el esperado binomial | 5 | por hacer |

**Al suplemento:**
- la correlación de Spearman entre rewiring y \|t\|;
- la dirección del cambio de grado en DE arriba y DE abajo;
- la tabla de genes de cada grupo;
- el enriquecimiento funcional de los genes "solo recableados".

## Figura 3. El rewiring concuerda con la biología conocida de AD

**Mensaje:** el rewiring sigue la vulnerabilidad regional predicha de antemano y alcanza a los genes de riesgo genético de AD.

| Panel | Qué muestra | Análisis | Estado |
|---|---|---|---|
| a | Predicción registrada: índice de rewiring de TCX contra CBE (10 réplicas cada una), con el p de Wilcoxon. Al lado, lo mismo con la J directa | 6 | por hacer (los datos ya salen de `rewiring_index`) |
| b | Enriquecimiento de genes de riesgo (Bellenguez et al., 2022) entre los genes recableados y entre los hubs ganados y perdidos, por región: observado contra la distribución nula | 7 | por hacer |
| c | Mapa de calor de los genes de riesgo (filas) × regiones (columnas): percentil de rewiring como color, con marcas de hub y de DE | 7 | por hacer |
| d | Red ego CT y AD de un gen de riesgo, como ejemplo (uno de los genes de riesgo, elegido por tener el mayor rewiring) | 4, 7 | por hacer |

## Figuras suplementarias

| # | Contenido | Estado |
|---|---|---|
| S1 | Demografía de los grupos (sexo, edad, APOE, RIN, donantes compartidos entre regiones) | existe (`plot_demographics`) |
| S2 | Diagrama completo del diseño (análisis 1–7) | existe (`docs/figures/experimental_design`) |
| S3 | Validación del marco con redes simuladas: con rewiring plantado, el índice recupera el orden plantado; sin rewiring, el índice es ≈ 1 en todas las regiones | existe (`results/simulated/*/metrics/`) |
| S4 | J directa AD vs CT por región, y proporción de mujeres lograda en cada red | parcial |
| S5 | Resultados de la Figura 2 réplica por réplica (estabilidad) | por hacer |
| S6 | Enriquecimiento funcional de los genes "solo recableados" y de los vecinos ganados y perdidos de las redes ego | por hacer |

## Si los resultados no salen como esperamos

| Resultado | Qué cambia en el paper |
|---|---|
| El índice es ≈ 1 en todas las regiones (no hay rewiring) | El mensaje pasa a ser que AD cambia el nivel de expresión pero no la coordinación entre genes. La Figura 1 sigue igual. La Figura 2 se enfoca en los genes DE y sus redes ego |
| Hay rewiring, pero igual en todas las regiones | El mensaje pasa a ser que el rewiring de AD es común a todo el cerebro. La Figura 2 busca qué genes lo comparten |
| No se cumple TCX > CBE | Se reporta igual (la predicción está registrada). La discusión trata por qué el rewiring no sigue la vulnerabilidad, por ejemplo por la composición celular |
| Rewiring y DE se solapan casi por completo | La red no añade información sobre qué genes cambian. La aportación se queda en la comparación entre regiones (Figura 1) |
| No hay enriquecimiento de genes de riesgo | Se reporta la tabla descriptiva (3c) y el bajo poder de la prueba |

## Decisiones pendientes

1. **Enfoque:** biología primero, con el marco metodológico como apoyo (propuesta), o metodología primero.
2. **Revista objetivo:** define el número de figuras, el tamaño de los paneles y cuánto va al suplemento.
3. **Panel 1a:** esquema dibujado, o nada y que el diseño quede solo en S2.
4. **Base de datos para el enriquecimiento funcional** (S6): GO Biological Process u otra.
5. **Criterio para elegir los genes de ejemplo** de las redes ego (2c y 3d). Debe fijarse antes de ver los resultados, para no elegir el ejemplo más llamativo.
