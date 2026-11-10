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
| 0 | ¿Qué se compara? | Cinco regiones, una red CT y diez redes AD por región, comparadas con la J | 1 | 1a |
| 1 | ¿Se pueden comparar redes de regiones con distinto n? | No directamente: J sube con n aunque no haya diferencia biológica. Medimos esa dependencia (J_ref) y la validamos en grupos independientes | 3 | 1b |
| 2 | ¿Cuánto rewiring hay en cada región? | Índice de rewiring por región, con el n igualado | 1, 3 | 1c |
| 3 | ¿Esto es biología de AD o ruido? | Se cumple la predicción registrada: TCX tiene más rewiring que CBE | 6 | 1c |
| 4 | ¿Qué cambia en la estructura de la red, además de los enlaces? | Grado, comunidades y densidad | 1 | 1d, 1e |
| 5 | ¿Los genes que se recablean son los genes DE? | En parte sí y en parte no. Hay genes recableados que no son DE | 4 | 2a |
| 6 | ¿Qué hacen los genes que solo detecta la red? | Sus funciones (enriquecimiento GO) | 4 | 2b |
| 7 | ¿Los hubs que cambian son DE? | Hubs compartidos, ganados y perdidos por región, contra el DE | 5 | 2c |

Pasos 0–4 → **Figura 1** · pasos 5–7 → **Figura 2**. La robustez (análisis 2, sensibilidad al umbral, estabilidad entre réplicas) y los genes de riesgo (análisis 7) van al suplemento.

**Dónde están las figuras.** Cada figura tiene su carpeta en `results/figures/` (por ejemplo `results/figures/figure1/`), con la figura completa (`figure1.pdf`) y un archivo por panel (`1a.pdf`, `1b.pdf`…), en PDF y JPEG. Los paneles hechos a mano se guardan en la misma carpeta: el pipeline nunca borra la carpeta, solo escribe sus propios archivos.

Los pies de figura son una propuesta en inglés, como el texto de las figuras. Los números entre corchetes se completan con los resultados reales.

## Figura 1. El rewiring de las redes de co-expresión en AD es específico de región y sigue la vulnerabilidad regional

**Mensaje:** con una comparación justa entre regiones, AD recablea las redes de forma distinta en cada región, y más en la región que se afecta antes (TCX) que en la preservada (CBE).

Se genera con la regla `plot_network_metrics`, salvo el panel a. La figura completa (`figure1.pdf`) deja un espacio vacío para el panel a, que se inserta al armar la versión final.

### Paneles

**a. Esquema del estudio (hecho a mano, `1a.pdf`).** Un cerebro con las cinco regiones marcadas. De cada región sale un par de redes de ejemplo, CT → AD, con la J entre ellas. Le muestra al lector, antes de cualquier número, qué es una red de co-expresión, qué se compara (CT contra AD dentro de cada región) y qué mide la J.

> Al modificarlo, para que coincida con el diseño: la versión actual muestra el hipocampo (HCN), que no está en el diseño, y le falta la corteza cingulada anterior (ACC). Las etiquetas deben ser las del resto de las figuras (CT en vez de CTL, TCX en vez de TC, CBE en vez de CRB). Los valores de J deben ser los resultados reales del análisis 1, o indicarse como ilustrativos.

**b. Curva de J de referencia, J_ref(n).** Eje X: n (personas por red). Eje Y: J entre dos redes hechas con dos mitades del mismo grupo, sin ninguna persona en común. Puntos negros: calibración con DLPFC AD (n = 60, 100 y 140), con las 5 réplicas como puntos chicos y su promedio como punto grande. Línea: el modelo J_ref(n) = a + b · log(n). Círculos blancos: validación con otros tres grupos (DLPFC CT con n = 70, ACC AD con n = 85 y PCC AD con n = 75). **Cómo se lee:** J sube con n aunque no haya ninguna diferencia biológica; si los círculos blancos caen sobre la línea, la curva vale para cualquier grupo y solo depende del n.

**c. Índice de rewiring por región (panel principal).** Eje Y: J(AD, CT) / J_ref(n). Un punto por réplica AD (10 por región), con el promedio y su desviación estándar en negro. Las regiones se agrupan por estudio (ROSMAP y Mayo). Línea punteada en 1: AD y CT se parecen tanto como dos redes del mismo grupo (sin rewiring). Encima de TCX y CBE, la predicción registrada con su p de Wilcoxon. **Cómo se lee:** cuanto más abajo de 1, más rewiring; las regiones se comparan por la altura de sus puntos.

**d. Estructura de las redes, AD contra CT.** Dos medidas de similitud entre cada red AD y la red CT de su región: la correlación de Spearman del grado de los genes (¿los genes muy conectados en CT lo siguen siendo en AD?) y la NMI de las comunidades de Leiden (¿los genes se agrupan igual?). Un punto por réplica, con el promedio en negro. **Cómo se lee:** valores más bajos indican que AD cambia más la estructura de la red, no solo enlaces sueltos.

**e. Densidad de las redes.** Grado medio de la red CT (un punto) y de las 10 redes AD de cada región. **Cómo se lee:** si AD y CT tienen densidades parecidas, las diferencias de los paneles c y d no se deben a que una red tenga más enlaces que la otra.

### Pie de figura propuesto

> **Figure 1. Alzheimer's disease rewires gene co-expression networks in a region-specific manner that follows regional vulnerability.**
> **a**, Study design. Co-expression networks were inferred with ARACNe-AP for control (CT) and AD donors in five brain regions from two cohorts (ROSMAP: DLPFC, ACC and PCC; Mayo: TCX and CBE). In each region, one network was built with all CT donors and ten networks with random subsamples of AD donors of the same size and sex ratio as CT. Networks are compared with the edge Jaccard index (J). **b**, Reference similarity J_ref(n): J between two networks inferred from two disjoint halves of the same group, as a function of the number of samples per network (n). Filled circles, calibration group (DLPFC AD; n = 60, 100 and 140; small dots, 5 replicates; large dots, mean); line, fitted model J_ref(n) = a + b·log(n); open circles, independent validation groups. **c**, Rewiring index, J(AD, CT) / J_ref(n), for each region, with AD and CT networks of equal size (n = 141, 94, 99, 68 and 65 donors for DLPFC, ACC, PCC, TCX and CBE). Dots, AD subsamples (10 per region); black, mean ± s.d. Dashed line, no rewiring (AD and CT as similar as two networks of the same group). The prediction registered before network inference (TCX more rewired than CBE) was tested with a two-sided Wilcoxon rank-sum test over the 10 subsamples of each region [p = …]. **d**, Similarity of network structure between AD and CT: Spearman correlation of gene degree and normalized mutual information (NMI) of Leiden communities, computed on genes with at least one edge in both networks. **e**, Mean degree of the CT network and of the 10 AD networks of each region.

## Figura 2. La red detecta cambios que la expresión diferencial no ve

**Mensaje:** el rewiring y la expresión diferencial se solapan solo en parte. Los genes recableados que no son DE tienen funciones propias, y los hubs que cambian se comparan con el DE en cada región.

Se genera con la regla `plot_networks_de`.

### Paneles

**a. Genes recableados: ¿son DE?** Una barra por región, agrupadas por estudio. Cada barra son los genes recableados de la región (el 10% con mayor rewiring ajustado), divididos en DE ("rewired and DE", rojo) y no DE ("rewired only", azul). Línea negra: la fracción DE esperada por azar entre genes con la misma expresión y el mismo grado. Encima de cada barra, el odds ratio (OR) y su q. **Cómo se lee:** todas las barras suman 100% por definición; lo que importa es si la parte roja queda por encima de la línea negra (los recableados son DE más de lo esperado, OR > 1). La parte azul son los genes que la red detecta y el DE no. La línea del esperado es necesaria porque el número de genes DE cambia mucho entre regiones.

**b. Funciones de los genes "solo recableados".** Términos de GO Biological Process sobrerrepresentados entre los genes recableados que no son DE: los 5 más significativos de cada región. Filas: términos; columnas: regiones. Tamaño del punto: número de genes; color: −log10(q). Los círculos grises vacíos marcan términos que no son significativos en esa región. **Cómo se lee:** qué procesos cambian su coordinación en AD sin cambiar su nivel de expresión, y si esos procesos se repiten entre regiones o son propios de una.

**c. Hubs que cambian: ¿son DE?** Tres facetas: hubs compartidos, ganados en AD y perdidos en AD. En cada una, por región: el porcentaje de hubs que son DE (punto lleno) contra el esperado por azar entre genes con la misma expresión (círculo vacío). Encima, el número de hubs, y un asterisco si q < 0.05. **Cómo se lee:** una línea larga entre el círculo vacío y el punto lleno indica que ese tipo de hub es DE más (o menos) de lo esperado.

### Pie de figura propuesto

> **Figure 2. Network rewiring and differential expression capture partly different genes.**
> **a**, Rewired genes of each region (top 10% by adjusted rewiring), split into differentially expressed (DE; red) and not DE (blue). Rewiring of a gene was measured as 1 − Jaccard index between its neighbors in AD and in CT, averaged over the 10 AD networks, for genes with at least 5 neighbors, and adjusted for gene degree and mean expression (residuals of a linear model). DE genes: adjusted P < 0.05 and |fold change| > 1.2 in AD vs CT. Black line, DE fraction expected by chance; OR, odds ratio of rewired × DE; q, Benjamini–Hochberg-adjusted P across regions from 1,000 permutations of the DE labels within deciles of mean expression and degree. **b**, Gene Ontology Biological Process terms over-represented among rewired, not DE genes (top 5 per region; universe, genes that passed the degree filter; q, Benjamini–Hochberg-adjusted P within region). Dot size, number of genes; color, −log10(q); open grey circles, q ≥ 0.05. **c**, Fraction of DE genes among shared hubs, hubs gained in AD and hubs lost in AD (filled), and the fraction expected by chance for genes of the same expression (open). Hubs, top 1% of genes by PageRank; AD hub, hub in at least 8 of the 10 AD networks; not hub in AD, hub in at most 2. n, number of hubs; asterisk, q < 0.05 (permutation test as in a, within expression deciles; Benjamini–Hochberg across regions and hub types).

## Figuras suplementarias

### Figura S1. Demografía de los grupos

Se genera con la regla `plot_demographics`, en `results/figures/figureS1/`. Una columna de facetas por región, agrupadas por estudio; en cada una, CT y AD.

**Paneles:**
- **a.** Número de donantes por grupo.
- **b.** Proporción de mujeres y hombres.
- **c.** Edad al morir, en tres grupos (< 80, 80–89 y 90+ años).
- **d.** Número de alelos APOE ε4 (0, 1 o 2).
- **e.** Integridad del RNA (RIN), un punto por muestra con su caja.
- **f.** Intervalo post mortem (horas), igual que e.
- **g.** Combinaciones de regiones muestreadas por donante (gráfica UpSet), por estudio: cuántos donantes tienen muestra en cada combinación de regiones.

**Para qué sirve:** muestra que AD y CT difieren en sexo y edad (por eso el submuestreo iguala la proporción de sexos), que Mayo tiene grupos más chicos que ROSMAP, y que las regiones de un estudio comparten parte de los donantes.

> **Figure S1. Demographics of the AD and control groups.**
> For each region and diagnosis: **a**, number of donors; **b**, sex; **c**, age at death (ages of 90 or more are censored as 90+); **d**, number of APOE ε4 alleles; **e**, RNA integrity number (RIN); **f**, post-mortem interval (PMI). Dots, samples; boxes, median and interquartile range. **g**, Number of donors with samples in each combination of regions, per cohort. Data from the AMP-AD RNA-seq Harmonization Study (ROSMAP and Mayo); AD and CT as defined in the harmonized metadata.

### Figura S2. Diagrama del diseño

Es `docs/figures/experimental_design`, generado con `workflow/scripts/experimental_design_plot.R` fuera del pipeline. Muestra los datos, la selección de una muestra por donante, la inferencia de las redes y los siete análisis en sus tres bloques (A: rewiring AD vs CT; B: redes y DE; C: validación biológica), con lo que entra a cada uno.

> **Figure S2. Overview of the study design.**
> Residualized RNA-seq counts from the AMP-AD Harmonization Study were reduced to one sample per donor (highest RIN). Networks were inferred with ARACNe-AP (100 bootstraps; edges consolidated at P < 0.05) for three subsampling designs: (1) CT complete vs AD subsampled to the size of CT (5 regions × 10 subsamples); (2) AD and CT both subsampled to n = 90 (3 ROSMAP regions × 10 subsamples); (3) two disjoint halves of one group, to estimate the reference similarity J_ref(n). Block A compares AD and CT networks (rewiring index); block B relates rewiring and hubs to differential expression; block C tests a prediction registered before network inference (TCX more rewired than CBE) and the overlap with AD risk genes.

### Figura S3. Robustez del índice de rewiring

Se genera con la regla `plot_network_metrics`, en `results/figures/figureS3/`.

**a. J directa AD vs CT por región**, sin dividir entre J_ref. Un punto por réplica, con el promedio y su desviación estándar; encima, la predicción registrada probada con la J directa. **Para qué sirve:** TCX y CBE tienen casi el mismo n (68 y 65), así que la predicción también se puede probar sin usar la curva J_ref, que se calibró con ROSMAP.

**b. Índice de rewiring con n = 90 en todas las redes** (análisis 2), solo en las 3 regiones de ROSMAP. Un punto por réplica. El subtítulo dice en cuántas de las 10 réplicas el orden de las regiones es el mismo que en la Figura 1c, y la τ de Kendall promedio. **Para qué sirve:** comprueba que el orden de las regiones no depende de que cada región tenga un n distinto.

> **Figure S3. The regional pattern of rewiring does not depend on the reference curve or on sample size.**
> **a**, Edge Jaccard index between AD and CT networks (J, without normalization by J_ref) for each region, with AD and CT networks of equal size. Dots, AD subsamples (10 per region); black, mean ± s.d. The registered prediction (TCX more rewired than CBE) tested on J with a two-sided Wilcoxon rank-sum test [p = …]. **b**, Rewiring index with AD and CT both subsampled to n = 90 donors in the three ROSMAP regions (10 subsamples per region). The region order matched that of Fig. 1c in [x] of 10 subsamples (mean Kendall τ = [x]).

### Figura S4. Robustez de la comparación entre rewiring y DE

Se genera con la regla `plot_networks_de`, en `results/figures/figureS4/`.

**a. Otras definiciones de "recableado".** El OR de recableados × DE cuando "recableado" es el top 5%, 10% o 20%. Una línea por región; punto lleno si q < 0.05. **Para qué sirve:** comprueba que el resultado de la Figura 2a no depende del umbral de 10% elegido.

**b. Cada réplica AD por separado.** El OR de cada una de las 10 réplicas AD (puntos), calculado con el rewiring de esa réplica sola, junto al OR del promedio de las 10, que es la prueba principal (rombo). **Para qué sirve:** comprueba que el resultado no depende de qué donantes AD salieron en el sorteo.

> **Figure S4. Robustness of the overlap between rewired and DE genes.**
> **a**, Odds ratio of rewired × DE genes when rewired genes are defined as the top 5%, 10% or 20% by adjusted rewiring. Filled points, q < 0.05 (Benjamini–Hochberg across regions within each threshold). **b**, Odds ratio computed from each AD subsample alone (dots; top 10%), and from the rewiring averaged over the 10 subsamples (diamonds; main analysis, Fig. 2a). Dashed line, no association (odds ratio = 1).

## Tablas suplementarias

| # | Contenido | Archivo |
|---|---|---|
| T1 | Genes recableados robustos por región (top 10% en ≥ 8 de 10 réplicas), ordenados por rewiring ajustado: los genes con más rewiring | `rewiring_robust.tsv` |
| T2 | Términos GO completos de los genes "solo recableados" | `rewiring_go.tsv` |
| T3 | Lista de hubs ganados y perdidos por región | `hubs_genes.tsv` |
| T4 | Genes de riesgo de AD (análisis 7): pruebas contra recableados y hubs, y tabla por gen y región | `risk_genes_tests.tsv`, `risk_genes_table.tsv` |

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

## Decisiones tomadas

- **Enriquecimiento funcional:** GO Biological Process (`config/config.yaml`, `go`).
- **Genes de riesgo:** solo los 55 genes de los 42 loci nuevos de Bellenguez et al. (2022); los 33 loci ya conocidos no se incluyen.
- **Red consenso:** se quitó del diseño, porque ningún análisis la usa.
- **Validación con redes simuladas:** prueba interna, sin figura en el paper.
- **Esquema del estudio en la Figura 1:** sí, como panel 1a, hecho a mano (`results/figures/figure1/1a.pdf`).
