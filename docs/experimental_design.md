# Diseño experimental

## Pregunta

¿Cambia la conectividad de las redes de co-expresión entre AD y controles (CT) de forma distinta en cada región del cerebro? ¿Ese patrón de rewiring por región se sostiene después de controlar el tamaño de muestra, la composición por sexo, el ruido de estimación y la identidad de los donantes?

## Datos

- **Fuente:** AMP-AD RNA-seq Harmonization Study (Synapse, `syn26720676`).
- **Expresión:** conteos residualizados `(diagnosis)` de `sageseqr`. Están corregidos por sexo, edad, PMI, RIN, métricas de secuenciación, batch e individuo, y conservan el efecto de diagnóstico × región. El modelo exacto está en `resources/metadata/*_Configuration_file.yml`.
- **Definición de AD/CT:** la armonizada por Sage (Braak, CERAD y diagnóstico cognitivo). Se descartan OTHER (ROSMAP), y PSP y PATH_AGE (Mayo).
- **MSBB no se usa:** sus archivos residualizados no están disponibles por correcciones en la metadata.

### Unidad de análisis: una muestra por donante

Hay 8 donantes con dos muestras en la misma región (7 en DLPFC y 1 en PCC), secuenciadas en batches distintos. Dejarlas a ambas contaría dos veces a la misma persona en una red. **De cada donante se conserva la muestra con mayor RIN**; si hay empate, la de `specimenID` que va primero alfabéticamente.

### Tamaños de grupo (individuos)

| Estudio | Región | AD | CT | ¿≥100 en ambos? |
|---|---|---|---|---|
| ROSMAP | DLPFC | 289 | 141 | sí |
| ROSMAP | ACC | 171 | 94 | casi |
| ROSMAP | PCC | 150 | 99 | casi |
| Mayo | TCX | 80 | 68 | no |
| Mayo | CBE | 79 | 65 | no |

## Qué es una red en este diseño

Cada red es **una corrida completa de ARACNe-AP sobre un grupo de personas**:

1. Recibe la matriz de expresión de esas personas.
2. Hace bootstraps: re-sortea con reemplazo esas mismas personas y arma una red en cada vuelta.
3. Consolida: conserva los enlaces que aparecen en suficientes bootstraps, con **p < 0.05**.
4. Entrega **una lista de enlaces**: esa es la red que se usa para J, ρ y NMI.

El error se controla en **dos niveles**:

| Nivel | Quién lo hace | Qué controla |
|---|---|---|
| Dentro de una red | ARACNe-AP (bootstraps + consolidación) | Que los enlaces sean estables **para esas personas** |
| Entre redes | Submuestreo + réplicas (este diseño) | Que el resultado no dependa de **cuáles personas** se eligieron, del n, del sexo ni de los donantes |

El bootstrap de ARACNe solo re-sortea a las personas que recibió, así que no puede estimar qué pasaría con otro grupo de personas. Para eso están las réplicas de este diseño. Como cada red ya sale estabilizada por ARACNe, **bastan 10 réplicas**.

**Condición para todo el diseño:** todas las redes usan los mismos parámetros de ARACNe-AP (mismo número de bootstraps, misma lista de factores de transcripción, mismos genes, mismo umbral de MI y consolidación con p < 0.05). Si algún parámetro cambia, cambia la densidad de la red y la J deja de ser comparable.

## Principios del diseño

1. **Igualar el n entre las redes que se comparan.** Con más muestras, ARACNe tiene más poder y detecta más enlaces. Si una red tiene el doble de muestras que otra, tendrá más enlaces, y la J bajará aunque no exista ningún rewiring biológico. Como en todas las comparaciones las dos redes tienen el mismo n, se puede usar directamente la red consolidada con p < 0.05, sin recortar a top-k.
2. **Igualar la proporción de sexos.** En ROSMAP, el grupo AD tiene más mujeres (~70%) que el CT (~55%). Al submuestrear, las mujeres y los hombres se sortean por separado para que **ambos fenotipos tengan la proporción de mujeres de CT**. Cuando falta un sexo, el otro completa el n y el log lo reporta. No se balancea 50/50 porque eso bajaría el n.
3. **Réplicas.** Cada sorteo se repite **10 veces**. Así se ve si un resultado depende de qué personas salieron en el sorteo.
4. **Reproducibilidad.** Las semillas son fijas, se definen en `config/config.yaml` y son distintas para cada archivo. Con la misma configuración se obtienen exactamente las mismas submuestras.
5. **Las comparaciones principales son dentro de cada estudio.** ROSMAP y Mayo difieren en cohorte, protocolo y calidad del RNA: Mayo tiene RIN mucho más alto. Esas diferencias no se pueden separar de la región.

## Análisis

### 1. AD vs CT con el mismo n (análisis principal) · `matched`

**Pregunta:** ¿cuánto rewiring hay entre AD y CT en cada región?

1. **CT** se usa completo y da 1 red (`CT_all`).
2. **AD** se submuestrea al n de CT, con la proporción de sexos de CT, y da 1 red por réplica (`AD_rep01` … `AD_rep10`).
3. **Comparación:** cada red AD contra la red CT. Son 10 valores de J (y de ρ y NMI) por región.

CT no tiene réplicas porque se usa completo: no hay otras personas CT que elegir. Su estabilidad la garantiza ARACNe.

| Región | CT (completo) | AD (submuestreado) |
|---|---|---|
| DLPFC | 141 | 141 de 289 |
| ACC | 94 | 94 de 171 |
| PCC | 99 | 99 de 150 |
| TCX | 68 | 68 de 80 |
| CBE | 65 | 65 de 79 |

**Redes:** 5 CT + 5 × 10 AD = **55**.

**Resultado:** J promedio AD vs CT por región con su dispersión, y de ahí el **ranking de regiones**.

**Redes completas.** No se construyen como análisis aparte. Las redes CT de este diseño ya son las completas. Para describir la red AD (hubs, módulos) se usa una **red consenso**: los enlaces presentes en la mayoría de las réplicas AD.

### 2. ¿El ranking de regiones depende del n? (robustez) · `equal_n` → `rosmap_n90`

**Pregunta:** en el análisis 1, DLPFC tiene n = 141 y ACC n = 94. ¿El orden de las regiones cambia si todas tienen el mismo n?

1. En cada región de ROSMAP se sortean **90 AD y 90 CT**, con la proporción de sexos de CT. Así salen `AD_rep01` y `CT_rep01`, etc. Aquí también CT se sortea, porque hay que bajarlo a 90.
2. **Comparación:** en cada réplica k, `AD_repk` contra `CT_repk`, y se ordenan las 3 regiones.

Se usa n = 90, justo por debajo de ACC CT (94), para que las réplicas varíen y para quedar cerca del mínimo de ARACNe.

**Redes:** 3 regiones × 2 × 10 = **60**.

**Resultado:** el % de réplicas que conservan el orden de regiones del análisis 1, y la τ de Kendall contra ese orden.

**Descartado: `all_n60`.** Era la misma prueba con las 5 regiones y n = 60. Se quitó por tres razones: n = 60 está por debajo de lo que pide ARACNe; la comparación entre ROSMAP y Mayo está confundida por cohorte y calidad del RNA; y la dependencia del n ya la controlan `rosmap_n90` y el índice de rewiring.

### 3. Piso de reproducibilidad: curva J(n) · `curve`

**Pregunta:** si se comparan dos redes del **mismo** grupo, hechas con personas distintas, ¿qué J se obtiene según el n? Ese es el piso: la J que se obtiene sin ninguna diferencia biológica. Sin el piso, una J AD vs CT de 0.087 no se puede interpretar.

**Por qué el bootstrap de ARACNe no da este piso:** dos bootstraps de ARACNe comparten la mayoría de las personas. AD y CT, en cambio, **nunca comparten personas**. El piso tiene que medirse igual: con dos redes sin ninguna persona en común.

**Calibración con DLPFC AD (289):**
1. Para n = 100, se sortean 200 personas y se parten en **A (100) y B (100)**, sin nadie en ambas.
2. ARACNe con A da la red A, y ARACNe con B da la red B. Se calcula J(A, B).
3. Se repite 10 veces y se promedia. Esa es la **J esperada con n = 100**.
4. Se hace igual con **n = 60 y n = 140**. Los tres puntos forman la curva J_esperada(n), que cubre todo el rango de n del diseño (65–141).

**Validación, un punto por grupo:** DLPFC CT con n = 70, ACC AD con n = 85 y PCC AD con n = 75. Si caen sobre la curva, se sostiene el supuesto de que la reproducibilidad depende sobre todo del n y no de la región o el fenotipo.

**Redes:** 6 tamaños × 2 (A y B) × 10 = **120**.

**Índice de rewiring** = J(AD, CT) / J_esperada(n), con el n de cada comparación:
- análisis 1: J_esperada(141), J_esperada(94), etc., según la región;
- análisis 2: J_esperada(90).

Un índice ≈ 1 significa que AD y CT se parecen tanto como dos redes del mismo fenotipo, es decir, no hay rewiring. Un índice ≪ 1 significa que hay rewiring.

### 4. Diseño pareado (sensibilidad) · `paired` — POSPUESTO

> Por ahora no se corre. Se activará si los análisis 1–3 muestran diferencias entre regiones que convenga separar del efecto de los donantes. La regla ya existe: se activa en `config/config.yaml` (`paired`).

**Pregunta:** si dos regiones tienen distinto rewiring, ¿es por el tejido o porque sus muestras vienen de personas distintas?

1. Se usan solo los donantes que tienen **ambas regiones**.
2. **CT:** todos los donantes CT compartidos, con su muestra de cada región. Da 1 red por región.
3. **AD:** se sortean al n de CT, con la proporción de sexos de CT. **Las mismas personas** dan una red en cada región, y se repite 10 veces.
4. **Comparación:** en cada réplica, el rewiring de la región A contra el de la región B, con exactamente las mismas personas.

Pares: **DLPFC–PCC** (93 CT, 147 AD) y **DLPFC–ACC** (89 CT, 163 AD). ACC–PCC (70 CT) y las tres regiones juntas (66 CT) quedan muy por debajo de lo que pide ARACNe.

**Redes:** 2 pares × 2 regiones × (1 CT + 10 AD) = **44**.

Al interpretar: dos regiones de los mismos donantes se parecen en parte por ser las mismas personas (genética, estado agónico). Esto ayuda a atribuir las diferencias al tejido, pero sube la similitud absoluta entre las dos regiones.

## Métricas de comparación entre redes

Se calculan sobre las redes consolidadas por ARACNe-AP (p < 0.05):

- **J de Jaccard** de enlaces.
- **ρ de Spearman** del grado de los genes comunes.
- **NMI** entre las particiones en comunidades (algoritmo por definir).

Control opcional: repetir la J con top-k enlaces, para descartar que las diferencias de densidad entre AD y CT dominen el resultado.

## Resumen

| Análisis | Pregunta | Redes | En el paper |
|---|---|---|---|
| 1. `matched` | ¿Cuánto rewiring hay por región? | 55 | principal |
| 2. `rosmap_n90` | ¿El ranking depende del n? | 60 | suplemento |
| 3. `curve` | ¿Qué J da el puro ruido? | 120 | principal (índice) y suplemento (curva) |
| 4. `paired` (pospuesto) | ¿Es el tejido o son las personas? | (44) | suplemento, si hace falta |
| **Total activo** | | **235** | |

## Cómputo

**Referencia:** una red con 100 bootstraps tarda **~1 día**. Con ese valor, las 235 redes activas (análisis 1–3) suman **~235 días de cómputo** si se corren una por una.

El tiempo real depende de cuántas corridas puedan ir en paralelo:

| Corridas en paralelo | Tiempo total |
|---|---|
| 1 | ~235 días |
| 4 | ~59 días |
| 10 | ~24 días |
| 20 | ~12 días |
| 50 | ~5 días |

**Orden de ejecución.** Primero lo que necesita el texto principal: `matched` y `curve`, 175 redes. Después `rosmap_n90`, 60 redes.

Notas:
- Lo que más tarda es el cálculo de información mutua, y aumenta con el número de muestras. Las redes de n = 60–90 deberían tardar menos que las de 141, así que ~1 día por red es un techo.
- En ARACNe-AP cada bootstrap es una corrida independiente: se calcula el umbral, se corren los bootstraps y se consolidan. En Snakemake, cada bootstrap puede ser un job aparte y repartirse entre los núcleos o nodos del servidor. El cómputo total no cambia, pero el tiempo de espera sí.

## Limitaciones

- **Mayo** (65–80 por grupo) está por debajo de lo recomendado para ARACNe-AP. Se presenta como replicación o exploración.
- **Las conclusiones más sólidas son dentro de ROSMAP,** y sobre todo en DLPFC, la única región con ≥100 en ambos fenotipos.
- **En Mayo, el balance de sexos es aproximado,** porque faltan hombres AD. Los logs de submuestreo reportan la proporción lograda en cada red.
- **La comparación entre ROSMAP y Mayo** mezcla la región con la cohorte y la calidad del RNA.

## Implementación

Los parámetros (réplicas, semilla, tamaños, pares) están en `config/config.yaml`, en la sección `subsampling`.

| Paso | Regla | Salida |
|---|---|---|
| Una muestra por donante | `select_samples` | `results/samples/{study}/{region}_{dx}.tsv` |
| Matrices por grupo | `split_expression` | `results/expression/{study}/{region}_{dx}.tsv` |
| 1. AD vs CT con el mismo n | `subsample_matched` | `results/subsamples/matched/{study}/{region}.tsv` |
| 2. Robustez del n | `subsample_equal_n` | `results/subsamples/equal_n/{set}/{study}/{region}.tsv` |
| 3. Curva J(n) | `subsample_curve` | `results/subsamples/curve/{study}/{region}_{dx}_n{n}.tsv` |
| 4. Pareado | `subsample_paired` | `results/subsamples/paired/{study}/{region_a}-{region_b}.tsv` |

Cada archivo de submuestras tiene una fila por muestra y la columna `network` indica a qué red pertenece: `CT_all`, `AD_rep01`, `rep01_A`, etc. Cada red se infiere después con ARACNe-AP a partir de la matriz de su grupo y esa lista de muestras.
