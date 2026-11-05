# Diseño experimental

## Pregunta

¿Cambia la conectividad de las redes de co-expresión entre AD y controles (CT) de forma distinta en cada región del cerebro? ¿Ese patrón de rewiring por región se sostiene después de controlar el tamaño de muestra, la composición por sexo y el ruido de estimación propio de inferir redes?

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

Las redes se infieren con **ARACNe-AP**, que pide unas 100 muestras por fenotipo o el mayor número posible. Por eso el diseño busca **igualar el n perdiendo el mínimo de muestras**.

## Principios del diseño

1. **Igualar el n entre las redes que se comparan.** Con más muestras, ARACNe tiene más poder y detecta más enlaces. Si una red tiene el doble de muestras que otra, tendrá más enlaces, y la J de Jaccard bajará aunque no exista ningún rewiring biológico.
2. **Igualar la proporción de sexos.** En ROSMAP, el grupo AD tiene más mujeres (~70%) que el CT (~55%). La residualización quitó el efecto lineal del sexo, pero no la diferencia en la composición de los grupos. Al submuestrear, las mujeres y los hombres se sortean por separado para que **ambos fenotipos tengan la proporción de mujeres de CT**. Cuando falta un sexo, el otro completa el n y el log lo reporta. No se balancea 50/50 porque eso bajaría el n.
3. **Réplicas.** Cada sorteo se repite R veces: 5 para probar el pipeline y 30 para el paper. Así se ve si un resultado depende de qué personas salieron en el sorteo.
4. **Reproducibilidad.** Las semillas son fijas, se definen en `config/config.yaml` y son distintas para cada archivo. Con la misma configuración se obtienen exactamente las mismas submuestras.
5. **Las comparaciones principales son dentro de cada estudio.** ROSMAP y Mayo difieren en cohorte, protocolo y calidad del RNA: Mayo tiene RIN mucho más alto. Esas diferencias no se pueden separar de la región.

## Análisis

### 1. AD vs CT con el mismo n (análisis principal) · `matched`

- En cada región, **CT se usa completo** y **AD se submuestrea al n de CT**, con la proporción de sexos de CT y R réplicas.
- Se compara cada red AD con la red CT y se obtiene la distribución de J, ρ de grado y NMI por región.

| Región | CT (completo) | AD (submuestreado) |
|---|---|---|
| DLPFC | 141 | 141 de 289 |
| ACC | 94 | 94 de 171 |
| PCC | 99 | 99 de 150 |
| TCX | 68 | 68 de 80 |
| CBE | 65 | 65 de 79 |

**Redes completas.** No se construyen como análisis aparte. Las redes CT de este diseño ya son las completas. Para describir la red AD (hubs, módulos) se usa una **red consenso**: los enlaces presentes en la mayoría de las réplicas. Es más robusta que una red AD única con más muestras que la CT.

### 2. ¿El ranking de regiones depende del n? (robustez) · `equal_n`

En el análisis 1, cada región tiene un n distinto, de 141 a 65. Para descartar que el orden de las regiones por rewiring sea un efecto del n, se repite el análisis con **el mismo n en todas las redes**. AD y CT se submuestrean con la proporción de sexos de CT y R réplicas.

- **`rosmap_n90`:** las 3 regiones de ROSMAP con n = 90, justo por debajo de ACC CT (94). Es la comparación principal entre regiones, cerca del mínimo de ARACNe.
- **`all_n60`:** las 5 regiones con n = 60, por debajo de CBE CT (65). Prueba si el ranking aguanta al n de Mayo. Está por debajo de lo recomendado para ARACNe, así que se presenta solo como robustez.

Se usa un n un poco menor que el del grupo más chico para que haya variación entre réplicas. Si CBE CT tiene 65 y se piden 65, todas las réplicas serían iguales.

**Se reporta:** el % de réplicas que conserva el orden de regiones del análisis 1 y la τ de Kendall entre cada réplica y ese orden.

### 3. Piso de reproducibilidad: curva J(n) · `curve`

Dos redes del **mismo fenotipo** inferidas con personas distintas no son idénticas, por el ruido de estimación, y con n chico se parecen todavía menos. Sin saber qué J se obtiene cuando no hay diferencia biológica, una J AD vs CT de 0.087 no se puede interpretar.

Partir cada grupo en dos mitades no sirve: las mitades quedan muy chicas para ARACNe (32–47 en la mayoría de los grupos) y además tienen otro n que la comparación. En su lugar se construye una **curva de calibración**:

1. **Calibración con DLPFC AD (289):** dos submuestras **disjuntas**, sin donantes en común, de n = 60, 80, 100 y 140, con R réplicas. Se infiere una red por submuestra y se calcula la J entre ambas. Eso da la J esperada sin rewiring para cada n.
2. **Validación:** se repite en los grupos donde caben dos submuestras disjuntas: DLPFC CT (n = 60, 70), ACC AD (n = 60, 85) y PCC AD (n = 60, 75). Si caen sobre la curva, se sostiene el supuesto de que la reproducibilidad depende sobre todo del n y no de la región o el fenotipo.
3. **Índice de rewiring** = J(AD, CT) / J_esperada(n), donde n es el de la comparación:
   - ≈ 1: AD y CT se parecen tanto como dos redes del mismo fenotipo. No hay rewiring.
   - ≪ 1: AD y CT difieren más de lo que explica el ruido. Hay rewiring.

### 4. Diseño pareado (sensibilidad) · `paired`

Muchos donantes tienen varias regiones. En ROSMAP, 163 tienen DLPFC, ACC y PCC; en Mayo, 118 tienen TCX y CBE. Si dos regiones difieren, puede ser por el tejido o porque sus muestras son de personas distintas. Este análisis usa **solo los donantes que tienen ambas regiones**, así que en cada réplica las dos regiones tienen exactamente las mismas personas.

- Pares con suficientes CT: **DLPFC–PCC** (93 CT) y **DLPFC–ACC** (89 CT). ACC–PCC (70 CT) y las tres regiones juntas (66 CT) quedan muy por debajo de lo que pide ARACNe.
- CT se usa con todos los donantes compartidos. AD se submuestrea al n de CT con la proporción de sexos de CT y R réplicas, con los mismos donantes en ambas regiones.
- Al interpretar: dos regiones de los mismos donantes se parecen en parte por ser las mismas personas (genética, estado agónico). Esto ayuda a atribuir las diferencias al tejido, pero sube la similitud absoluta entre las dos regiones.

## Métricas de comparación entre redes

- **J de Jaccard** de enlaces. Para que no dependa de la densidad, las redes se comparan con el mismo número de enlaces (top-k).
- **ρ de Spearman** del grado de los genes comunes.
- **NMI** entre las particiones en comunidades.

Los detalles (k, algoritmo de comunidades, umbrales de ARACNe-AP) están por definir.

## Número de redes (R = 30)

| Análisis | Redes aprox. |
|---|---|
| 1. `matched` | 5 CT + 150 AD |
| 2. `equal_n` | 180 (`rosmap_n90`) + 300 (`all_n60`) |
| 3. `curve` | 2 × 30 × 10 tamaños = 600 |
| 4. `paired` | 2 pares × 2 regiones × (1 CT + 30 AD) ≈ 124 |

Esto corre solo en el servidor. Primero se prueba todo el pipeline con R = 5.

## Limitaciones

- **Mayo** (65–80 por grupo) está por debajo de lo recomendado para ARACNe-AP. Se presenta como replicación o exploración.
- **Las conclusiones más sólidas son dentro de ROSMAP,** y sobre todo en DLPFC, la única región con ≥100 en ambos fenotipos.
- **En Mayo, el balance de sexos es aproximado,** porque faltan hombres AD. Los logs de submuestreo reportan la proporción lograda en cada red.
- **La comparación entre ROSMAP y Mayo** mezcla la región con la cohorte y la calidad del RNA.

## Implementación

Los parámetros (R, semilla, tamaños, pares) están en `config/config.yaml`, en la sección `subsampling`.

| Paso | Regla | Salida |
|---|---|---|
| Una muestra por donante | `select_samples` | `results/samples/{study}/{region}_{dx}.tsv` |
| Matrices por grupo | `split_expression` | `results/expression/{study}/{region}_{dx}.tsv` |
| 1. AD vs CT con el mismo n | `subsample_matched` | `results/subsamples/matched/{study}/{region}.tsv` |
| 2. Robustez del n | `subsample_equal_n` | `results/subsamples/equal_n/{set}/{study}/{region}.tsv` |
| 3. Curva J(n) | `subsample_curve` | `results/subsamples/curve/{study}/{region}_{dx}_n{n}.tsv` |
| 4. Pareado | `subsample_paired` | `results/subsamples/paired/{study}/{region_a}-{region_b}.tsv` |

Cada archivo de submuestras tiene una fila por muestra y la columna `network` indica a qué red pertenece: `CT_all`, `AD_rep01`, `rep01_A`, etc. Cada red se infiere después con ARACNe-AP a partir de la matriz de su grupo y esa lista de muestras.
