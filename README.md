# Simulación del dogma central de la biología molecular

Práctica 1 de **Bioinformática** · Grado en Ciencia e Ingeniería de Datos

Simulador en Python que sigue el flujo completo de la información genética a partir de una molécula de ADN:

```
ADN ── replicación ──► ADN ── transcripción ──► ARNm ── traducción ──► proteína
```

El programa muestra por pantalla el resultado de cada etapa, de modo que se puede seguir la transformación de la información desde el ADN hasta la secuencia de aminoácidos.

## Qué simula el programa

### 1. Replicación del ADN (`replicacion.py`)

Replicación semiconservativa: cada molécula hija conserva una hebra parental y recibe una hebra nueva.

| Paso | Qué ocurre | Moléculas que intervienen |
|---|---|---|
| Apertura de la doble hélice | Se reduce la tensión y se separan las dos hebras | Topoisomerasa, helicasa |
| Cadena líder | Síntesis continua sobre el molde, a partir de un cebador de ARN | Primasa, ADN polimerasa |
| Cadena rezagada | Síntesis discontinua en **fragmentos de Okazaki**, cada uno con su cebador de ARN | Primasa, ADN polimerasa |
| Maduración | Se sustituyen los cebadores por ADN y se unen los fragmentos | ADN polimerasa, ligasa |

Cada fragmento de Okazaki abarca 15 nucleótidos (5 de cebador y 10 de ADN), salvo el último, que puede ser más corto. El programa muestra los tres primeros fragmentos y el número total.

### 2. Transcripción (`transcripcion.py`)

La **ARN polimerasa** lee la hebra molde (3'→5') y sintetiza el ARNm (5'→3') respetando la complementariedad `A→U`, `T→A`, `C→G`, `G→C`. El ARNm resultante tiene la misma secuencia que la hebra codificante, con U en lugar de T.

La transcripción se realiza sobre la hebra codificante de la molécula hija obtenida en la replicación, de modo que las etapas quedan encadenadas (ADN → ADN → ARNm).

### 3. Traducción (`traduccion.py`)

El **ribosoma** localiza el primer codón de inicio `AUG` y lee el ARNm de tres en tres bases. Para cada codón se muestra el **anticodón del ARNt** correspondiente y el aminoácido que se incorpora, hasta encontrar un codón de parada (`UAA`, `UAG` o `UGA`). Si falta la parada o la longitud no es múltiplo de 3, el programa avisa. La tabla del código genético procede de Biopython y se puede elegir con `--tabla`.

## Estructura del repositorio

```
Simulacion-del-dogma-central-de-la-biologia-molecular/
├── data/
│   ├── ejemplo.fasta                 # secuencia de demostración (formato FASTA)
│   ├── ejemplo.raw                   # la misma secuencia en formato RAW
│   ├── insulina_NM_000207.3.fasta    # caso real: insulina humana (NCBI)
│   └── insulina_NM_000207.3.raw
├── test/
│   ├── test_replicacion.py
│   ├── test_transcripcion.py
│   └── test_traduccion.py
├── .gitignore
├── .python-version
├── LICENSE
├── README.md
├── io_utils.py                       # lectura y validación de secuencias
├── main.py                           # punto de entrada (línea de comandos)
├── pyproject.toml
├── replicacion.py                    # ADN → ADN
├── simulador.ipynb                   # notebook con la teoría y la ejecución paso a paso
├── traduccion.py                     # ARNm → proteína
├── transcripcion.py                  # ADN → ARNm
└── uv.lock
```

## Requisitos e instalación

- Python (versión indicada en `.python-version` / `pyproject.toml`)
- Biopython (tabla del código genético)

Con [uv](https://docs.astral.sh/uv/) (el repositorio incluye `uv.lock`):

```bash
git clone <url-del-repositorio>
cd Simulacion-del-dogma-central-de-la-biologia-molecular
uv sync
```

Alternativa con `pip`:

```bash
python -m venv .venv
source .venv/bin/activate        # En Windows: .venv\Scripts\activate
pip install biopython ipykernel
```

## Uso

Desde la raíz del repositorio:

```bash
# Secuencia de ejemplo por defecto (data/ejemplo.fasta)
uv run python main.py

# Una secuencia concreta, en formato FASTA o RAW
uv run python main.py data/insulina_NM_000207.3.fasta
uv run python main.py data/ejemplo.raw

# Otra tabla del código genético (por ejemplo, la 2: mitocondrial de vertebrados)
uv run python main.py data/ejemplo.fasta --tabla 2
```

Si no usas `uv`, sustituye `uv run python` por `python`.

| Argumento | Descripción | Por defecto |
|---|---|---|
| `ruta` | Archivo FASTA o RAW con la secuencia de ADN | `data/ejemplo.fasta` |
| `--tabla` | Número de tabla del código genético (numeración de NCBI) | `1` (estándar) |

## Formatos de entrada

El formato se distingue por la extensión del archivo.

- **FASTA** (`.fasta`): una línea de cabecera que empieza por `>` y a continuación la secuencia, en una o varias líneas. Solo se admiten archivos con **una única secuencia**. El identificador mostrado es la cabecera completa.
- **RAW** (`.raw`): solo la secuencia, sin cabecera. El identificador es el nombre del archivo sin extensión.

La secuencia debe estar formada únicamente por `A`, `T`, `G` y `C` (se aceptan minúsculas, que se convierten a mayúsculas). En caso contrario, el programa muestra un error y termina sin simular.

Ejemplo de `data/ejemplo.fasta`:

```
>gen_ejemplo secuencia codificante de demostracion
TTGACATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCTTAACG
```

## Ejemplo de salida

Resultado de `python main.py` con la secuencia de ejemplo:

```text
Secuencia 'gen_ejemplo secuencia codificante de demostracion': 51 nucleótidos

=== REPLICACIÓN DEL ADN ===

1. ADN parental:
  5' TTGACATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCTTAACG 3'
  3' AACTGTACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATCGAATTGC 5'

2. Apertura de la doble hélice [Topoisomerasa, Helicasa]:
  La topoisomerasa reduce la tensión y la helicasa separa las dos hebras.
  Hebra parental 1: 5' TTGACATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCTTAACG 3'
  Hebra parental 2: 3' AACTGTACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATCGAATTGC 5'

3. Cadena líder [Primasa, ADN polimerasa]:
  Cebador: uugac
  Nueva hebra: 5' uugacATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCTTAACG 3'
  Síntesis continua por ADN polimerasa.

4. Cadena rezagada [Primasa, ADN polimerasa]:
  Okazaki 1 [1-15]: cebador guaac
  Okazaki 2 [16-30]: cebador gacuu
  Okazaki 3 [31-45]: cebador uaucg
  Total de fragmentos: 4
  Síntesis discontinua por ADN polimerasa.
  La ADN polimerasa elimina los cebadores y los sustituye por ADN.
  La ligasa une los fragmentos de Okazaki.

5. Moléculas hijas:
  Hija 1:
    5' TTGACATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCTTAACG 3'  [nueva, líder]
    3' AACTGTACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATCGAATTGC 5'  [parental]
  Hija 2:
    5' TTGACATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCTTAACG 3'  [parental]
    3' AACTGTACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATCGAATTGC 5'  [nueva, rezagada]

=== TRANSCRIPCIÓN DE ADN A ARNm [ARN polimerasa] ===
  La ARN polimerasa lee la hebra molde (3'->5') y sintetiza el ARNm (5'->3').
  Complementariedad molde -> ARNm:  A->U  T->A  C->G  G->C

  ADN Codificante: 5' TTGACATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCTTAACG 3'
  ADN Molde:       3' AACTGTACCGGTAACATTACCCGGCGACTTTCCCACGGGCTATCGAATTGC 5'
  ARNm:            5' UUGACAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGCUUAACG 3'

=== TRADUCCIÓN DE ARNm A PROTEÍNA [Ribosoma, ARNt] ===
  Inicio: AUG en la posición 6

  Codones y anticodones (anticodón en 3'->5'):
  1. Codón: AUG  Anticodón: UAC  Aminoácido: Met (M)
  2. Codón: GCC  Anticodón: CGG  Aminoácido: Ala (A)
  3. Codón: AUU  Anticodón: UAA  Aminoácido: Ile (I)
  4. Codón: GUA  Anticodón: CAU  Aminoácido: Val (V)
  5. Codón: AUG  Anticodón: UAC  Aminoácido: Met (M)
  6. Codón: GGC  Anticodón: CCG  Aminoácido: Gly (G)
  7. Codón: CGC  Anticodón: GCG  Aminoácido: Arg (R)
  STOP: UGA en la posición 27

  Proteína (7 aa): MAIVMGR
  En tres letras: Met-Ala-Ile-Val-Met-Gly-Arg
```

Con secuencias largas (como la insulina), por pantalla solo se dibujan los primeros 60 nucleótidos y las primeras filas de codones, pero **todos los cálculos usan la secuencia completa**.

## Tests

Los tests usan `unittest` y comprueban cada módulo con valores calculados a mano y contra Biopython. Desde la raíz del repositorio:

```bash
uv run python -m unittest discover -s test
```

## Notebook explicativo

`simulador.ipynb` explica la teoría de cada etapa y ejecuta el código real del repositorio, tanto con la secuencia de ejemplo como con la insulina humana (`NM_000207.3`, descargada de NCBI). Hay que abrirlo con la raíz del repositorio como directorio de trabajo, porque lee los archivos de `data/` con rutas relativas.

> `NM_000207.3` es un ARNm y en el archivo aparece escrito con T. En esta práctica se toma como representación de la hebra codificante de ADN, que tiene la misma secuencia que el ARNm.

## Licencia
This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for more details.
