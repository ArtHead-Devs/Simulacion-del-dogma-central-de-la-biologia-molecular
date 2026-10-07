"""Simulación de la traducción de ARNm a proteína."""

from Bio.Data import CodonTable
from Bio.SeqUtils import seq3


COMPLEMENTO_ARN = {"A": "U", "U": "A", "C": "G", "G": "C"}
MAX_FILAS = 12
MAX_AA = 20


def obtener_tabla(tabla_id: int = 1):
    """
    Devuelve la tabla del código genético de Biopython (versión ARN).

    Args:
        - tabla_id (int): Número de tabla de NCBI (1 = estándar).

    Returns:
        - Bio.Data.CodonTable.NCBICodonTableRNA: Tabla con los atributos forward_table (codón -> aminoácido), stop_codons
        y start_codons.
    """
    return CodonTable.unambiguous_rna_by_id[tabla_id]


def anticodon(codon: str):
    """
    Obtiene el anticodón del ARNt que reconoce un codón.

    El anticodón es complementario y antiparalelo al codón, por lo que se escribe en sentido 3'->5' alineado con el
    codón. Por ejemplo, el anticodón de "AUG" es "UAC".

    Args:
        - codon (str): Codón de ARNm (tres bases A, U, C, G).

    Returns:
        - str: Anticodón (3'->5').
    """
    return "".join(COMPLEMENTO_ARN[base] for base in codon)


def traducir_arnm(arnm: str, tabla_id: int = 1):
    """
    Traduce el ARNm desde el primer AUG hasta un codón de parada.

    Lee codones completos de tres en tres a partir del primer AUG y se detiene en la primera parada, que no se incluye.
    Si no hay parada, traduce todos los codones completos hasta el final.

    Args:
        - arnm (str): ARNm (5'->3'). Se convierte a mayúsculas.
        - tabla_id (int): Número de tabla del código genético de Biopython (1 = estándar).

    Returns:
        - str: Proteína con una letra por aminoácido, o "" si no hay AUG.
    """
    tabla = obtener_tabla(tabla_id)
    arnm = arnm.upper()
    inicio = arnm.find("AUG")

    if inicio == -1:
        return ""

    proteina = []

    for i in range(inicio, len(arnm) - 2, 3):
        codon = arnm[i:i + 3]

        if codon in tabla.stop_codons:
            break

        proteina.append(tabla.forward_table[codon])

    return "".join(proteina)


def mostrar_traduccion(arnm: str, tabla_id: int = 1):
    """
    Muestra el proceso de traducción y devuelve la proteína.

    Imprime la posición del AUG inicial, una fila por codón con su anticodón y su aminoácido (solo los MAX_FILAS
    primeros), el codón de parada o los avisos si falta, y la proteína completa en una letra y en tres letras (solo los
    MAX_AA primeros aminoácidos).

    Args:
        - arnm (str): ARNm (5'->3'). Se convierte a mayúsculas.
        - tabla_id (int): Número de tabla del código genético de Biopython (1 = estándar).

    Returns:
        - str: Proteína con una letra por aminoácido, o "" si no hay AUG.
    """
    tabla = obtener_tabla(tabla_id)
    arnm = arnm.upper()
    proteina = traducir_arnm(arnm, tabla_id)
    inicio = arnm.find("AUG")

    print("\n=== TRADUCCIÓN DE ARNm A PROTEÍNA [Ribosoma, ARNt] ===")

    if inicio == -1:
        print("  No hay codón de inicio AUG.")
        return ""

    print(f"  Inicio: AUG en la posición {inicio + 1}")
    print("\n  Codones y anticodones (anticodón en 3'->5'):")

    for numero, aminoacido in enumerate(proteina[:MAX_FILAS], 1):
        i = inicio + 3 * (numero - 1)
        codon = arnm[i:i + 3]
        print(
            f"  {numero}. Codón: {codon}  "
            f"Anticodón: {anticodon(codon)}  "
            f"Aminoácido: {seq3(aminoacido)} ({aminoacido})")
    if len(proteina) > MAX_FILAS:
        print(f"  ... y {len(proteina) - MAX_FILAS} codones más")

    fin = inicio + 3 * len(proteina)
    parada = arnm[fin:fin + 3]
    if parada in tabla.stop_codons:
        print(f"  STOP: {parada} en la posición {fin + 1}")
    else:
        print("  AVISO: no hay codón de parada.")
        sobran = (len(arnm) - inicio) % 3
        if sobran:
            print(f"  AVISO: la longitud desde el AUG no es múltiplo de 3 (sobran {sobran} nt).")

    puntos = "-..." if len(proteina) > MAX_AA else ""
    print(f"\n  Proteína ({len(proteina)} aa): {proteina}")
    print("  En tres letras:", "-".join(seq3(aa) for aa in proteina[:MAX_AA]) + puntos)

    return proteina
