"""Simulación de la traducción de ARNm a proteína."""

import textwrap

from Bio.Data import CodonTable
from Bio.SeqUtils import seq3

COMPLEMENTO_ARN = {"A": "U", "U": "A", "C": "G", "G": "C"}
FILAS = 12
ANCHO = 60


def obtener_tabla(tabla_id):
    """Devuelve la tabla del código genético (ARN) de Biopython."""
    return CodonTable.unambiguous_rna_by_id[tabla_id]

def anticodon(codon):
    """Devuelve el anticodón del ARNt, escrito 3'->5'."""
    return "".join(COMPLEMENTO_ARN[base] for base in codon)


def buscar_marco(arnm, tabla):
    """Lee el ARNm desde el primer AUG hasta el codón de parada.

    Devuelve (inicio, codones, parada). El inicio es -1 si no hay AUG y la
    parada es None si el ARNm acaba sin codón de parada. Las posiciones
    se cuentan desde 0.
    """
    inicio = arnm.find("AUG")
    if inicio == -1:
        return -1, [], None
    codones = []
    for i in range(inicio, len(arnm) - 2, 3):
        codon = arnm[i:i + 3]
        if codon in tabla.stop_codons:
            return inicio, codones, i
        codones.append(codon)
    return inicio, codones, None


def traducir_arnm(arnm, tabla_id=1):
    """Traduce el ARNm (5'->3') y devuelve la proteína en una letra.

    Devuelve "" si no hay codón de inicio AUG.
    """
    tabla = obtener_tabla(tabla_id)
    _, codones, _ = buscar_marco(arnm, tabla)
    return "".join(tabla.forward_table[codon] for codon in codones)


def es_cds(arnm, tabla):
    """Comprueba si el ARNm es una CDS completa.

    Una CDS empieza en AUG, tiene un número entero de codones, acaba en un
    codón de parada y no tiene paradas intermedias.
    """
    if len(arnm) < 6 or len(arnm) % 3 != 0:
        return False
    codones = [arnm[i:i + 3] for i in range(0, len(arnm), 3)]
    paradas = set(tabla.stop_codons)
    return (
        codones[0] == "AUG"
        and codones[-1] in paradas
        and not paradas & set(codones[:-1])
    )


def mostrar_elongacion(codones, inicio, tabla):
    """Muestra una fila por codón leído (solo las primeras FILAS)."""
    print("\nElongación (un aminoácido por codón):")
    print(f"  {'N':>3}  {'Pos':>4}  {'Codón':<6} {'Anticodón 3-5':<15} "
          f"{'AA':<4} 1L")
    for numero, codon in enumerate(codones[:FILAS], 1):
        letra = tabla.forward_table[codon]
        posicion = inicio + 3 * (numero - 1) + 1
        print(f"  {numero:>3}  {posicion:>4}  {codon:<6} "
              f"{anticodon(codon):<15} {seq3(letra):<4} {letra}")
    if len(codones) > FILAS:
        print(f"  ... y {len(codones) - FILAS} más")

def mostrar_terminacion(arnm, inicio, parada):
    """Muestra el codón de parada o los avisos si no lo hay."""
    print("\nTerminación [Factor de liberación]:")
    if parada is not None:
        codon = arnm[parada:parada + 3]
        print(f"  Codón de parada {codon} en la posición {parada + 1}")
        return
    print("  AVISO: el ARNm acaba sin codón de parada.")
    sobran = (len(arnm) - inicio) % 3
    if sobran:
        print(f"  AVISO: desde el AUG sobran {sobran} nt (no múltiplo de 3).")


def mostrar_proteina(proteina):
    """Muestra la proteína en una letra y en tres letras."""
    print(f"\nProteína ({len(proteina)} aa):")
    for linea in textwrap.wrap(proteina, ANCHO):
        print("  " + linea)
    print("En tres letras:")
    for linea in textwrap.wrap("-".join(seq3(aa) for aa in proteina), ANCHO):
        print("  " + linea)


def mostrar_cds(arnm, inicio, parada, tabla):
    """Indica si desde el AUG hasta la parada hay una CDS completa."""
    if parada is None:
        print("\nCDS completa: no (falta el codón de parada)")
        return
    region = arnm[inicio:parada + 3]
    veredicto = "sí" if es_cds(region, tabla) else "no"
    print(f"\nCDS completa: {veredicto} ({len(region)} nt)")


def mostrar_traduccion(arnm, tabla_id=1):
    """Muestra el proceso de traducción y devuelve la proteína."""
    tabla = obtener_tabla(tabla_id)
    inicio, codones, parada = buscar_marco(arnm, tabla)
    print("\n=== TRADUCCIÓN DE ARNm A PROTEÍNA [Ribosoma, ARNt] ===\n")
    print(f"Código genético: tabla {tabla_id} ({tabla.names[0]})")
    if inicio == -1:
        print("\nNo hay codón de inicio AUG: no empieza la traducción.")
        return ""
    print(f"\nIniciación: primer AUG en la posición {inicio + 1} "
          f"({inicio} nt antes sin traducir)")
    mostrar_elongacion(codones, inicio, tabla)
    mostrar_terminacion(arnm, inicio, parada)
    proteina = traducir_arnm(arnm, tabla_id)
    mostrar_proteina(proteina)
    mostrar_cds(arnm, inicio, parada, tabla)
    return proteina