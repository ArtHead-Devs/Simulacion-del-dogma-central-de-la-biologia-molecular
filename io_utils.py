"""Lectura y validación de secuencias de ADN en formato FASTA o RAW."""

import os

from Bio import SeqIO

BASES_VALIDAS = set("ATGC")


def leer_fasta(ruta):
    """Lee la primera secuencia de un archivo FASTA.

    Devuelve (id, secuencia) o (None, None) si no hay registros.
    """
    registros = list(SeqIO.parse(ruta, "fasta"))

    if not registros:
        return None, None

    registro = registros[0]
    return registro.id, str(registro.seq).upper()


def leer_raw(ruta):
    """Lee una secuencia en texto plano, sin cabecera.

    Devuelve (nombre_del_archivo, secuencia).
    """
    with open(ruta, "r") as archivo:
        secuencia = "".join(archivo.read().split()).upper()

    nombre = os.path.splitext(os.path.basename(ruta))[0]
    return nombre, secuencia


def validar_adn(secuencia):
    """Comprueba que la secuencia no esté vacía y solo tenga A, T, G, C."""
    if not secuencia:
        return False

    return set(secuencia) <= BASES_VALIDAS


def cargar_secuencia(ruta):
    """Lee un archivo FASTA o RAW según su extensión y valida el ADN.

    Devuelve (id, secuencia) si todo es correcto, o (None, None) si no.
    """
    if not os.path.exists(ruta):
        print(f"Error: no se encuentra el archivo '{ruta}'.")
        return None, None

    if ruta.lower().endswith(".raw"):
        seq_id, secuencia = leer_raw(ruta)
    else:
        seq_id, secuencia = leer_fasta(ruta)

    if not validar_adn(secuencia):
        print(f"Error: '{ruta}' no contiene una secuencia de ADN válida.")
        return None, None

    return seq_id, secuencia