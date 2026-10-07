"""Lectura y validación de secuencias de ADN en formato FASTA o RAW."""

import os

from pathlib import Path

BASES_VALIDAS = set("ATGC")


def leer_fasta(ruta: str):
    """
    Lee la primera secuencia de un archivo FASTA.

    Args:
        - ruta (str): Ruta al archivo FASTA.

    Returns:
        - (str, str): Tupla con el ID de la secuencia y la secuencia en mayúsculas, o (None, None) si no se pudo leer.
    """
    ruta = Path(ruta)

    if ruta.suffix.lower() != ".fasta":
        print(f"Error: '{ruta}' no es un archivo FASTA.")
        return None, None

    with ruta.open("r") as archivo:
        seq_id = None
        secuencia = ""

        for line in archivo:
            line = line.strip()

            if line.startswith(">"):
                if seq_id is None:
                    seq_id = line[1:]
            else:
                secuencia += line.upper()

    if seq_id is None or not secuencia:
        return None, None

    return seq_id, secuencia


def leer_raw(ruta: str):
    """
    Lee una secuencia en texto plano, sin cabecera.

    Args:
        - ruta (str): Ruta al archivo RAW.

    Returns:
        - (str, str): Tupla con el nombre del archivo (sin extensión) y la secuencia en mayúsculas.
    """
    ruta = Path(ruta)

    with ruta.open("r") as archivo:
        secuencia = "".join(archivo.read().split()).upper()

    nombre = ruta.stem

    return nombre, secuencia


def validar_adn(secuencia: str):
    """
    Comprueba que la secuencia no esté vacía y solo tenga A, T, G, C.

    Args:
        - secuencia (str): Secuencia de ADN a validar.

    Returns:
        - bool: True si la secuencia es válida, False en caso contrario.
    """
    if not secuencia:
        return False

    return set(secuencia) <= BASES_VALIDAS


def cargar_secuencia(ruta: str):
    """
    Lee un archivo FASTA o RAW según su extensión y valida el ADN.

    Args:
        - ruta (str): Ruta al archivo de secuencia.

    Returns:
        - (str, str): Tupla con el ID de la secuencia y la secuencia en mayúsculas, o (None, None) si hubo un error.
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