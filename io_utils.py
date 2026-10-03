"""
io_utils.py
-----------
Módulo de entrada/salida para el simulador del Dogma Central.

Permite leer secuencias de ADN desde archivos FASTA o RAW,
y comprobar que solo contengan las bases válidas A, T, G, C.
"""

import os
from Bio import SeqIO


# Bases de ADN permitidas
BASES_VALIDAS = set("ATGC")


def leer_fasta(ruta):
    """Lee la primera secuencia de un archivo FASTA y devuelve (id, secuencia)."""
    registros = list(SeqIO.parse(ruta, "fasta"))

    if not registros:
        print(f"El archivo '{ruta}' no contiene secuencias FASTA válidas.")
        return None, None

    registro = registros[0]
    secuencia = str(registro.seq).upper()
    print(f"Secuencia FASTA leída: '{registro.id}' ({len(secuencia)} nucleótidos)")
    return registro.id, secuencia


def leer_raw(ruta):
    """Lee una secuencia de texto plano (sin cabecera) y devuelve (nombre_archivo, secuencia)."""
    with open(ruta, "r") as f:
        # Juntamos todas las líneas no vacías en una sola cadena
        secuencia = "".join(linea.strip() for linea in f if linea.strip()).upper()

    nombre = os.path.splitext(os.path.basename(ruta))[0]

    if not secuencia:
        print(f"El archivo '{ruta}' está vacío.")
        return nombre, None

    print(f"Secuencia RAW leída: '{nombre}' ({len(secuencia)} nucleótidos)")
    return nombre, secuencia


def validar_adn(secuencia):
    """
    Comprueba que la secuencia solo tenga las bases A, T, G, C.
    Devuelve True si es válida, False si contiene caracteres extraños.
    """
    bases_encontradas = set(secuencia)
    invalidas = bases_encontradas - BASES_VALIDAS

    if invalidas:
        print(f"Aviso: la secuencia tiene caracteres no esperados: {sorted(invalidas)}")
        print("Solo se aceptan las bases A, T, G, C.")
        return False

    print("Secuencia válida: solo contiene bases A, T, G, C.")
    return True


def cargar_secuencia(ruta):
    """
    Función principal. Lee el archivo (FASTA o RAW según extensión),
    valida la secuencia y la devuelve lista para usar.

    Devuelve (id, secuencia) si todo va bien, o (None, None) si hay algún problema.
    """
    print(f"\nCargando archivo: {ruta}")

    if not os.path.exists(ruta):
        print(f"Error: no se encuentra el archivo '{ruta}'.")
        return None, None

    # Elegir formato según extensión
    if ruta.endswith(".raw"):
        seq_id, secuencia = leer_raw(ruta)
    else:
        seq_id, secuencia = leer_fasta(ruta)

    if secuencia is None:
        return None, None

    # Validar que la secuencia sea ADN puro
    es_valida = validar_adn(secuencia)

    if not es_valida:
        return None, None

    print(f"Listo. Secuencia '{seq_id}' cargada correctamente.\n")
    return seq_id, secuencia


