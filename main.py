"""Simulador del dogma central de la biología molecular."""

import argparse
from pathlib import Path

from io_utils import cargar_secuencia
from replicacion import mostrar_replicacion
from transcripcion import mostrar_transcripcion
from traduccion import mostrar_traduccion

RUTA_POR_DEFECTO = str(Path(__file__).parent / "data" / "ejemplo.fasta")


def simular(ruta: str, tabla_id: int = 1):
    """
    Ejecuta ADN -> ADN, ADN -> ARN y ARN -> proteína.


    Args:
        - ruta (str): Ruta al archivo FASTA o RAW con la secuencia de ADN.
        - tabla_id (int): Número de tabla del código genético de Biopython (1 = estándar).

    Returns:
        - str | None: Proteína con una letra por aminoácido, o None si la secuencia no se pudo cargar.
    """
    seq_id, codificante = cargar_secuencia(ruta)
    if codificante is None:
        return None

    print(f"Secuencia '{seq_id}': {len(codificante)} nucleótidos")
    hija_1, hija_2 = mostrar_replicacion(codificante)
    arnm = mostrar_transcripcion(hija_1[0])
    proteina = mostrar_traduccion(arnm, tabla_id)
    return proteina


def main():
    """
    Lee los argumentos de la línea de comandos y lanza la simulación.

    Args:
        - ruta (opcional): Archivo FASTA o RAW con la secuencia de ADN. Por defecto, data/ejemplo.fasta.
        - --tabla (opcional): Número de tabla del código genético de Biopython. Por defecto, 1 (estándar).
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ruta", nargs="?", default=RUTA_POR_DEFECTO,
                        help="archivo FASTA o RAW con la secuencia de ADN")
    parser.add_argument("--tabla", type=int, default=1,
                        help="tabla del código genético de Biopython")
    argumentos = parser.parse_args()
    simular(argumentos.ruta, argumentos.tabla)


if __name__ == "__main__":
    main()