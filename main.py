"""
main.py - Punto de entrada del simulador del Dogma Central.
Uso: python main.py <archivo.fasta|archivo.raw>
"""

import sys
from io_utils import cargar_secuencia
from replicacion import replicar_adn
from transcripcion import transcribir_adn

if len(sys.argv) < 2:
    print("Uso: python main.py <archivo.fasta | archivo.raw>")
    sys.exit(1)

seq_id, secuencia = cargar_secuencia(sys.argv[1])

if secuencia is None:
    sys.exit(1)

hija1, hija2 = replicar_adn(secuencia)

arnm = transcribir_adn(secuencia)

print(f"  ID          : {seq_id}")
print(f"  Codificante : 5'-{secuencia}-3'")
print(f"  Hija 1      : 5'-{hija1}-3'")
print(f"  Hija 2      : 5'-{hija2}-3'")
print(f"  ARNm        : 5'-{arnm}-3'")
