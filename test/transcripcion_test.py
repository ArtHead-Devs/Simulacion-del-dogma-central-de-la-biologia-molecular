import os
import sys
import unittest

from Bio.Seq import Seq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # para encontrar transcripcion.py
import transcripcion

EJEMPLO = "TTGACATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCTTAACG"


class TestTranscripcion(unittest.TestCase):

    def test_valor_calculado_a_mano(self):
        self.assertEqual(transcripcion.transcribir_adn("ATGGCCTAA"), "AUGGCCUAA")

    def test_coincide_con_biopython(self):
        self.assertEqual(transcripcion.transcribir_adn(EJEMPLO), str(Seq(EJEMPLO).transcribe()))

    def test_el_arn_no_tiene_timina(self):
        self.assertNotIn("T", transcripcion.transcribir_adn(EJEMPLO))


if __name__ == "__main__":
    unittest.main()