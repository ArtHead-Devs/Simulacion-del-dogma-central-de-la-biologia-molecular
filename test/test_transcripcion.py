import unittest

from Bio.Seq import Seq

import transcripcion

EJEMPLO = "TTGACATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGCTTAACG"


class TestTranscripcion(unittest.TestCase):

    def test_valor_calculado_a_mano(self):
        arnm = transcripcion.transcribir_adn("ATGGCCTAA")
        self.assertEqual(arnm, "AUGGCCUAA")

    def test_coincide_con_biopython(self):
        arnm = transcripcion.transcribir_adn(EJEMPLO)
        self.assertEqual(arnm, str(Seq(EJEMPLO).transcribe()))

    def test_el_arn_no_tiene_timina(self):
        self.assertNotIn("T", transcripcion.transcribir_adn(EJEMPLO))

    def test_el_arnm_es_la_codificante_con_uracilo(self):
        arnm = transcripcion.transcribir_adn(EJEMPLO)
        self.assertEqual(arnm, EJEMPLO.replace("T", "U"))


if __name__ == "__main__":
    unittest.main()