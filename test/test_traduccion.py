import unittest

from Bio.Seq import Seq

import traduccion

ARNM = "UUGACAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGCUUAACG"
TABLA = traduccion.obtener_tabla(1)


class TestTraduccion(unittest.TestCase):

    def test_valor_calculado_a_mano(self):
        self.assertEqual(traduccion.traducir_arnm("AUGGCCUAA"), "MA")

    def test_sin_codon_de_inicio_no_hay_proteina(self):
        self.assertEqual(traduccion.traducir_arnm("UUUGGGCCC"), "")

    def test_sin_codon_de_parada_traduce_hasta_el_final(self):
        self.assertEqual(traduccion.traducir_arnm("AUGGCCAUU"), "MAI")

    def test_la_tabla_del_codigo_genetico_importa(self):
        self.assertEqual(traduccion.traducir_arnm("AUGAUAUAA", 1), "MI")
        self.assertEqual(traduccion.traducir_arnm("AUGAUAUAA", 2), "MM")

    def test_rechaza_tablas_con_codones_ambiguos(self):
        with self.assertRaises(ValueError):
            traduccion.traducir_arnm("AUGGCCUAA", 27)

    def test_coincide_con_biopython(self):
        desde_aug = ARNM[ARNM.find("AUG"):]
        util = desde_aug[:len(desde_aug) // 3 * 3]
        esperada = str(Seq(util).translate(to_stop=True))
        self.assertEqual(traduccion.traducir_arnm(ARNM), esperada)
        self.assertEqual(esperada, "MAIVMGR")

    def test_el_anticodon_es_el_complementario(self):
        self.assertEqual(traduccion.anticodon("AUG"), "UAC")

    def test_cds_valida(self):
        self.assertTrue(traduccion.es_cds("AUGGCCUAA", TABLA))

    def test_cds_invalidas(self):
        casos = ["AUGGCC", "AUGGCCUA", "GCCUAA", "AUGUAAGCCUAA", ""]
        for caso in casos:
            self.assertFalse(traduccion.es_cds(caso, TABLA), caso)


if __name__ == "__main__":
    unittest.main()