"""公开契约与错误语义测试：类型/数值边界、可定位错误、无状态性、上限边界。"""

from __future__ import annotations

import math
import time
import unittest

from re_digest import (
    MAX_RECOGNITION_LEN,
    MAX_SEQUENCE_LEN,
    DigestError,
    Enzyme,
    InvalidEnzymeError,
    InvalidSequenceError,
    UnsupportedDigestError,
    digest,
)


def _blunt4():
    return Enzyme("Hae", "GGCC", 2, 2)


class NumericContractTests(unittest.TestCase):
    def test_bool_offset_rejected(self):
        with self.assertRaises(TypeError):
            Enzyme("b", "GG", True, 0)  # type: ignore[arg-type]
        with self.assertRaises(TypeError):
            Enzyme("b", "GG", 0, False)  # type: ignore[arg-type]

    def test_float_nan_inf_rejected(self):
        for bad in (math.nan, math.inf, -math.inf, 1.0):
            with self.assertRaises(TypeError):
                Enzyme("f", "GG", bad, 0)  # type: ignore[arg-type]

    def test_bool_circular_rejected(self):
        # 只接受真正的 bool；0/1 等 int 不得隐式充当
        with self.assertRaises(TypeError):
            digest("AAGGCCAA", _blunt4(), circular=1)  # type: ignore[arg-type]

    def test_non_string_sequence_rejected(self):
        with self.assertRaises(TypeError):
            digest(["A", "A"], _blunt4())  # type: ignore[arg-type]

    def test_non_enzyme_rejected(self):
        with self.assertRaises(TypeError):
            digest("AAGGCCAA", ("GGCC", 2, 2))  # type: ignore[arg-type]


class SequenceValidationTests(unittest.TestCase):
    def test_empty_rejected(self):
        with self.assertRaises(InvalidSequenceError) as cm:
            digest("", _blunt4())
        self.assertIsNone(cm.exception.index)

    def test_lowercase_rejected_with_location(self):
        with self.assertRaises(InvalidSequenceError) as cm:
            digest("AAGGccAA", _blunt4())
        self.assertEqual(cm.exception.index, 4)

    def test_illegal_char_with_location(self):
        with self.assertRaises(InvalidSequenceError) as cm:
            digest("ACGTXYZ", _blunt4())
        self.assertEqual(cm.exception.index, 4)

    def test_length_boundary_5000_accepted_and_fast(self):
        # 长度恰好 5000，末端构造一个 GGCC 位点（切在 4998）
        seq = ("A" * (MAX_SEQUENCE_LEN - 4)) + "GGCC"
        self.assertEqual(len(seq), MAX_SEQUENCE_LEN)
        t0 = time.perf_counter()
        res = digest(seq, _blunt4())
        self.assertLess(time.perf_counter() - t0, 2.0)
        self.assertEqual(len(res.fragments), 2)

    def test_length_5001_rejected(self):
        seq = "A" * (MAX_SEQUENCE_LEN + 1)
        with self.assertRaises(InvalidSequenceError):
            digest(seq, _blunt4())

    def test_linear_min_length_one(self):
        # n=1，平端单碱基酶；线性无内部键可切 -> 证明无解
        res = digest("A", Enzyme("A", "A", 0, 0))
        self.assertEqual(len(res.fragments), 1)
        self.assertEqual(res.fragments[0].shape, "linear")
        self.assertEqual(res.positions_scanned, 1)


class EnzymeValidationTests(unittest.TestCase):
    def test_empty_recognition_rejected(self):
        with self.assertRaises(InvalidEnzymeError) as cm:
            Enzyme("x", "", 0, 0)
        self.assertEqual(cm.exception.field, "recognition")

    def test_recognition_too_long_rejected(self):
        with self.assertRaises(InvalidEnzymeError) as cm:
            Enzyme("x", "A" * (MAX_RECOGNITION_LEN + 1), 0, 0)
        self.assertEqual(cm.exception.field, "recognition")

    def test_recognition_len_12_accepted(self):
        e = Enzyme("x", "A" * MAX_RECOGNITION_LEN, 0, 0)
        self.assertEqual(len(e.recognition), MAX_RECOGNITION_LEN)

    def test_bad_iupac_char_with_field(self):
        with self.assertRaises(InvalidEnzymeError) as cm:
            Enzyme("x", "GGQCC", 2, 2)
        self.assertEqual(cm.exception.field, "recognition")
        self.assertIn("2", str(cm.exception))

    def test_huge_offset_rejected(self):
        with self.assertRaises(InvalidEnzymeError) as cm:
            Enzyme("x", "GG", 5001, 0)
        self.assertEqual(cm.exception.field, "top_offset")

    def test_stagger_property(self):
        self.assertEqual(Enzyme("a", "GG", 1, 5).stagger, 4)
        self.assertEqual(Enzyme("a", "GG", 5, 1).stagger, -4)
        self.assertEqual(Enzyme("a", "GG", 2, 2).stagger, 0)


class TopologyAndSemanticsTests(unittest.TestCase):
    def test_unsupported_is_digest_error(self):
        enz = Enzyme("G2", "G", 0, 2)
        with self.assertRaises(DigestError):
            digest("AGGGA", enz)
        with self.assertRaises(UnsupportedDigestError):
            digest("AGGGA", enz)

    def test_digest_is_stateless_repeated_calls_identical(self):
        enz = _blunt4()
        r1 = digest("AAGGCCAAGGCC", enz)
        r2 = digest("AAGGCCAAGGCC", enz)
        self.assertEqual(r1, r2)
        # Enzyme/结果均为不可变值对象
        with self.assertRaises(Exception):
            enz.recognition = "AA"  # type: ignore[misc]

    def test_proven_no_site_distinct_from_budget(self):
        # 穷举算法无预算耗尽概念：positions_scanned 覆盖全部合法起点，
        # 无命中即证明无解（而不是“没搜到”）。
        res = digest("ACGTACGTAC", _blunt4())
        self.assertEqual(res.positions_scanned, 10 - 4 + 1)
        self.assertEqual(res.cut_events, ())
        self.assertEqual(res.site_matches, ())
        self.assertEqual(len(res.fragments), 1)

    def test_ligation_compatibility_rules(self):
        eco = Enzyme("Eco", "GAATTC", 1, 5)
        kpn = Enzyme("Kpn", "GGTACC", 5, 1)
        # 两个相同 EcoRI 位点
        r = digest("AAGAATTCAAGAATTCAA", eco)
        f = r.fragments
        self.assertTrue(f[0].right_end.compatible_with(f[1].left_end))
        # 5' 与 3' 突出不互通
        rk = digest("AAGGTACCAA", kpn)
        self.assertFalse(
            f[0].right_end.compatible_with(rk.fragments[0].right_end)
        )
        # 非互补 5' 突出端不互通
        enz2 = Enzyme("Bam-like", "GGATCC", 1, 5)  # 突出 GATC
        rb = digest("AAGGATCCAA", enz2)
        self.assertFalse(
            f[0].right_end.compatible_with(rb.fragments[0].right_end)
        )

    def test_circular_single_cut_n_equals_stagger_rejected(self):
        # n=4，单一位点，h=4：取模后上下切口都到 0，打开后不存在双链区段
        enz = Enzyme("x", "GAAA", 0, 4)
        with self.assertRaises(UnsupportedDigestError):
            digest("GAAA", enz, circular=True)

    def test_circular_recognition_longer_than_molecule_no_site(self):
        # 环状 n=3 < L=6：识别窗不可能落在 6 个不同碱基上，证明无位点
        res = digest("GAA", Enzyme("Eco", "GAATTC", 1, 5), circular=True)
        self.assertEqual(res.cut_events, ())
        self.assertEqual(res.positions_scanned, 3)
        self.assertEqual(res.fragments[0].shape, "circular")


if __name__ == "__main__":
    unittest.main()
