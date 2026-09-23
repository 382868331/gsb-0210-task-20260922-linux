"""手工参考测试：所有期望值均由人工按坐标约定独立推导，不调用被测核心生成。

被测核心只产出结果，期望值（序列、长度、边界、突出端、片段数等）全部是
本文件内手写常量。
"""

from __future__ import annotations

import unittest

from re_digest import (
    Enzyme,
    digest,
)


def eco():
    # G^AATTC / CTTAA^G：上链偏移 1，下链偏移 5（相对匹配起点），5' 突出 4nt
    return Enzyme("EcoRI-like", "GAATTC", 1, 5)


def hae3():
    # GG^CC 平端
    return Enzyme("HaeIII-like", "GGCC", 2, 2)


def kpn():
    # GGTAC^C / C^CATGG：上链偏移 5，下链偏移 1，3' 突出 4nt
    return Enzyme("KpnI-like", "GGTACC", 5, 1)


def bsa():
    # 非回文类型 IIS：GGTCTC，上链偏移 1，下链偏移 5（反向命中时切口在另一侧）
    return Enzyme("BsaI-like", "GGTCTC", 1, 5)


class HandReferenceLinearTests(unittest.TestCase):
    def test_linear_sticky_5prime(self):
        # 上链 0..9: A A G A A T T C A A；位点 s=2，上切3/下切7
        result = digest("AAGAATTCAA", eco())
        self.assertEqual(result.length, 10)
        self.assertEqual(len(result.cut_events), 1)
        ev = result.cut_events[0]
        self.assertEqual((ev.top_cut, ev.bottom_cut), (3, 7))
        # 回文位点：正、反向各命中一次，但物理事件去重为 1
        self.assertEqual(len(result.site_matches), 2)
        self.assertEqual(result.positions_scanned, 5)  # n-L+1

        self.assertEqual(len(result.fragments), 2)
        f0, f1 = result.fragments

        self.assertEqual(f0.shape, "linear")
        self.assertEqual(f0.top_sequence, "AAG")
        self.assertEqual(f0.bottom_sequence_5to3, "AATTCTT")
        self.assertEqual(f0.left_end.kind, "blunt")
        self.assertEqual(f0.left_end.overhang, "")
        self.assertEqual(f0.right_end.kind, "5prime")
        self.assertEqual(f0.right_end.overhang, "AATT")
        self.assertEqual((f0.top_start, f0.top_end), (0, 3))
        self.assertEqual((f0.bottom_start, f0.bottom_end), (0, 7))

        self.assertEqual(f1.top_sequence, "AATTCAA")
        self.assertEqual(f1.bottom_sequence_5to3, "TTG")
        self.assertEqual(f1.left_end.kind, "5prime")
        self.assertEqual(f1.left_end.overhang, "AATT")
        self.assertEqual(f1.right_end.kind, "blunt")
        self.assertEqual((f1.top_start, f1.top_end), (3, 10))
        self.assertEqual((f1.bottom_start, f1.bottom_end), (7, 10))

        # 同一 EcoRI 切口两侧互补，可连接
        self.assertTrue(f0.right_end.compatible_with(f1.left_end))
        self.assertFalse(f0.right_end.compatible_with(f0.left_end))

    def test_linear_blunt(self):
        # AA GG|CC AA，切在 4
        result = digest("AAGGCCAA", hae3())
        self.assertEqual(len(result.cut_events), 1)
        self.assertEqual((result.cut_events[0].top_cut,
                          result.cut_events[0].bottom_cut), (4, 4))
        f0, f1 = result.fragments
        self.assertEqual(f0.top_sequence, "AAGG")
        self.assertEqual(f0.bottom_sequence_5to3, "CCTT")
        self.assertEqual(f1.top_sequence, "CCAA")
        self.assertEqual(f1.bottom_sequence_5to3, "TTGG")
        for f in (f0, f1):
            self.assertEqual(f.left_end.kind, "blunt")
            self.assertEqual(f.right_end.kind, "blunt")
            self.assertEqual(f.left_end.overhang, "")
        # 任意平端互连
        self.assertTrue(f0.right_end.compatible_with(f1.left_end))
        self.assertTrue(f0.left_end.compatible_with(f1.right_end))

    def test_linear_sticky_3prime(self):
        # A A G G T A C|C A A；上切7/下切3
        result = digest("AAGGTACCAA", kpn())
        ev = result.cut_events[0]
        self.assertEqual((ev.top_cut, ev.bottom_cut), (7, 3))
        f0, f1 = result.fragments
        self.assertEqual(f0.top_sequence, "AAGGTAC")
        self.assertEqual(f0.bottom_sequence_5to3, "CTT")
        self.assertEqual(f0.right_end.kind, "3prime")
        self.assertEqual(f0.right_end.overhang, "GTAC")
        self.assertEqual(f1.top_sequence, "CAA")
        self.assertEqual(f1.bottom_sequence_5to3, "TTGGTAC")
        self.assertEqual(f1.left_end.kind, "3prime")
        self.assertEqual(f1.left_end.overhang, "GTAC")
        self.assertTrue(f0.right_end.compatible_with(f1.left_end))

    def test_linear_reverse_orientation_enzyme(self):
        # 上链含 GAGACC（GGTCTC 的反向互补），s=2；切 上3/下7
        result = digest("AAGAGACCAA", bsa())
        starts = [(m.start, m.orientation) for m in result.site_matches]
        self.assertEqual(starts, [(2, "reverse")])
        ev = result.cut_events[0]
        self.assertEqual((ev.top_cut, ev.bottom_cut), (3, 7))
        f0, f1 = result.fragments
        self.assertEqual(f0.top_sequence, "AAG")
        self.assertEqual(f0.bottom_sequence_5to3, "GTCTCTT")
        self.assertEqual(f1.top_sequence, "AGACCAA")
        self.assertEqual(f1.bottom_sequence_5to3, "TTG")
        # 右末端下链 5' 突出 GTCT，左末端上链 5' 突出 AGAC，互为反向互补
        self.assertEqual(f0.right_end.overhang, "GTCT")
        self.assertEqual(f1.left_end.overhang, "AGAC")
        self.assertTrue(f0.right_end.compatible_with(f1.left_end))

    def test_linear_repeat_sites(self):
        # GAATTC GAATTC：s=0 与 s=6；每次正反向都命中 -> 4 命中 / 2 物理事件
        result = digest("GAATTCGAATTC", eco())
        self.assertEqual(len(result.site_matches), 4)
        self.assertEqual(
            [(e.top_cut, e.bottom_cut) for e in result.cut_events],
            [(1, 5), (7, 11)],
        )
        self.assertEqual(len(result.fragments), 3)
        tops = [f.top_sequence for f in result.fragments]
        bots = [f.bottom_sequence_5to3 for f in result.fragments]
        self.assertEqual(tops, ["G", "AATTCG", "AATTC"])
        self.assertEqual(bots, ["AATTC", "AATTCG", "G"])
        # 相邻切口端均可两两连接
        f0, f1, f2 = result.fragments
        self.assertTrue(f0.right_end.compatible_with(f1.left_end))
        self.assertTrue(f1.right_end.compatible_with(f2.left_end))

    def test_linear_no_site_is_proven_absence(self):
        # 穷举扫描：无位点 = 证明无解；整分子作为一个平端线性片段
        result = digest("AAAAAAAA", hae3())
        self.assertEqual(result.cut_events, ())
        self.assertEqual(result.site_matches, ())
        self.assertEqual(result.positions_scanned, 5)
        self.assertEqual(len(result.fragments), 1)
        f = result.fragments[0]
        self.assertEqual(f.shape, "linear")
        self.assertEqual(f.top_sequence, "AAAAAAAA")
        self.assertEqual(f.bottom_sequence_5to3, "TTTTTTTT")
        self.assertEqual(f.left_end.kind, "blunt")
        self.assertEqual(f.right_end.kind, "blunt")

    def test_linear_cut_at_terminal_bond_is_not_a_cut(self):
        # 偏移把切点放到分子末端键上（坐标 0 或 n）：整事件不切
        # GGTCTC 偏移 (1,7)：s=3 -> 上切4（内）、下切10（=n，末端）-> 整事件丢弃
        enz = Enzyme("outside", "GGTCTC", 1, 7)
        result = digest("AAAGGTCTCA", enz)
        self.assertEqual(result.cut_events, ())
        self.assertEqual(len(result.fragments), 1)

    def test_negative_offsets(self):
        # GG 偏移 (-1,1)：s=2 -> 上切1/下切3，5'突出 2nt
        # 5' A A | G G A A 3'
        # 3' T T C C | T T 5'
        enz = Enzyme("neg", "GG", -1, 1)
        result = digest("AAGGAA", enz)
        ev = result.cut_events[0]
        self.assertEqual((ev.top_cut, ev.bottom_cut), (1, 3))
        f0, f1 = result.fragments
        self.assertEqual(f0.top_sequence, "A")
        self.assertEqual(f0.bottom_sequence_5to3, "CTT")
        self.assertEqual(f0.right_end.kind, "5prime")
        self.assertEqual(f0.right_end.overhang, "CT")  # 下链 5'->3'
        self.assertEqual(f1.left_end.overhang, "AG")
        self.assertTrue(f0.right_end.compatible_with(f1.left_end))

    def test_overlapping_staggered_cuts_unsupported(self):
        # 酶 G 偏移 (0,2)；AGGGA 中 s=1 ->(1,3)、s=2 ->(2,4)、s=3 ->(3,5)
        # 相邻上链切口间距 1 <= 交错长度 2
        from re_digest import UnsupportedDigestError
        enz = Enzyme("G2", "G", 0, 2)
        with self.assertRaises(UnsupportedDigestError):
            digest("AGGGA", enz)


class HandReferenceCircularTests(unittest.TestCase):
    def test_circular_single_sticky_cut(self):
        # n=10，位点 s=2，单切打开环
        result = digest("AAGAATTCTT", eco(), circular=True)
        self.assertEqual(len(result.cut_events), 1)
        ev = result.cut_events[0]
        self.assertEqual((ev.top_cut, ev.bottom_cut), (3, 7))
        self.assertEqual(len(result.fragments), 1)
        f = result.fragments[0]
        self.assertEqual(f.shape, "linear")  # 切开后是线性分子
        self.assertEqual(f.top_sequence, "AATTCTTAAG")
        self.assertEqual(f.bottom_sequence_5to3, "AATTCTTAAG")
        self.assertEqual(f.left_end.kind, "5prime")
        self.assertEqual(f.right_end.kind, "5prime")
        self.assertEqual(f.left_end.overhang, "AATT")
        self.assertEqual(f.right_end.overhang, "AATT")
        self.assertTrue(f.left_end.compatible_with(f.right_end))
        self.assertEqual(result.positions_scanned, 10)  # 环状扫描 n 个起点

    def test_circular_no_site_is_uncut_circle(self):
        result = digest("AAAAAAAA", hae3(), circular=True)
        self.assertEqual(len(result.fragments), 1)
        f = result.fragments[0]
        self.assertEqual(f.shape, "circular")
        self.assertIsNone(f.left_end)
        self.assertIsNone(f.right_end)
        self.assertEqual(f.top_sequence, "AAAAAAAA")
        self.assertEqual(f.bottom_sequence_5to3, "TTTTTTTT")

    def test_circular_cross_zero_blunt_recognition_and_cut(self):
        # n=7: 位点 GGCC 跨零 s=5（G,G,C,C 位于 5,6,0,1），切在 7%7=0
        result = digest("CCAAAGG", hae3(), circular=True)
        self.assertEqual(len(result.site_matches), 2)  # 回文：正反各一
        self.assertEqual(len(result.cut_events), 1)
        self.assertEqual((result.cut_events[0].top_cut,
                          result.cut_events[0].bottom_cut), (0, 0))
        f = result.fragments[0]
        self.assertEqual(f.top_sequence, "CCAAAGG")
        self.assertEqual(f.bottom_sequence_5to3, "CCTTTGG")
        self.assertEqual(f.left_end.kind, "blunt")
        self.assertEqual(f.right_end.kind, "blunt")
        self.assertEqual(f.left_end.top_boundary, 0)
        self.assertEqual(f.right_end.top_boundary, 0)

    def test_circular_cross_zero_sticky_overhang_read(self):
        # n=11，位点 s=8 跨零（G,A,A,T,T,C 在 8,9,10,0,1,2）
        # 上切9，下切 13%11=2；左末端突出端读取需跨零
        result = digest("TTCCCCCCGAA", eco(), circular=True)
        ev = result.cut_events[0]
        self.assertEqual((ev.top_cut, ev.bottom_cut), (9, 2))
        f = result.fragments[0]
        self.assertEqual(f.top_sequence, "AATTCCCCCCG")
        self.assertEqual(f.bottom_sequence_5to3, "AATTCGGGGGG")
        self.assertEqual(f.left_end.overhang, "AATT")   # 位置 9,10,0,1
        self.assertEqual(f.right_end.overhang, "AATT")
        self.assertTrue(f.left_end.compatible_with(f.right_end))

    def test_circular_cross_zero_reverse_enzyme(self):
        # n=10，反向模式 GAGACC 跨零 s=8（G,A,G,A,C,C 在 8,9,0,1,2,3）
        result = digest("GACCTTTTGA", bsa(), circular=True)
        self.assertEqual(
            [(m.start, m.orientation) for m in result.site_matches],
            [(8, "reverse")],
        )
        self.assertEqual(
            (result.cut_events[0].top_cut, result.cut_events[0].bottom_cut),
            (9, 3),
        )
        f = result.fragments[0]
        self.assertEqual(f.left_end.overhang, "AGAC")    # 位置 9,0,1,2
        self.assertEqual(f.right_end.overhang, "GTCT")
        self.assertTrue(f.left_end.compatible_with(f.right_end))

    def test_circular_modulo_collapse_unsupported(self):
        # n=2，酶 G 偏移 (0,2)：模 2 后上下切口都到 0，交错结构被取模吞掉
        enz = Enzyme("G2", "G", 0, 2)
        from re_digest import UnsupportedDigestError
        with self.assertRaises(UnsupportedDigestError):
            digest("GC", enz, circular=True)

    def test_circular_repeat_two_events_two_fragments(self):
        # 两个 EcoRI 位点的环：s=0 与 s=5（串联），切 1/5、7/11%12=11
        result = digest("GAATTCGAATTC", eco(), circular=True)
        self.assertEqual(len(result.cut_events), 2)
        self.assertEqual(
            [(e.top_cut, e.bottom_cut) for e in result.cut_events],
            [(1, 5), (7, 11)],
        )
        f0, f1 = result.fragments
        # 上链弧长 6 与 6；下链弧分别 5..11 与 11..5（跨零），均为 6nt
        self.assertEqual(f0.top_sequence, "AATTCG")
        self.assertEqual(f1.top_sequence, "AATTCG")
        self.assertEqual(f0.bottom_sequence_5to3, "AATTCG")
        self.assertEqual(f1.bottom_sequence_5to3, "AATTCG")
        self.assertTrue(f0.left_end.compatible_with(f1.right_end))


class IupacHandTests(unittest.TestCase):
    def test_degenerate_code_match(self):
        # GNCT（N 集合；非自反向互补），s=1，平端偏移 (1,1)，切在 2
        enz = Enzyme("GNCT-blunt", "GNCT", 1, 1)
        result = digest("GGACTAAA", enz)
        self.assertEqual(
            [(m.start, m.orientation) for m in result.site_matches],
            [(1, "forward")],
        )
        self.assertEqual(len(result.cut_events), 1)
        self.assertEqual(result.cut_events[0].top_cut, 2)
        f0, f1 = result.fragments
        self.assertEqual(f0.top_sequence, "GG")
        self.assertEqual(f1.top_sequence, "ACTAAA")

    def test_odd_self_complement_degenerate_site_two_orientations_distinct(self):
        # GANTC 是长度 5 的自反向互补码：偏移 2 时正、反方向切在中央 N 的
        # 两侧（键 s+2 与 s+3），是两个不同物理平端事件，不去重。
        enz = Enzyme("GANTC-blunt", "GANTC", 2, 2)
        result = digest("GACTCAAAA", enz)
        self.assertEqual(
            sorted((m.start, m.orientation) for m in result.site_matches),
            [(0, "forward"), (0, "reverse")],
        )
        self.assertEqual(
            [(e.top_cut, e.bottom_cut) for e in result.cut_events],
            [(2, 2), (3, 3)],
        )
        self.assertEqual(
            [f.top_sequence for f in result.fragments], ["GA", "C", "TCAAAA"]
        )

    def test_purine_code_linear_terminal_bond(self):
        # R=A/G 在 AGAG 的每个位置都正向匹配；偏移 0：
        # s=0 的切点是末端键 0（整事件丢弃）；s=1,2,3 切内部键 1,2,3
        enz = Enzyme("R0", "R", 0, 0)
        result = digest("AGAG", enz)
        self.assertEqual(
            [(m.start, m.orientation) for m in result.site_matches],
            [(1, "forward"), (2, "forward"), (3, "forward")],
        )
        self.assertEqual([e.top_cut for e in result.cut_events], [1, 2, 3])

    def test_purine_code_circular_keeps_zero(self):
        enz = Enzyme("R0", "R", 0, 0)
        result = digest("AGAG", enz, circular=True)
        self.assertEqual([e.top_cut for e in result.cut_events], [0, 1, 2, 3])
        self.assertEqual(len(result.fragments), 4)


if __name__ == "__main__":
    unittest.main()
