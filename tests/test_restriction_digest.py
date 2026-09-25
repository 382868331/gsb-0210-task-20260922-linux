"""独立验收测试。

期望值全部手工计算或以本文件内的独立小函数生成，
不调用 restriction_digest 的核心算法来产生期望值；
与核心的往返一致性检查仅作为补充测试（*RoundTrip*）。
"""

import math
import unittest

from restriction_digest import (
    MAX_BASES,
    DigestError,
    DigestUnsupportedError,
    Enzyme,
    EnzymeDefinitionError,
    SequenceError,
    compatible,
    digest,
    recognize_sites,
)


# 测试文件自己的独立参考实现（与核心代码无关，用于核对，不用于“喂”期望值）。
def ref_revcomp(s: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[c] for c in reversed(s))


ECORI = Enzyme("EcoRI", "GAATTC", 1, 5)     # G^AATTC / CTTAA^G，5' 突出 AATT
SMAI = Enzyme("SmaI", "CCCGGG", 3, 3)      # CCC^GGG，平端
NSII = Enzyme("NsiI-like", "ATGCAT", 5, 1)  # ATGCA^T / T^ACGTA，3' 突出 TGCA


class TestBluntEnds(unittest.TestCase):
    """手算平端：SmaI CCC^GGG 切 AAACCCGGGTT。"""

    def setUp(self):
        self.res = digest("AAACCCGGGTT", SMAI)

    def test_sites(self):
        self.assertEqual(self.res.cut_count, 1)
        site = self.res.cut_sites[0]
        self.assertEqual((site.start, site.strand,
                          site.top_boundary, site.bottom_boundary),
                         (3, "+", 6, 6))

    def test_fragments(self):
        f1, f2 = self.res.fragments
        self.assertEqual((f1.top_strand, f1.bottom_strand),
                         ("AAACCC", ref_revcomp("AAACCC")))
        self.assertEqual((f2.top_strand, f2.bottom_strand),
                         ("GGGTT", ref_revcomp("GGGTT")))
        self.assertEqual((f1.paired_length, f2.paired_length), (6, 5))
        for f in self.res.fragments:
            self.assertEqual(f.left_end.kind, "blunt")
            self.assertEqual(f.right_end.kind, "blunt")
            self.assertEqual(f.left_end.overhang, "")

    def test_blunt_ligation(self):
        f1, f2 = self.res.fragments
        self.assertTrue(compatible(f1.right_end, f2.left_end))


class TestStickyFivePrime(unittest.TestCase):
    """手算 5' 黏端：EcoRI 切 GGGAATTCCC（位点起点 2：上链切 3，下链切 7）。"""

    def setUp(self):
        self.res = digest("GGGAATTCCC", ECORI)

    def test_boundaries(self):
        site = self.res.cut_sites[0]
        # 起点 2：上链 2+1=3，下链 2+5=7
        self.assertEqual((site.top_boundary, site.bottom_boundary), (3, 7))

    def test_left_fragment(self):
        f = self.res.fragments[0]
        self.assertEqual(f.top_strand, "GGG")
        # 下链按 5'→3'，对齐区 top[0:7]=GGGAATT
        self.assertEqual(f.bottom_strand, ref_revcomp("GGGAATT"))
        self.assertEqual(f.paired_length, 3)
        self.assertEqual(f.left_end.kind, "blunt")
        self.assertEqual(f.right_end.kind, "5_prime")
        self.assertEqual(f.right_end.overhang, "AATT")
        self.assertEqual(f.right_end.overhang_length, 4)

    def test_right_fragment(self):
        f = self.res.fragments[1]
        self.assertEqual(f.top_strand, "AATTCCC")
        self.assertEqual(f.bottom_strand, ref_revcomp("CCC"))
        self.assertEqual(f.paired_length, 3)
        self.assertEqual(f.left_end.kind, "5_prime")
        self.assertEqual(f.left_end.overhang, "AATT")
        self.assertEqual(f.right_end.kind, "blunt")

    def test_compatible_join(self):
        f1, f2 = self.res.fragments
        self.assertTrue(compatible(f1.right_end, f2.left_end))
        # 黏端不能与分子外侧的平端拼接
        self.assertFalse(compatible(f1.right_end, f2.right_end))
        self.assertFalse(compatible(f1.left_end, f2.left_end))
        # 分子两端的平端彼此可配
        self.assertTrue(compatible(f1.left_end, f2.right_end))


class TestStickyThreePrime(unittest.TestCase):
    """手算 3' 黏端：ATGCA^T / T^ACGTA 切 CCATGCATGG。"""

    def setUp(self):
        self.res = digest("CCATGCATGG", NSII)

    def test_boundaries(self):
        # 起点 2：上链 2+5=7，下链 2+1=3
        site = self.res.cut_sites[0]
        self.assertEqual((site.top_boundary, site.bottom_boundary), (7, 3))

    def test_fragments(self):
        f1, f2 = self.res.fragments
        self.assertEqual(f1.top_strand, "CCATGCA")
        self.assertEqual(f1.bottom_strand, ref_revcomp("CCA"))
        self.assertEqual(f1.right_end.kind, "3_prime")
        self.assertEqual(f1.right_end.overhang, "TGCA")
        self.assertEqual(f2.top_strand, "TGG")
        self.assertEqual(f2.bottom_strand, ref_revcomp("TGCATGG"))
        self.assertEqual(f2.left_end.kind, "3_prime")
        self.assertEqual(f2.left_end.overhang, "TGCA")
        self.assertTrue(compatible(f1.right_end, f2.left_end))


class TestCircularSingleCut(unittest.TestCase):
    """单切环状：平端与黏端各一例，切后仍为一个线性化片段。"""

    def test_blunt(self):
        res = digest("TTTCCCGGGAAA", SMAI, topology="circular")
        self.assertEqual(res.cut_count, 1)
        self.assertEqual(len(res.fragments), 1)
        f = res.fragments[0]
        self.assertFalse(f.circular)
        # 从上链切口 6 起绕一圈
        self.assertEqual(f.top_start, 6)
        self.assertEqual(f.top_length, 12)
        self.assertEqual(f.top_strand, "GGGAAATTTCCC")
        self.assertEqual(f.bottom_strand, ref_revcomp("GGGAAATTTCCC"))
        self.assertEqual(f.paired_length, 12)
        self.assertEqual(f.left_end.kind, "blunt")
        self.assertEqual(f.right_end.kind, "blunt")
        # 线性化片段的左右端来自同一个物理切口，可重新环化
        self.assertTrue(compatible(f.right_end, f.left_end))

    def test_sticky_religation(self):
        res = digest("GGGAATTCCCAA", ECORI, topology="circular")
        self.assertEqual(res.cut_count, 1)
        f = res.fragments[0]
        self.assertEqual((f.top_start, f.top_length), (3, 12))
        self.assertEqual(f.top_strand, "AATTCCCAAGGG")
        self.assertEqual(f.bottom_strand, ref_revcomp("CCCAAGGGAATT"))
        self.assertEqual(f.left_end.kind, "5_prime")
        self.assertEqual(f.right_end.kind, "5_prime")
        self.assertEqual(f.left_end.overhang, "AATT")
        self.assertEqual(f.right_end.overhang, "AATT")
        self.assertTrue(compatible(f.right_end, f.left_end))


class TestCrossZeroRecognition(unittest.TestCase):
    """跨坐标 0 的识别位点（环长 12，位点在 start=10：位置 10,11,0,1,2,3）。"""

    SEQ = "ATTC" + "G" * 7 + "A"  # 12 bp；仅 start=10 处绕环为 GAATTC

    def test_sequence_design(self):
        self.assertEqual(len(self.SEQ), 12)
        # 独立核验：唯一命中起点是 10，且绕环字符确为 GAATTC
        n = len(self.SEQ)
        windows = {
            s: "".join(self.SEQ[(s + i) % n] for i in range(6))
            for s in range(n)
        }
        self.assertEqual([s for s, w in windows.items() if w == "GAATTC"], [10])

    def test_ecori_wrapping(self):
        res = digest(self.SEQ, ECORI, topology="circular")
        site = res.cut_sites[0]
        self.assertEqual(site.start, 10)
        # 上链 10+1=11；下链 10+5=15 ≡ 3 (mod 12)
        self.assertEqual((site.top_boundary, site.bottom_boundary), (11, 3))
        f = res.fragments[0]
        self.assertEqual(f.top_strand, "AATTCGGGGGGG")
        self.assertEqual(f.bottom_strand, ref_revcomp("CGGGGGGGAATT"))
        self.assertEqual(f.left_end.overhang, "AATT")
        self.assertEqual(f.right_end.overhang, "AATT")
        self.assertTrue(compatible(f.right_end, f.left_end))


class TestReverseStrandEnzyme(unittest.TestCase):
    """非回文酶的反向互补命中坐标转换。

    酶 GGAT 切 (1,2)；序列 TTTATCCAAA 在 start=3 只有反向互补串 ATCC。
    反向命中：上链 s+L-b = 3+4-2 = 5，下链 s+L-a = 3+4-1 = 6。
    """

    def setUp(self):
        self.enz = Enzyme("FakeAsym", "GGAT", 1, 2)
        self.seq = "TTTATCCAAA"

    def test_reverse_coords(self):
        res = digest(self.seq, self.enz)
        self.assertEqual(res.cut_count, 1)
        site = res.cut_sites[0]
        self.assertEqual((site.strand, site.start,
                          site.top_boundary, site.bottom_boundary),
                         ("-", 3, 5, 6))

    def test_fragments_and_ends(self):
        res = digest(self.seq, self.enz)
        f1, f2 = res.fragments
        self.assertEqual(f1.top_strand, "TTTAT")
        self.assertEqual(f1.bottom_strand, ref_revcomp("TTTATC"))
        self.assertEqual(f1.right_end.kind, "5_prime")
        self.assertEqual(f1.right_end.overhang, "G")  # revcomp(C)=G
        self.assertEqual(f2.top_strand, "CCAAA")
        self.assertEqual(f2.bottom_strand, ref_revcomp("CAAA"))
        self.assertEqual(f2.left_end.overhang, "C")
        self.assertTrue(compatible(f1.right_end, f2.left_end))

    def test_forward_hit_on_other_sequence(self):
        # 同酶正向命中 GGAT：上链 s+1，下链 s+2
        res = digest("AAGGATTTTT", self.enz)
        site = res.cut_sites[0]
        self.assertEqual((site.strand, site.top_boundary,
                          site.bottom_boundary), ("+", 3, 4))

    def test_palindrome_counted_once(self):
        # EcoRI 回文：同一位点不能因正反向各扫一遍而计成两次断裂。
        res = digest("GGGAATTCCC", ECORI)
        self.assertEqual(res.cut_count, 1)
        self.assertEqual(len(res.cut_sites), 1)


class TestRepeatedSites(unittest.TestCase):
    """重复位点：线性 GAATTCGAATTC，两个相同切口与三段。"""

    def setUp(self):
        self.res = digest("GAATTCGAATTC", ECORI)

    def test_counts(self):
        self.assertEqual(self.res.cut_count, 2)
        self.assertEqual([s.start for s in self.res.cut_sites], [0, 6])
        self.assertEqual(len(self.res.fragments), 3)

    def test_middle_fragment(self):
        _f0, f1, _f2 = self.res.fragments
        self.assertEqual(f1.top_strand, "AATTCG")
        self.assertEqual(f1.bottom_strand, ref_revcomp("CGAATT"))
        self.assertEqual(f1.left_end.overhang, "AATT")
        self.assertEqual(f1.right_end.overhang, "AATT")
        # 中间片段两端互相可配（同型 5' 突出）
        self.assertTrue(compatible(f1.left_end, f1.right_end) or
                        compatible(f1.right_end, f1.left_end))

    def test_end_fragments(self):
        f0, _f1, f2 = self.res.fragments
        self.assertEqual(f0.top_strand, "G")
        self.assertEqual(f0.bottom_strand, ref_revcomp("GAATT"))
        self.assertEqual(f2.top_strand, "AATTC")
        self.assertEqual(f2.bottom_strand, "G")

    def test_duplicate_physical_cut_dedup(self):
        # 两个酶定义相同识别串与相同偏移 -> 两个命中证据，一个物理断裂。
        other = Enzyme("EcoRI-clone", "GAATTC", 1, 5)
        res = digest("GGGAATTCCC", [ECORI, other])
        self.assertEqual(len(res.cut_sites), 2)
        self.assertEqual(res.cut_count, 1)
        self.assertEqual(len(res.fragments), 2)


class TestNoCut(unittest.TestCase):
    """不切：线性得到一个两端平端的完整片段；环状得到无末端环状片段。"""

    def test_linear_uncut(self):
        seq = "ACGTACGTAC"
        res = digest(seq, ECORI)
        self.assertEqual(res.cut_sites, ())
        self.assertEqual(len(res.fragments), 1)
        f = res.fragments[0]
        self.assertFalse(f.circular)
        self.assertEqual((f.top_strand, f.bottom_strand),
                         (seq, ref_revcomp(seq)))
        self.assertEqual(f.paired_length, 10)
        self.assertEqual(f.left_end.kind, "blunt")
        self.assertEqual(f.right_end.kind, "blunt")

    def test_circular_uncut(self):
        seq = "ACGTACGT"
        res = digest(seq, ECORI, topology="circular")
        self.assertEqual(res.cut_sites, ())
        self.assertEqual(len(res.fragments), 1)
        f = res.fragments[0]
        self.assertTrue(f.circular)
        self.assertIsNone(f.left_end)
        self.assertIsNone(f.right_end)
        self.assertEqual((f.top_strand, f.bottom_strand),
                         (seq, ref_revcomp(seq)))
        self.assertFalse(compatible(f.left_end, f.right_end))


class TestIUPAC(unittest.TestCase):
    """识别串 IUPAC 简并：RATT 匹配 AATT（R=A/G），平端 (2,2)。"""

    def setUp(self):
        self.enz = Enzyme("Deg", "RATT", 2, 2)

    def test_match_a(self):
        res = digest("CCAATTCC", self.enz)
        self.assertEqual(res.cut_count, 1)
        self.assertEqual(res.cut_sites[0].start, 2)
        f1, f2 = res.fragments
        self.assertEqual(f1.top_strand, "CCAA")
        self.assertEqual(f2.top_strand, "TTCC")

    def test_match_g(self):
        res = digest("CCGATTCC", self.enz)
        self.assertEqual(res.cut_count, 1)

    def test_no_match(self):
        res = digest("CCTATTCC", self.enz)
        self.assertEqual(res.cut_count, 0)

    def test_n_wildcard(self):
        # GANTC 为回文简并模式（正反向同构，只产生一次物理断裂）。
        enz = Enzyme("Any", "GANTC", 2, 2)
        res = digest("CCGATTCGG", enz)  # start=2: GATTC 命中 GANTC
        self.assertEqual(res.cut_count, 1)
        self.assertEqual(res.cut_sites[0].start, 2)
        self.assertEqual(res.cut_sites[0].strand, "+")
        f1, f2 = res.fragments
        self.assertEqual((f1.top_strand, f2.top_strand), ("CCGA", "TTCGG"))


class TestUnsupported(unittest.TestCase):
    """交错切口无法做简单归属时必须明确拒绝。"""

    def test_overlapping_staggered_cuts_circular(self):
        # 同一位点两个酶、上链切口相差 1 且都带 3~4nt 突出：
        # 相邻切口间距 1 容不下突出 -> 非简单片段。
        e1 = Enzyme("E1", "GAATTC", 1, 5)
        e2 = Enzyme("E2", "GAATTC", 2, 5)
        with self.assertRaises(DigestUnsupportedError):
            digest("GAATTCAAAA", [e1, e2], topology="circular")

    def test_stagger_too_large_for_circle(self):
        enz = Enzyme("Wide", "AAA", 0, 12)  # |delta|=12 >= 环长
        with self.assertRaises(DigestUnsupportedError):
            digest("AAAAAAAAAA", enz, topology="circular")

    def test_cut_outside_linear(self):
        enz = Enzyme("Outside", "GAATTC", 1, 12)  # 下链切口超出 n=10
        with self.assertRaises(DigestUnsupportedError):
            digest("GAATTCAAAA", enz)

    def test_one_strand_cut_at_molecular_end(self):
        enz = Enzyme("Edge", "GAATTC", 0, 2)  # 起点 0：上链切口 0、下链 2
        with self.assertRaises(DigestUnsupportedError):
            digest("GAATTCAAAA", enz)

    def test_same_top_cut_conflicting_partner(self):
        e1 = Enzyme("E1", "GAATTC", 1, 5)
        e2 = Enzyme("E2", "GAATTC", 1, 4)  # 同上链切口、不同下链切口
        with self.assertRaises(DigestUnsupportedError):
            digest("GAATTCAAAA", [e1, e2])

    def test_overhangs_consume_whole_interval(self):
        # E1 切 GAATTC 于 (1,5)：上链 1、下链 5。
        # E2 识别 AATTCC 切 (2,6)，在 start=1：上链 3、下链 7。
        # 两个上链切口间距 2，却各拖 4nt 5'突出，配对区 = 2-4 < 0。
        e1 = Enzyme("E1", "GAATTC", 1, 5)
        e2 = Enzyme("E2", "AATTCC", 2, 6)
        with self.assertRaises(DigestUnsupportedError):
            digest("GAATTCCAAAA", [e1, e2])


class TestValidation(unittest.TestCase):
    """可定位的输入错误；拒绝 NaN/Infinity 与 bool 整数。"""

    def test_bad_base_location(self):
        with self.assertRaises(SequenceError) as cm:
            digest("ACGTNACGT", ECORI)
        self.assertIn("位置 4", str(cm.exception))

    def test_lowercase_not_silently_fixed(self):
        with self.assertRaises(SequenceError):
            digest("gggaattccc", ECORI)

    def test_length_bounds(self):
        with self.assertRaises(SequenceError):
            digest("", ECORI)
        with self.assertRaises(SequenceError):
            digest("A" * (MAX_BASES + 1), ECORI)

    def test_recognition_bounds(self):
        with self.assertRaises(EnzymeDefinitionError):
            Enzyme("X", "A" * 13, 1, 1)
        with self.assertRaises(EnzymeDefinitionError):
            Enzyme("X", "", 0, 0)
        with self.assertRaises(EnzymeDefinitionError):
            Enzyme("X", "GAATUX", 1, 1)

    def test_reject_bool_for_ints(self):
        with self.assertRaises(EnzymeDefinitionError):
            Enzyme("X", "GAATTC", True, 5)  # type: ignore[arg-type]
        with self.assertRaises(EnzymeDefinitionError):
            Enzyme("X", "GAATTC", 1, False)  # type: ignore[arg-type]

    def test_reject_nan_infinity(self):
        for bad in (math.nan, math.inf, -math.inf):
            with self.assertRaises(EnzymeDefinitionError):
                Enzyme("X", "GAATTC", bad, 5)  # type: ignore[arg-type]

    def test_reject_plain_float(self):
        with self.assertRaises(EnzymeDefinitionError):
            Enzyme("X", "GAATTC", 1.0, 5)  # type: ignore[arg-type]

    def test_bad_topology(self):
        with self.assertRaises(ValueError):
            digest("GAATTC", ECORI, topology="circle")  # type: ignore[arg-type]

    def test_empty_enzyme_iterable(self):
        with self.assertRaises(EnzymeDefinitionError):
            digest("GAATTC", [])

    def test_non_enzyme_member(self):
        with self.assertRaises(EnzymeDefinitionError):
            digest("GAATTC", ["not an enzyme"])  # type: ignore[list-item]

    def test_bad_sequence_type(self):
        with self.assertRaises(SequenceError):
            digest(b"GAATTC", ECORI)  # type: ignore[arg-type]
        with self.assertRaises(SequenceError):
            digest(123, ECORI)  # type: ignore[arg-type]


class TestRoundTripInvariants(unittest.TestCase):
    """补充性质：不依赖手算期望值的结构不变量。"""

    def test_strand_length_balance(self):
        for seq in ("GGGAATTCCC", "GAATTCGAATTC", "AAACCCGGGTT", "CCATGCATGG"):
            res = digest(seq, [ECORI, SMAI, NSII])
            for f in res.fragments:
                # 链长差 = 两端突出所在链的带符号贡献：
                # 左端 5' 突出在上链(+)，3' 突出在下链(-)；右端相反。
                signed = 0
                total_overhang = 0
                for side, end in (("left", f.left_end), ("right", f.right_end)):
                    total_overhang += end.overhang_length
                    if end.kind == "5_prime":
                        signed += end.overhang_length if side == "left" else -end.overhang_length
                    elif end.kind == "3_prime":
                        signed += -end.overhang_length if side == "left" else end.overhang_length
                self.assertEqual(f.top_length - f.bottom_length, signed)
                # 突出核苷酸总数守恒：两链长度之和 = 2×配对区 + 两端突出
                self.assertEqual(
                    f.top_length + f.bottom_length,
                    2 * f.paired_length + total_overhang,
                )

    def test_bottom_is_5to3_complement_circular(self):
        res = digest("GGGAATTCCCAA", ECORI, topology="circular")
        f = res.fragments[0]
        # 自己独立算：下链对齐区从下链切口 7 起绕环 12 nt
        aligned = "CCCAAGGGAATT"
        self.assertEqual(f.bottom_strand, ref_revcomp(aligned))

    def test_total_top_material_preserved_circular(self):
        seq = "GGGAATTCCCAA"
        res = digest(seq, ECORI, topology="circular")
        self.assertEqual(sum(f.top_length for f in res.fragments), len(seq))

    def test_total_top_material_preserved_linear(self):
        seq = "xGAATTCyGAATTCz".replace("x", "AA").replace("y", "CC").replace("z", "TT")
        res = digest(seq, ECORI)
        self.assertEqual(sum(f.top_length for f in res.fragments), len(seq))
        # 线性分子两端始终是平端
        self.assertEqual(res.fragments[0].left_end.kind, "blunt")
        self.assertEqual(res.fragments[-1].right_end.kind, "blunt")


class TestSelectiveDigest(unittest.TestCase):
    """确定性部分酶切：识别全部位点，只启用 ID 集合中的物理位点。"""

    SEQ = "AAAAGAATTCTTTT"  # 14 bp，唯一位点 start=4

    def test_recognition_site_keeps_palindromic_direction_sources(self):
        recognized = recognize_sites(self.SEQ, ECORI)
        self.assertEqual(len(recognized.sites), 1)
        site = recognized.sites[0]
        self.assertEqual(
            (site.enzyme, site.start, site.matched_sequence, site.strand_sources),
            ("EcoRI", 4, "GAATTC", ("+", "-"))
        )

    def test_empty_enabled_set_does_not_cut_linear(self):
        recognized = recognize_sites(self.SEQ, ECORI)
        res = digest(self.SEQ, ECORI, enabled_site_ids=set())
        self.assertEqual(res.cut_sites, ())
        self.assertEqual(len(res.fragments), 1)
        self.assertEqual(res.fragments[0].top_strand, self.SEQ)
        self.assertEqual(res.recognition_sites, recognized.sites)

    def test_empty_enabled_set_preserves_circular_molecule(self):
        res = digest(self.SEQ, ECORI, topology="circular", enabled_site_ids=())
        self.assertEqual(res.cut_sites, ())
        fragment = res.fragments[0]
        self.assertTrue(fragment.circular)
        self.assertIsNone(fragment.left_end)
        self.assertIsNone(fragment.right_end)
        self.assertEqual(fragment.top_length, len(self.SEQ))

    def test_only_one_of_two_conflicting_bottom_cuts_is_enabled(self):
        # 两酶识别同一 GAATTC：上链都切在 3，下链分别切在 5 和 7。
        e1 = Enzyme("A-two", "GAATTC", 1, 3)
        e2 = Enzyme("B-four", "GAATTC", 1, 5)
        recognized = recognize_sites(self.SEQ, [e1, e2])
        site1, site2 = recognized.sites

        res_two = digest(self.SEQ, [e1, e2], enabled_site_ids=[site1.site_id])
        left, right = res_two.fragments
        self.assertEqual((res_two.cut_sites[0].top_boundary,
                          res_two.cut_sites[0].bottom_boundary), (5, 7))
        self.assertEqual((left.top_strand, left.bottom_strand),
                         ("AAAAG", ref_revcomp("AAAAGAA")))
        self.assertEqual(left.right_end.overhang, "TT")
        self.assertEqual((right.top_strand, right.bottom_strand),
                         ("AATTCTTTT", ref_revcomp("TTCTTTT")))
        self.assertEqual(right.left_end.overhang, "AA")
        self.assertEqual(sum(f.top_length for f in res_two.fragments), len(self.SEQ))
        self.assertEqual(sum(f.bottom_length for f in res_two.fragments), len(self.SEQ))
        self.assertEqual(
            sum(f.paired_length * 2 + f.left_end.overhang_length +
                f.right_end.overhang_length for f in res_two.fragments),
            2 * len(self.SEQ)
        )

        res_four = digest(self.SEQ, [e1, e2], enabled_site_ids={site2.site_id})
        left, right = res_four.fragments
        self.assertEqual((res_four.cut_sites[0].top_boundary,
                          res_four.cut_sites[0].bottom_boundary), (5, 9))
        self.assertEqual((left.top_strand, left.bottom_strand),
                         ("AAAAG", ref_revcomp("AAAAGAATT")))
        self.assertEqual(left.right_end.overhang, "AATT")
        self.assertEqual((right.top_strand, right.bottom_strand),
                         ("AATTCTTTT", ref_revcomp("CTTTT")))

        with self.assertRaises(DigestUnsupportedError):
            digest(self.SEQ, [e1, e2])

    def test_identical_dual_direction_site_cuts_once(self):
        res = digest(self.SEQ, ECORI, enabled_site_ids=[
            recognize_sites(self.SEQ, ECORI).sites[0].site_id
        ])
        self.assertEqual(len(res.cut_sites), 1)
        self.assertEqual(len(res.fragments), 2)

    def test_degenerate_site_with_different_directional_cuts_is_still_rejected(self):
        enzyme = Enzyme("AmbiguousTwo", "RN", 0, 1)
        recognized = recognize_sites("AT", enzyme)
        site0 = recognized.sites[0]
        self.assertEqual(site0.strand_sources, ("+", "-"))
        with self.assertRaises(DigestUnsupportedError):
            digest("AT", enzyme, enabled_site_ids=[site0.site_id])

    def test_foreign_or_stale_id_rejected_entirely(self):
        current = recognize_sites(self.SEQ, ECORI)
        stale = recognize_sites("TAGAATTCCT", ECORI).sites[0].site_id
        with self.assertRaises(DigestError):
            digest(self.SEQ, ECORI, enabled_site_ids=[stale])
        changed_enzyme = recognize_sites(self.SEQ, SMAI)
        self.assertNotIn(current.sites[0].site_id,
                         {site.site_id for site in changed_enzyme.sites})

    def test_circular_identity_is_rotation_invariant_and_cross_zero_cuts(self):
        seq = TestCrossZeroRecognition.SEQ
        recognized = recognize_sites(seq, ECORI, topology="circular")
        site = recognized.sites[0]
        rotated = seq[2:] + seq[:2]
        rotated_recognized = recognize_sites(rotated, ECORI, topology="circular")
        rotated_site = rotated_recognized.sites[0]
        self.assertEqual(site.site_id, rotated_site.site_id)
        self.assertEqual((site.start, rotated_site.start), (10, 8))

        res = digest(seq, ECORI, topology="circular",
                     enabled_site_ids=[site.site_id])
        self.assertEqual((res.cut_sites[0].top_boundary,
                          res.cut_sites[0].bottom_boundary), (11, 3))
        fragment = res.fragments[0]
        self.assertEqual((fragment.top_start, fragment.top_length), (11, 12))
        self.assertEqual(fragment.bottom_length, 12)
        self.assertEqual(fragment.top_strand, "AATTCGGGGGGG")
        self.assertEqual(fragment.bottom_strand, ref_revcomp("CGGGGGGGAATT"))
        self.assertEqual(fragment.paired_length, 8)
        self.assertEqual(
            fragment.left_end.overhang_length + fragment.right_end.overhang_length,
            fragment.top_length + fragment.bottom_length - 2 * fragment.paired_length
        )

    def test_default_is_identical_to_all_enabled(self):
        enzymes = [ECORI, SMAI]
        recognized = recognize_sites("AAACCCGGGTTCGAATTCAA", enzymes)
        default = digest("AAACCCGGGTTCGAATTCAA", enzymes)
        all_enabled = digest(
            "AAACCCGGGTTCGAATTCAA",
            enzymes,
            enabled_site_ids=[site.site_id for site in recognized.sites],
        )
        self.assertEqual(default.cut_sites, all_enabled.cut_sites)
        self.assertEqual(default.fragments, all_enabled.fragments)
        self.assertEqual(
            sum(f.top_length for f in default.fragments),
            sum(f.bottom_length for f in default.fragments)
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
