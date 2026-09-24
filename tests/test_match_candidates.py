"""候选环状切位图筛查（match_candidates）的独立验收测试。

期望值全部手工计算，或由本文件内的独立小函数（ref_circular_lengths）
显式列窗生成，不调用 restriction_digest 的核心算法来产生期望值。
"""

import math
import unittest

from restriction_digest import (
    Enzyme,
    EnzymeDefinitionError,
    FragmentLengthsError,
    SequenceError,
    match_candidates,
)


ECORI = Enzyme("EcoRI", "GAATTC", 1, 5)   # G^AATTC / CTTAA^G
SMAI = Enzyme("SmaI", "CCCGGG", 3, 3)     # CCC^GGG，平端

# 手工构造的候选环状序列（长度、位点、切口间距均手算）：
DUAL6 = "GAATTCGAATTC"        # 12 bp；EcoRI 位点 0/6 -> 上链切口 1/7 -> (6, 6)
DUAL6_ROT = "TTCGAATTCGAA"    # DUAL6 旋转 3 位的写法（同一物理环状分子）
DUAL7A = "GAATTCCGAATTCC"     # 14 bp；位点 0/7 -> 切口 1/8 -> (7, 7)
DUAL7B = "GAATTCAGAATTCA"     # 14 bp；位点 0/7 -> (7, 7)，与 DUAL7A 非旋转等价
UNC12 = "ACGTACGTACGT"        # 12 bp；无 EcoRI 位点，未切环 -> (12,)
ONE12 = "GGGAATTCCCAA"        # 12 bp；唯一位点 3 -> 单切环 -> (12,)
TRIPLE6 = "GAATTCGAATTCGAATTC"  # 18 bp；位点 0/6/12 -> (6, 6, 6)
BOTH = "GAATTCCCCGGG"         # 12 bp；EcoRI 切口 1、SmaI 切口 9 -> (4, 8)


def ref_circular_lengths(seq: str, recognition: str, top_cut: int) -> tuple[int, ...]:
    """独立小参考：显式列窗找识别串，取上链切口间距的多重集。

    只支持具体碱基识别串；与 restriction_digest 的核心代码无关。
    """
    n = len(seq)
    length = len(recognition)
    cuts = sorted(
        (s + top_cut) % n
        for s in range(n)
        if all(seq[(s + i) % n] == recognition[i] for i in range(length))
    )
    if len(cuts) <= 1:
        return (n,)  # 未切环与单切环都是全长单片段
    return tuple(sorted(
        (cuts[(i + 1) % len(cuts)] - cuts[i]) % n for i in range(len(cuts))
    ))


class TestMatching(unittest.TestCase):
    """正常匹配：精确多重集、旋转去重、重复片段、未切环。"""

    def test_basic_match_and_rotation_dedup(self):
        res = match_candidates([DUAL6, DUAL6_ROT, DUAL7A, UNC12], ECORI, [6, 6])
        self.assertEqual(res.candidate_count, 4)
        self.assertEqual(res.class_count, 3)  # DUAL6 与其旋转写法归为一类
        self.assertEqual(res.target_lengths, (6, 6))
        self.assertEqual(res.unsupported, ())
        self.assertEqual(len(res.matches), 1)
        m = res.matches[0]
        self.assertEqual(m.representative, DUAL6)  # 输入中该类首次出现的写法
        self.assertEqual(m.indices, (0, 1))
        self.assertEqual(m.length, 12)
        self.assertEqual(m.fragment_lengths, (6, 6))
        self.assertEqual(m.cut_count, 2)

    def test_repeated_fragments_multiplicity(self):
        # 重复片段按计数匹配：(6,) 与 (6,6,6) 都不能匹配产生 (6,6) 的候选。
        candidates = [DUAL6, TRIPLE6]
        self.assertEqual(
            match_candidates(candidates, ECORI, [6]).matches, ())
        self.assertEqual(
            match_candidates(candidates, ECORI, [6, 6, 6, 6]).matches, ())
        res = match_candidates(candidates, ECORI, [6, 6, 6])
        self.assertEqual(len(res.matches), 1)
        self.assertEqual(res.matches[0].representative, TRIPLE6)
        self.assertEqual(res.matches[0].fragment_lengths, (6, 6, 6))
        self.assertEqual(res.matches[0].cut_count, 3)
        res2 = match_candidates(candidates, ECORI, [6, 6])
        self.assertEqual([m.representative for m in res2.matches], [DUAL6])

    def test_two_distinct_classes_same_pattern(self):
        # 独立核验：DUAL7B 不是 DUAL7A 的旋转。
        self.assertNotIn(DUAL7B, DUAL7A + DUAL7A)
        res = match_candidates([DUAL7A, DUAL7B], ECORI, [7, 7])
        self.assertEqual(res.class_count, 2)
        self.assertEqual([m.representative for m in res.matches],
                         [DUAL7A, DUAL7B])  # 按输入首次出现顺序
        for m in res.matches:
            self.assertEqual(m.fragment_lengths, (7, 7))
            self.assertEqual(m.cut_count, 2)

    def test_uncut_and_single_cut_same_length_pattern(self):
        # 未切环与单切环的长度多重集都是 (12,)，只能靠 cut_count 区分。
        res = match_candidates([UNC12, ONE12, DUAL6], ECORI, [12])
        self.assertEqual(len(res.matches), 2)
        uncut, single = res.matches
        self.assertEqual(uncut.representative, UNC12)
        self.assertEqual(uncut.cut_count, 0)
        self.assertEqual(uncut.fragment_lengths, (12,))
        self.assertEqual(single.representative, ONE12)
        self.assertEqual(single.cut_count, 1)
        self.assertEqual(single.fragment_lengths, (12,))

    def test_no_match(self):
        res = match_candidates([DUAL6, UNC12], ECORI, [5, 7])
        self.assertEqual(res.matches, ())
        self.assertEqual(res.candidate_count, 2)
        self.assertEqual(res.class_count, 2)
        # 总长度一致但拆分不同也不匹配。
        self.assertEqual(
            match_candidates([DUAL6, UNC12], ECORI, [4, 8]).matches, ())

    def test_target_is_multiset_order_irrelevant(self):
        res = match_candidates([BOTH], [ECORI, SMAI], [8, 4])
        self.assertEqual(res.target_lengths, (4, 8))  # 归一化为升序
        self.assertEqual(len(res.matches), 1)
        self.assertEqual(res.matches[0].fragment_lengths, (4, 8))
        self.assertEqual(res.matches[0].cut_count, 2)

    def test_multi_enzyme_map(self):
        # 手算：BOTH 上 EcoRI 上链切口 1、SmaI 平端切口 9，间距 8 与 4。
        res = match_candidates([BOTH], [ECORI, SMAI], [4, 8])
        self.assertEqual(len(res.matches), 1)
        # 只给其中一种酶则图样不同：单 EcoRI -> (12,)，单 SmaI -> (12,)。
        self.assertEqual(
            match_candidates([BOTH], ECORI, [4, 8]).matches, ())
        self.assertEqual(
            len(match_candidates([BOTH], SMAI, [12]).matches), 1)

    def test_single_enzyme_argument(self):
        res = match_candidates([DUAL6], ECORI, [6, 6])
        self.assertEqual(len(res.matches), 1)

    def test_identical_candidates_one_class(self):
        res = match_candidates([UNC12, UNC12, UNC12], ECORI, [12])
        self.assertEqual(res.class_count, 1)
        self.assertEqual(len(res.matches), 1)
        self.assertEqual(res.matches[0].indices, (0, 1, 2))

    def test_all_rotations_dedup(self):
        rotations = [DUAL6[i:] + DUAL6[:i] for i in range(len(DUAL6))]
        res = match_candidates(rotations, ECORI, [6, 6])
        self.assertEqual(res.candidate_count, 12)
        self.assertEqual(res.class_count, 1)
        self.assertEqual(len(res.matches), 1)
        self.assertEqual(res.matches[0].indices, tuple(range(12)))
        self.assertEqual(res.matches[0].representative, DUAL6)

    def test_length_one_candidate(self):
        res = match_candidates(["A"], ECORI, [1])
        self.assertEqual(len(res.matches), 1)
        self.assertEqual(res.matches[0].fragment_lengths, (1,))
        self.assertEqual(res.matches[0].cut_count, 0)

    def test_deterministic_repeat(self):
        args = ([DUAL6, DUAL7A, UNC12], ECORI, [6, 6])
        self.assertEqual(match_candidates(*args), match_candidates(*args))

    def test_generator_inputs(self):
        # 候选与酶都可以是惰性生成器：内部预先归一化，不会被多轮评估耗尽。
        cands = (c for c in (DUAL6, DUAL6_ROT, UNC12))
        enzymes = (e for e in (ECORI,))
        res = match_candidates(cands, enzymes, (x for x in (6, 6)))
        self.assertEqual(res.candidate_count, 3)
        self.assertEqual(res.class_count, 2)
        self.assertEqual(len(res.matches), 1)
        self.assertEqual(res.matches[0].indices, (0, 1))


class TestUnsupportedCandidates(unittest.TestCase):
    """几何上不可归属的候选：收入 unsupported，不中断其余候选。"""

    def test_unsupported_collected_not_raised(self):
        # E1 切 (1,5)、E2 切 (2,5)：同一物理下链切口 5 被两个上链切口
        # 争夺，digest 明确报 DigestUnsupportedError。
        e1 = Enzyme("E1", "GAATTC", 1, 5)
        e2 = Enzyme("E2", "GAATTC", 2, 5)
        res = match_candidates(["GAATTCAAAA", "CCCCCCCCCC"], [e1, e2], [10])
        self.assertEqual(len(res.unsupported), 1)
        u = res.unsupported[0]
        self.assertEqual(u.representative, "GAATTCAAAA")
        self.assertEqual(u.indices, (0,))
        self.assertTrue(u.reason)
        # 另一个候选不受影响：无位点未切环，长度多重集 (10,) 精确匹配。
        self.assertEqual(len(res.matches), 1)
        self.assertEqual(res.matches[0].representative, "CCCCCCCCCC")
        self.assertEqual(res.matches[0].cut_count, 0)

    def test_unsupported_rotation_class_grouped(self):
        e1 = Enzyme("E1", "GAATTC", 1, 5)
        e2 = Enzyme("E2", "GAATTC", 2, 5)
        bad = "GAATTCAAAA"
        bad_rot = bad[3:] + bad[:3]  # 旋转写法，同样 unsupported
        res = match_candidates([bad, bad_rot], [e1, e2], [10])
        self.assertEqual(len(res.unsupported), 1)
        self.assertEqual(res.unsupported[0].indices, (0, 1))
        self.assertEqual(res.matches, ())


class TestValidation(unittest.TestCase):
    """可定位的输入错误；任何评估之前整体校验，不存在部分结果。"""

    def test_empty_fragment_lengths(self):
        with self.assertRaises(FragmentLengthsError):
            match_candidates([DUAL6], ECORI, [])

    def test_nonpositive_lengths(self):
        for bad in ([0], [-3], [6, 0]):
            with self.assertRaises(FragmentLengthsError):
                match_candidates([DUAL6], ECORI, bad)

    def test_bool_length_rejected(self):
        with self.assertRaises(FragmentLengthsError):
            match_candidates([DUAL6], ECORI, [True])

    def test_float_lengths_rejected(self):
        for bad in ([6.0], [2.5], [math.nan], [math.inf], [-math.inf]):
            with self.assertRaises(FragmentLengthsError):
                match_candidates([DUAL6], ECORI, bad)

    def test_non_int_length(self):
        with self.assertRaises(FragmentLengthsError):
            match_candidates([DUAL6], ECORI, ["6"])  # type: ignore[list-item]

    def test_fragment_lengths_not_iterable(self):
        with self.assertRaises(FragmentLengthsError):
            match_candidates([DUAL6], ECORI, 6)  # type: ignore[arg-type]

    def test_fragment_lengths_string_rejected(self):
        with self.assertRaises(FragmentLengthsError):
            match_candidates([DUAL6], ECORI, "66")  # type: ignore[arg-type]

    def test_bad_candidate_base_located(self):
        with self.assertRaises(SequenceError) as cm:
            match_candidates([DUAL6, "GGNAGG"], ECORI, [6, 6])
        self.assertIn("candidates[1]", str(cm.exception))

    def test_candidate_wrong_type(self):
        with self.assertRaises(SequenceError) as cm:
            match_candidates([123], ECORI, [6, 6])  # type: ignore[list-item]
        self.assertIn("candidates[0]", str(cm.exception))

    def test_candidates_bare_string_rejected(self):
        with self.assertRaises(SequenceError):
            match_candidates(DUAL6, ECORI, [6, 6])  # type: ignore[arg-type]

    def test_empty_candidates_allowed(self):
        res = match_candidates([], ECORI, [6, 6])
        self.assertEqual(res.candidate_count, 0)
        self.assertEqual(res.class_count, 0)
        self.assertEqual(res.matches, ())
        self.assertEqual(res.unsupported, ())

    def test_enzyme_errors_propagate(self):
        with self.assertRaises(EnzymeDefinitionError):
            match_candidates([DUAL6], [], [6, 6])
        with self.assertRaises(EnzymeDefinitionError):
            match_candidates([DUAL6], ["not an enzyme"], [6, 6])  # type: ignore[list-item]


class TestIndependentReference(unittest.TestCase):
    """与测试文件内独立小参考的交叉核对（不调用被测核心生成期望值）。"""

    def test_reference_fragment_lengths(self):
        cases = (
            (DUAL6, (6, 6)),
            (TRIPLE6, (6, 6, 6)),
            (UNC12, (12,)),
            (ONE12, (12,)),
        )
        for seq, expected in cases:
            # 独立参考先确认手算期望值本身。
            self.assertEqual(ref_circular_lengths(seq, "GAATTC", 1), expected)
            res = match_candidates([seq], ECORI, list(expected))
            self.assertEqual(len(res.matches), 1)
            self.assertEqual(res.matches[0].fragment_lengths, expected)

    def test_reference_multi_enzyme(self):
        # BOTH：EcoRI 切口 1、SmaI 切口 9 -> 独立列窗核验间距 (4, 8)。
        ecori_cuts = [1]
        smai_cuts = sorted(
            (s + 3) % 12
            for s in range(12)
            if all(BOTH[(s + i) % 12] == "CCCGGG"[i] for i in range(6))
        )
        self.assertEqual(smai_cuts, [9])
        cuts = sorted(ecori_cuts + smai_cuts)
        gaps = tuple(sorted(
            (cuts[(i + 1) % 2] - cuts[i]) % 12 for i in range(2)))
        self.assertEqual(gaps, (4, 8))
        res = match_candidates([BOTH], [ECORI, SMAI], [4, 8])
        self.assertEqual(len(res.matches), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
