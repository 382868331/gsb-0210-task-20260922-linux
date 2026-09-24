"""候选环状切位图 × 片段长度多重集匹配的独立验收测试。

期望值全部手工计算，或由本文件内独立编写的小规模全枚举（indep_*）
交叉核对；独立参考不调用 restriction_digest 的匹配/归一化实现来
产生期望值，只用来核对被测函数的结论。
"""

import math
import unittest
from collections import Counter

from restriction_digest import (
    CandidateError,
    CutMap,
    DigestResult,
    Enzyme,
    MatchResult,
    circular_fragment_lengths,
    digest,
    map_from_digest,
    match_candidates,
)


SMAI = Enzyme("SmaI", "CCCGGG", 3, 3)  # CCC^GGG，平端
ECORI = Enzyme("EcoRI", "GAATTC", 1, 5)


# ---------------------------------------------------------------------------
# 测试文件自己的独立参考（显式公式 / 全枚举，与被测核心分开编写）
# ---------------------------------------------------------------------------


def indep_gaps(n, cuts):
    """独立公式：切口升序时的循环间隙元组（未切环为 (n,)）。"""
    s = sorted(cuts)
    if not s:
        return (n,)
    out = [s[i + 1] - s[i] for i in range(len(s) - 1)]
    out.append(n - s[-1] + s[0])
    return tuple(out)


def indep_orbit(cuts, n):
    """一张切位图在所有环旋转下的位置集合（frozenset）轨道。"""
    base = frozenset(cuts)
    return {frozenset((c + d) % n for c in base) for d in range(n)}


class TestHandComputedLengths(unittest.TestCase):
    """手算：环长 12 的几张图的完全消化片段长度。"""

    def test_basic_gaps(self):
        # 切口 0,2,4,8：间距 2,2,4，绕回 12-8+0=4
        self.assertEqual(circular_fragment_lengths(12, (0, 2, 4, 8)), (2, 2, 4, 4))
        # 切口 0,2,6,8：间距 2,4,2，绕回 4
        self.assertEqual(circular_fragment_lengths(12, (0, 2, 6, 8)), (2, 4, 2, 4))

    def test_uncut_and_single_cut(self):
        # 未切环：整环一段
        self.assertEqual(circular_fragment_lengths(12, ()), (12,))
        # 单切口（任意位置）：线性化，仍是一整圈长度
        self.assertEqual(circular_fragment_lengths(12, (5,)), (12,))
        self.assertEqual(circular_fragment_lengths(12, (0,)), (12,))

    def test_cutmap_method_matches_function(self):
        m = CutMap(12, (0, 2, 4, 8))
        self.assertEqual(m.fragment_lengths(), (2, 2, 4, 4))
        self.assertEqual(m.cut_count, 4)
        self.assertFalse(m.is_uncut)
        self.assertTrue(CutMap(12, ()).is_uncut)


class TestNormalMatching(unittest.TestCase):
    """正常匹配：长度精确、顺序无关、重复片段。"""

    def test_two_nonequivalent_arrangements_same_multiset(self):
        # 两张图的循环排列不同（(2,2,4,4) vs (2,4,2,4)），
        # 多重集相同但互不旋转等价 -> 都必须保留为匹配。
        c1 = CutMap(12, (0, 2, 4, 8), label="A")
        c2 = CutMap(12, (0, 2, 6, 8), label="B")
        res = match_candidates([2, 2, 4, 4], [c1, c2])
        self.assertIsInstance(res, MatchResult)
        self.assertTrue(res.matched)
        self.assertEqual([m.label for m in res.matches], ["A", "B"])
        self.assertEqual(res.examined, 2)
        self.assertEqual(res.unique_classes, 2)

    def test_observed_order_and_type_independent(self):
        c1 = CutMap(12, (0, 2, 4, 8))
        r1 = match_candidates((4, 2, 4, 2), [c1])
        r2 = match_candidates(iter([2, 4, 2, 4]), [c1])
        self.assertEqual(len(r1.matches), 1)
        self.assertEqual(len(r2.matches), 1)

    def test_repeated_fragments(self):
        # 四等分环：四段全为 3
        m = CutMap(12, (0, 3, 6, 9), label="quarters")
        res = match_candidates([3, 3, 3, 3], [m, CutMap(12, (0, 4, 8))])
        self.assertEqual([m.label for m in res.matches], ["quarters"])

    def test_conservation_sum_is_necessary_not_sufficient(self):
        # 长度和同为 12 但多重集不同：不匹配（精确比较，不容差）。
        res = match_candidates([2, 2, 4, 4], [CutMap(12, (0, 3, 6, 9))])
        self.assertFalse(res.matched)
        self.assertEqual(res.matches, ())
        # 和不同（13）也不匹配
        res2 = match_candidates([2, 2, 4, 5], [CutMap(12, (0, 2, 4, 8))])
        self.assertFalse(res2.matched)

    def test_no_match_is_proven_empty(self):
        res = match_candidates([5, 7], [CutMap(12, (0, 2, 4, 8))])
        self.assertFalse(res.matched)
        self.assertEqual(res.matches, ())
        self.assertEqual(res.examined, 1)
        self.assertEqual(res.unique_classes, 1)


class TestRotationalDedup(unittest.TestCase):
    """旋转等价去重：保留先出现者，不影响计数字段。"""

    def test_rotated_duplicate_collapsed_to_first(self):
        # (1,3,5,9) 是 (0,2,4,8) 整体旋转 +1：间隙同为 (2,2,4,4)。
        first = CutMap(12, (0, 2, 4, 8), label="first")
        rotated = CutMap(12, (1, 3, 5, 9), label="rotated")
        res = match_candidates([2, 2, 4, 4], [first, rotated])
        self.assertEqual([m.label for m in res.matches], ["first"])
        self.assertEqual(res.examined, 2)
        self.assertEqual(res.unique_classes, 1)

    def test_reverse_order_keeps_other_representative(self):
        first = CutMap(12, (1, 3, 5, 9), label="first")
        rotated = CutMap(12, (0, 2, 4, 8), label="rotated")
        res = match_candidates([2, 2, 4, 4], [first, rotated])
        self.assertEqual([m.label for m in res.matches], ["first"])

    def test_nonmatching_duplicate_still_deduped_in_class_count(self):
        # 一对旋转等价但与观测不符的候选：不匹配，但只算一个等价类。
        a = CutMap(12, (0, 1, 4, 8))
        b = CutMap(12, (1, 2, 5, 9))  # a 旋转 +1
        res = match_candidates([3, 3, 3, 3], [a, b])
        self.assertFalse(res.matched)
        self.assertEqual(res.examined, 2)
        self.assertEqual(res.unique_classes, 1)

    def test_uncut_and_single_cut_are_distinct_classes(self):
        # 未切环与单切环长度多重集都是 (12,)：都匹配，但互不等价，
        # 各自作为独立等价类保留。
        uncut = CutMap(12, (), label="uncut")
        one = CutMap(12, (5,), label="one-cut")
        one_rot = CutMap(12, (7,), label="one-cut-rotated")  # 与 one 等价
        res = match_candidates([12], [uncut, one, one_rot])
        self.assertEqual([m.label for m in res.matches], ["uncut", "one-cut"])
        self.assertEqual(res.examined, 3)
        self.assertEqual(res.unique_classes, 2)


class TestMapFromDigest(unittest.TestCase):
    """复用旧切割规则：从环状 digest 结果生成候选切位图。"""

    def test_single_blunt_site(self):
        # "AAACCCGGGTTT"（12 bp 环）：SmaI 起点 3，平端切口 3+3=6。
        seq = "AAACCCGGGTTT"
        res = digest(seq, SMAI, topology="circular")
        cm = map_from_digest(res)
        self.assertEqual((cm.length, cm.cuts), (12, (6,)))
        self.assertEqual(cm.fragment_lengths(), (12,))

    def test_two_sites_gaps_hand_checked(self):
        # "CCCGGGAAACCCGGG"（15 bp 环）：位点起点 0 与 9 -> 切口 3 与 12，
        # 间隙 12-3=9 与 15-12+3=6。
        seq = "CCCGGGAAACCCGGG"
        self.assertEqual(len(seq), 15)
        cm = map_from_digest(digest(seq, SMAI, topology="circular"))
        self.assertEqual(cm.cuts, (3, 12))
        self.assertEqual(cm.fragment_lengths(), (9, 6))
        # 用测得的片段长度（顺序打乱）在候选池中匹配，只有这张图命中
        others = [
            CutMap(15, (0, 5, 10)),     # (5,5,5)
            CutMap(15, (0, 5)),         # (5,10)
            CutMap(15, ()),             # 未切 (15,)
        ]
        hit = match_candidates([6, 9], [*others, cm])
        self.assertEqual(hit.matches, (cm,))
        # 旋转等价的另一张表示图 (3,9) -> 间隙 (6,9)：与 cm 同类，
        # 且因 cm 先出现而被去重，不产生第二个匹配。
        equiv = CutMap(15, (3, 9))
        hit2 = match_candidates([6, 9], [cm, equiv])
        self.assertEqual(hit2.matches, (cm,))
        self.assertEqual(hit2.unique_classes, 1)

    def test_uncut_circle(self):
        cm = map_from_digest(digest("ACGTACGT", ECORI, topology="circular"))
        self.assertEqual(cm.cuts, ())
        self.assertTrue(cm.is_uncut)

    def test_rejects_linear_and_foreign(self):
        lin = digest("AAACCCGGGTT", SMAI)
        with self.assertRaises(CandidateError):
            map_from_digest(lin)
        with self.assertRaises(CandidateError):
            map_from_digest("not a result")  # type: ignore[arg-type]


class TestExhaustiveCrossCheck(unittest.TestCase):
    """小规模全枚举交叉核对：n=6 的全部 64 张切位图。

    独立枚举旋转轨道与间隙多重集，要求被测去重/匹配结论与枚举逐一相同。
    """

    N = 6

    def _all_subset_maps(self):
        maps = []
        for mask in range(1 << self.N):
            cuts = tuple(p for p in range(self.N) if mask & (1 << p))
            maps.append(CutMap(self.N, cuts, label=f"m{mask:02d}"))
        return maps

    def test_class_count_and_matches_against_enumeration(self):
        n = self.N
        maps = self._all_subset_maps()

        # 独立计算旋转轨道代表（固定顺序首次出现）与轨道总数。
        seen_orbits: set[frozenset] = set()
        indep_reps: list[str] = []
        for cm in maps:
            orbit = indep_orbit(cm.cuts, n)
            if orbit & seen_orbits:
                continue
            seen_orbits |= orbit
            indep_reps.append(cm.label)

        # 目标：取 mask=13={0,2,3} 的间隙多重集（独立计算）。
        target = indep_gaps(n, (0, 2, 3))
        target_counter = Counter(target)
        expected_match_labels = [
            label
            for label in indep_reps
            if Counter(indep_gaps(n, next(
                m.cuts for m in maps if m.label == label))) == target_counter
        ]

        res = match_candidates(target, maps)
        self.assertEqual(res.examined, len(maps))
        self.assertEqual(res.unique_classes, len(indep_reps))
        self.assertEqual([m.label for m in res.matches], expected_match_labels)
        # 手算 {0,2,3} 的间隙：2,1,绕回 3 -> (2,1,3)
        self.assertEqual(sorted(target), [1, 2, 3])
        # 该长度多重集在 n=6 上至少有两个不旋转等价的排列
        self.assertGreaterEqual(len(expected_match_labels), 2)

    def test_every_map_agrees_with_independent_gaps(self):
        for cm in self._all_subset_maps():
            self.assertEqual(
                cm.fragment_lengths(), indep_gaps(self.N, cm.cuts)
            )


class TestCandidateValidation(unittest.TestCase):
    """非法输入：可定位错误，原子失败（纯函数无半更新状态）。"""

    def test_bad_length(self):
        for bad in (0, -1):
            with self.assertRaises(CandidateError):
                CutMap(bad, ())
        for bad in (True, 1.5, math.nan, math.inf):
            with self.assertRaises(CandidateError):
                CutMap(bad, ())  # type: ignore[arg-type]
        with self.assertRaises(CandidateError):
            CutMap("12", ())  # type: ignore[arg-type]

    def test_bad_cuts(self):
        with self.assertRaises(CandidateError) as cm:
            CutMap(12, (12,))            # 上界越界（边界范围 0..n-1）
        self.assertIn("cuts[0]", str(cm.exception))
        with self.assertRaises(CandidateError):
            CutMap(12, (-1,))
        with self.assertRaises(CandidateError):
            CutMap(12, (1, 1))           # 重复
        with self.assertRaises(CandidateError) as cm2:
            CutMap(12, (3, 2))           # 非升序
        self.assertIn("cuts[1]", str(cm2.exception))
        with self.assertRaises(CandidateError):
            CutMap(12, (1, 2.0))         # type: ignore[arg-type]
        with self.assertRaises(CandidateError):
            CutMap(12, (1, False))       # type: ignore[arg-type]
        with self.assertRaises(CandidateError):
            CutMap(12, "012")            # type: ignore[arg-type]
        with self.assertRaises(CandidateError):
            CutMap(12, 3)                # type: ignore[arg-type]

    def test_bad_label(self):
        with self.assertRaises(CandidateError):
            CutMap(12, (), label=5)  # type: ignore[arg-type]

    def test_bad_observed(self):
        good = [CutMap(12, (0, 2, 4, 8))]
        with self.assertRaises(CandidateError):
            match_candidates([], good)                  # 空多重集
        with self.assertRaises(CandidateError) as cm:
            match_candidates([2, 0, 4, 4], good)        # 0 长片段
        self.assertIn("observed[1]", str(cm.exception))
        with self.assertRaises(CandidateError):
            match_candidates([2, -2], good)
        with self.assertRaises(CandidateError):
            match_candidates([2, 4.0], good)            # type: ignore[list-item]
        with self.assertRaises(CandidateError):
            match_candidates([2, True], good)           # type: ignore[list-item]
        for bad in (math.nan, math.inf):
            with self.assertRaises(CandidateError):
                match_candidates([bad], good)           # type: ignore[list-item]
        with self.assertRaises(CandidateError):
            match_candidates(12, good)                  # type: ignore[arg-type]

    def test_bad_candidate_member(self):
        with self.assertRaises(CandidateError):
            match_candidates([12], ["not a map"])  # type: ignore[list-item]

    def test_generator_with_invalid_map_raises_atomically(self):
        # 生成器先产出合法图，再产出构造即非法的图：整体失败，
        # 调用方拿不到任何部分结果（纯函数，异常即无返回）。
        def gen():
            yield CutMap(12, (0, 2, 4, 8))
            yield CutMap(12, (5, 5))  # 重复切口

        with self.assertRaises(CandidateError):
            match_candidates([2, 2, 4, 4], gen())

    def test_inputs_not_mutated(self):
        cuts = [0, 2, 4, 8]
        observed = [4, 2, 4, 2]
        cands = [CutMap(12, tuple(cuts))]
        match_candidates(observed, cands)
        self.assertEqual(cuts, [0, 2, 4, 8])
        self.assertEqual(observed, [4, 2, 4, 2])

    def test_result_is_readonly(self):
        res = match_candidates([12], [CutMap(12, ())])
        with self.assertRaises(Exception):
            res.examined = 99  # type: ignore[misc]
        with self.assertRaises(Exception):
            res.matches[0].cuts = (1,)  # type: ignore[misc]

    def test_result_fields_frozen_dataclass(self):
        # DigestResult/CutMap/MatchResult 均为 frozen dataclass
        with self.assertRaises(Exception):
            CutMap(12, (1,)).length = 9  # type: ignore[misc]
        self.assertIsInstance(map_from_digest(
            digest("AAACCCGGGTTT", SMAI, topology="circular")), CutMap)


if __name__ == "__main__":
    unittest.main(verbosity=2)
