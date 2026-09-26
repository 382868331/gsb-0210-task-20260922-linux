"""环状部分酶切期望的独立验收测试。

小例参考通过显式枚举实际切点的全部 0/1 子集计算，不调用生产代码中的
O(m²) 期望公式。
"""

import itertools
import unittest

from restriction_digest import (
    DigestError,
    Enzyme,
    digest,
    partial_digest_expectation,
)


AA_AFTER = Enzyme("AAfter", "A", 1, 1)
ECORI = Enzyme("EcoRI", "GAATTC", 1, 5)


def subset_reference(n, cuts_with_probability):
    """cuts_with_probability: ((上链坐标, 概率), ...)，返回 ((len,kind),count)。"""
    expected = {}
    positions = [position for position, _probability in cuts_with_probability]
    probabilities = [probability for _position, probability in cuts_with_probability]

    for selected in itertools.product((False, True), repeat=len(positions)):
        subset_probability = 1.0
        for chosen, probability in zip(selected, probabilities, strict=True):
            subset_probability *= probability if chosen else 1.0 - probability
        chosen_positions = sorted(
            position for position, chosen in zip(positions, selected, strict=True)
            if chosen
        )
        if not chosen_positions:
            pieces = [((n, "circular"), 1)]
        elif len(chosen_positions) == 1:
            pieces = [((n, "linear"), 1)]
        else:
            pieces = []
            for i, start in enumerate(chosen_positions):
                end = chosen_positions[(i + 1) % len(chosen_positions)]
                pieces.append((((end - start) % n, "linear"), 1))
        for key, count in pieces:
            expected[key] = expected.get(key, 0.0) + subset_probability * count
    return expected


class TestPartialDigestExpectation(unittest.TestCase):
    def test_three_adjacent_sites_matches_full_subsets(self):
        seq = "AAACCCCCCCCC"  # 12 bp；A 后切点在 1、2、3，彼此相邻
        probabilities = (0.2, 0.5, 0.8)
        result = partial_digest_expectation(
            seq, AA_AFTER, probabilities
        )

        self.assertEqual(len(result.merged_cut_sites), 3)
        self.assertEqual(
            [site.top_boundary for site in result.merged_cut_sites],
            [1, 2, 3],
        )
        expected = subset_reference(
            12,
            ((1, 0.2), (2, 0.5), (3, 0.8)),
        )
        actual = {
            (fragment.length, fragment.kind): fragment.expected_count
            for fragment in result.expected_fragments
        }
        self.assertEqual(set(actual), set(expected))
        for key, value in expected.items():
            self.assertAlmostEqual(actual[key], value)
        self.assertAlmostEqual(result.expected_total_length, 12)

    def test_all_probability_zero_keeps_only_uncut_circle(self):
        seq = "AAACCCCCCCCC"
        result = partial_digest_expectation(seq, AA_AFTER, 0.0)
        self.assertEqual(
            [(f.length, f.kind, f.expected_count) for f in result.expected_fragments],
            [(12, "circular", 1.0)],
        )

    def test_all_probability_one_gives_only_complete_digest_arcs(self):
        seq = "AAACCCCCCCCC"
        result = partial_digest_expectation(seq, AA_AFTER, 1.0)
        counts = {
            (fragment.length, fragment.kind): fragment.expected_count
            for fragment in result.expected_fragments
        }
        self.assertEqual(counts, {(1, "linear"): 2.0, (10, "linear"): 1.0})
        self.assertAlmostEqual(result.expected_total_length, 12)

    def test_two_adjacent_sites(self):
        seq = "AACCCCCCCC"  # 10 bp；切点在 1、2
        result = partial_digest_expectation(seq, AA_AFTER, (0.25, 0.75))
        expected = subset_reference(10, ((1, 0.25), (2, 0.75)))
        actual = {
            (fragment.length, fragment.kind): fragment.expected_count
            for fragment in result.expected_fragments
        }
        self.assertEqual(set(actual), set(expected))
        for key, value in expected.items():
            self.assertAlmostEqual(actual[key], value)
        self.assertAlmostEqual(result.expected_total_length, 10)

    def test_coincident_sites_merge_as_at_least_one_cut(self):
        seq = "GGGAATTCCCAA"
        clone = Enzyme("EcoRI-clone", "GAATTC", 1, 5)

        complete = digest(seq, [ECORI, clone], topology="circular")
        self.assertEqual(len(complete.cut_sites), 2)
        result = partial_digest_expectation(
            seq,
            [ECORI, clone],
            (0.5, 0.5),
        )

        self.assertEqual(len(result.merged_cut_sites), 1)
        merged = result.merged_cut_sites[0]
        self.assertEqual((merged.top_boundary, merged.bottom_boundary), (3, 7))
        self.assertAlmostEqual(merged.probability, 0.75)
        self.assertEqual(len(merged.evidence), 2)
        counts = {
            (fragment.length, fragment.kind): fragment.expected_count
            for fragment in result.expected_fragments
        }
        self.assertEqual(set(counts), {(12, "circular"), (12, "linear")})
        self.assertAlmostEqual(counts[(12, "circular")], 0.25)
        self.assertAlmostEqual(counts[(12, "linear")], 0.75)
        self.assertAlmostEqual(result.expected_total_length, 12)

    def test_callable_probabilities_align_with_identified_sites(self):
        seq = "AAACCCCCCCCC"

        def probability(site):
            return {0: 0.1, 1: 0.2, 2: 0.3}[site.start]

        result = partial_digest_expectation(seq, AA_AFTER, probability)
        probabilities = [site.probability for site in result.merged_cut_sites]
        self.assertEqual(len(probabilities), 3)
        for actual, expected in zip(probabilities, (0.1, 0.2, 0.3), strict=True):
            self.assertAlmostEqual(actual, expected)

    def test_invalid_probability_and_sequence_length(self):
        seq = "AACCCCCCCC"
        with self.assertRaises(DigestError):
            partial_digest_expectation(seq, AA_AFTER, (0.5,))
        with self.assertRaises(DigestError):
            partial_digest_expectation(seq, AA_AFTER, 1.5)
        with self.assertRaises(DigestError):
            partial_digest_expectation(seq, AA_AFTER, (0.5, float("nan")))

    def test_only_circular_topology_supported(self):
        with self.assertRaises(DigestError):
            partial_digest_expectation("AACCCCCCCC", AA_AFTER, 0.5, "linear")


if __name__ == "__main__":
    unittest.main(verbosity=2)
