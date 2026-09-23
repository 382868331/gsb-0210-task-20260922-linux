"""交叉验证：确定性随机用例对照独立预言机 tests.oracle（不调用被测核心造期望）。

另含物理不变量检查（仅作为补充）：上链弧串联复原原序列、下链亦然、
同一物理切口两侧末端必然可连接。
"""

from __future__ import annotations

import random
import unittest

from re_digest import Enzyme, UnsupportedDigestError, digest

try:  # 兼容 `python -m unittest discover -s tests`（tests 在 sys.path 上）
    from . import oracle
except ImportError:  # pragma: no cover
    import oracle


def _gap_supported(cuts, n, h, circ):
    """独立判定切口间距是否足够（与核心拒绝条件互为独立写法）。"""
    ts = sorted(t for t, _ in cuts)
    if not ts:
        return True
    if not circ:
        return all(b - a > h for a, b in zip(ts, ts[1:]))
    if len(ts) == 1:
        return n > h
    return all((b - a) % n > h for a, b in zip(ts, ts[1:] + ts[:1]))


class OracleCrossCheckTests(unittest.TestCase):
    def test_random_cases(self):
        rng = random.Random(20260923)
        enzymes = [
            ("G", 0, 0), ("A", 0, 2), ("C", 2, 0), ("AT", 0, 0),
            ("GC", 1, 3), ("GC", 3, 1), ("R", 0, 0), ("N", 1, 1),
            ("GATC", 1, 5), ("AAG", -1, 2), ("GT", 0, 1),
        ]
        checked = 0
        for trial in range(400):
            circ = trial % 2 == 0
            n = rng.randint(1, 40)
            seq = "".join(rng.choice("ACGT") for _ in range(n))
            rec, ot, ob = rng.choice(enzymes)
            if not circ and n < len(rec):
                continue
            cuts, hits = oracle.find_cuts(seq, rec, ot, ob, circ)
            h = abs(ob - ot)
            enz = Enzyme("rand", rec, ot, ob)
            if not _gap_supported(cuts, n, h, circ):
                with self.assertRaises(UnsupportedDigestError):
                    digest(seq, enz, circular=circ)
                continue
            res = digest(seq, enz, circular=circ)
            # 切口与命中
            self.assertEqual(
                [(e.top_cut, e.bottom_cut) for e in res.cut_events], cuts
            )
            self.assertEqual(
                [(m.start, m.orientation) for m in res.site_matches], hits
            )
            # 片段双链弧段
            got_arcs = [(f.top_sequence, f.bottom_sequence_5to3)
                        for f in res.fragments]
            self.assertEqual(got_arcs, oracle.fragment_arcs(seq, cuts, circ))
            # 突出端
            want_oh = oracle.expected_overhangs(seq, cuts, ot, ob, circ)
            for ev, frag_pair in zip(res.cut_events,
                                    self._event_fragments(res, circ)):
                lk, lo, rk, ro = want_oh[(ev.top_cut, ev.bottom_cut)]
                lf, rf = frag_pair
                self.assertEqual((lf.kind, lf.overhang), (lk, lo))
                self.assertEqual((rf.kind, rf.overhang), (rk, ro))
                # 同切口两侧必互补可连（平端自然互连）
                self.assertTrue(lf.compatible_with(rf))
            checked += 1
        self.assertGreater(checked, 100)  # 确认确实跑了大量有效用例

    @staticmethod
    def _event_fragments(res, circ):
        """对每个物理切口给出 (其左侧片段的左末端?, 实际配对)。

        返回每个切口的 (离开片段左端, 到达片段右端) 末端对。
        """
        frags = res.fragments
        if not res.cut_events:
            return []
        if circ:
            # 片段 i 左端 = 事件 i；片段 (i-1) 右端 = 事件 i
            return [
                (frags[i].left_end, frags[(i - 1) % len(frags)].right_end)
                for i in range(len(frags))
            ]
        pairs = []
        for i, ev in enumerate(res.cut_events):
            # 线性：事件 i 左侧是片段 i 的右端，右侧是片段 i+1 的左端
            pairs.append((frags[i + 1].left_end, frags[i].right_end))
        return pairs

    def test_physical_reassembly_invariants(self):
        rng = random.Random(424242)
        for trial in range(120):
            circ = trial % 3 == 0
            n = rng.randint(4, 30)
            seq = "".join(rng.choice("ACGT") for _ in range(n))
            rec, ot, ob = rng.choice([("GC", 1, 3), ("GC", 3, 1),
                                      ("AT", 0, 0), ("G", 0, 0)])
            enz = Enzyme("inv", rec, ot, ob)
            cuts, _ = oracle.find_cuts(seq, rec, ot, ob, circ)
            if not _gap_supported(cuts, n, abs(ob - ot), circ):
                continue
            res = digest(seq, enz, circular=circ)
            if circ and not cuts:
                self.assertEqual(res.fragments[0].shape, "circular")
                continue
            # 上链弧按序串联复原原序列（循环移位意义下）
            top_concat = "".join(f.top_sequence for f in res.fragments)
            bot_concat = "".join(oracle.rc_dna(f.bottom_sequence_5to3)
                                 for f in res.fragments)
            if circ:
                self.assertEqual(len(top_concat), n)
                self.assertIn(seq, top_concat + top_concat)
                self.assertIn(seq, bot_concat + bot_concat)
            else:
                self.assertEqual(top_concat, seq)
                self.assertEqual(bot_concat, seq)
            # 每个片段：两端都黏/都平 -> 两链等长；恰一黏 -> 链长差 = h
            h = abs(ob - ot)
            for f in res.fragments:
                sticky_ends = int(f.left_end.kind != "blunt") + int(
                    f.right_end.kind != "blunt"
                )
                if sticky_ends in (0, 2):
                    self.assertEqual(
                        len(f.top_sequence), len(f.bottom_sequence_5to3)
                    )
                else:
                    self.assertEqual(
                        abs(len(f.top_sequence)
                            - len(f.bottom_sequence_5to3)),
                        h,
                    )


if __name__ == "__main__":
    unittest.main()
