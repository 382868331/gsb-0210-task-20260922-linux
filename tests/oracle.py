"""独立暴力预言机：不导入、不调用 re_digest，用最直白的方式独立计算。

用途：与被测核心的结果交叉核对（切口坐标、片段弧段、突出端序列）。
期望值由本文件的独立实现现场算出，而非由被测核心生成后回灌。
"""

from __future__ import annotations

_IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT",
    "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    "N": "ACGT",
}
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def rc_dna(s: str) -> str:
    return "".join(_COMP[c] for c in reversed(s))


def rc_pattern(p: str) -> str:
    # 码集合互补：逐码独立计算，不查被测核心的互补表
    out = []
    for ch in reversed(p):
        comp_set = {_COMP[b] for b in _IUPAC[ch]}
        for code, bases in _IUPAC.items():
            if set(bases) == comp_set:
                out.append(code)
                break
    return "".join(out)


def _matches(seq: str, pattern: str, s: int, circ: bool) -> bool:
    n = len(seq)
    for k, ch in enumerate(pattern):
        base = seq[(s + k) % n] if circ else seq[s + k]
        if base not in _IUPAC[ch]:
            return False
    return True


def find_cuts(seq, rec, ot, ob, circ):
    """返回去重排序后的 (上切, 下切) 列表，以及全部命中 (start, dir)。"""
    n = len(seq)
    L = len(rec)
    rrec = rc_pattern(rec)
    events: dict[tuple[int, int], set[str]] = {}
    hits: list[tuple[int, str]] = []
    starts = range(n) if circ else range(n - L + 1)
    for s in starts:
        if _matches(seq, rec, s, circ):
            t, b = s + ot, s + ob
            ok = True
            if circ:
                t, b = t % n, b % n
            elif not (0 < t < n and 0 < b < n):
                ok = False
            if ok:
                events.setdefault((t, b), set()).add(f"F{s}")
                hits.append((s, "forward"))
        if _matches(seq, rrec, s, circ):
            t, b = s + L - ob, s + L - ot
            ok = True
            if circ:
                t, b = t % n, b % n
            elif not (0 < t < n and 0 < b < n):
                ok = False
            if ok:
                events.setdefault((t, b), set()).add(f"R{s}")
                hits.append((s, "reverse"))
    return sorted(events), sorted(hits)


def arc(seq: str, start: int, length: int, circ: bool) -> str:
    n = len(seq)
    if length == 0:
        return ""
    if circ:
        return "".join(seq[(start + k) % n] for k in range(length))
    return seq[start:start + length]


def fragment_arcs(seq, cuts, circ):
    """独立计算每个片段的 (上链弧, 下链5->3弧)。cuts 为排序后的 (t,b)。"""
    n = len(seq)
    out = []
    if not cuts:
        if circ:
            return [(seq, rc_dna(seq))]
        return [(seq, rc_dna(seq))]
    if circ:
        for i, (t, b) in enumerate(cuts):
            t2, b2 = cuts[(i + 1) % len(cuts)]
            tl = n if len(cuts) == 1 else (t2 - t) % n
            bl = n if len(cuts) == 1 else (b2 - b) % n
            out.append((arc(seq, t, tl, True), rc_dna(arc(seq, b, bl, True))))
    else:
        bounds = [(0, 0)] + cuts + [(n, n)]
        for (a1, b1), (a2, b2) in zip(bounds, bounds[1:]):
            out.append((arc(seq, a1, a2 - a1, False),
                        rc_dna(arc(seq, b1, b2 - b1, False))))
    return out


def expected_overhangs(seq, cuts, ot, ob, circ):
    """独立计算各物理切口从两侧看的突出端。

    返回 { (t,b): (left_kind, left_oh, right_kind, right_oh) }，blunt 为 ("blunt","")。
    左=沿顺时针方向离开该切口的片段左端；右=沿顺时针到达该切口的片段右端。
    """
    n = len(seq)
    d = ob - ot
    h = abs(d)
    res = {}
    for (t, b) in cuts:
        if d == 0:
            res[(t, b)] = ("blunt", "", "blunt", "")
        elif d > 0:
            # 左端：上链 5' 突出，读 t..t+h；右端：下链 5' 突出
            res[(t, b)] = (
                "5prime", arc(seq, t, h, circ),
                "5prime", rc_dna(arc(seq, t, h, circ)),
            )
        else:
            # 左端：下链 3' 突出（上链位置 t-h..t，5'->3' 书写）
            # 右端：上链 3' 突出（读 t-h..t）
            res[(t, b)] = (
                "3prime", rc_dna(arc(seq, t - h, h, circ)),
                "3prime", arc(seq, t - h, h, circ),
            )
    return res
