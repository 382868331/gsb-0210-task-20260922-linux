"""限制性内切酶环状片段归属核心（标准库实现）。

坐标约定（本模块公开契约的一部分）
----------------------------------
1. 碱基沿上链 5'->3' 编号 0..n-1；切口 ``c`` 表示“碱基 c 之前”的磷酸二酯键。
2. 酶的上链/下链切偏移均相对“匹配起点”（匹配串第一个碱基所在位置 s）：
       上链切口 = s + top_offset
       下链切口 = s + bottom_offset
   例：EcoRI 样酶 G^AATTC / CTTAA^G，top_offset=1, bottom_offset=5。
3. top_offset < bottom_offset 产生 5' 突出端；top_offset > bottom_offset 产生
   3' 突出端；相等为平端。突出端长度 d = bottom_offset - top_offset。
4. 线性 DNA：只有两个链的切口坐标都严格位于 (0, n) 时该识别事件才真正切断
   （任一切点落在分子末端之外 = 切到空气，整事件不切，这是物理事实而非纠正
   输入）。环状 DNA：坐标对 n 取模，0 与 n 为同一物理键。
5. 反向互补命中：若反向互补模式在上链 j..j+L-1 命中，则
       上链切口 = j + L - bottom_offset
       下链切口 = j + L - top_offset
6. 相同 (上链切口, 下链切口) 的物理双切事件去重（回文酶的正/反向命中即如此）。
7. 当相邻上链切口的环上间距 <= |d| 时，切口互相套叠，无法归属为彼此独立的
   常规双链片段，抛出 :class:`UnsupportedDigestError`，绝不任意配对。

搜索为穷举（O(n*L)，n<=5000, L<=12），因此“无命中”是证明无解，不存在
搜索预算耗尽这一模式。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal

# ----------------------------------------------------------------------------
# 公开常量
# ----------------------------------------------------------------------------

MAX_SEQUENCE_LEN = 5000
"""允许的 DNA 长度（碱基），含 5000。"""

MAX_RECOGNITION_LEN = 12
"""识别串最大长度。"""

MAX_OFFSET_ABS = 5000
"""切偏移绝对值的接口护栏（>= 最大分子长度时任何情况下都不可能切中）。"""

# IUPAC 核酸简并码 -> 可匹配的具体碱基集合
IUPAC_BASES: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),       # 嘌呤
    "Y": frozenset("CT"),       # 嘧啶
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

# IUPAC 码互补（A<->T, C<->G, R<->Y, K<->M, B<->V, D<->H, S/W/N 自互补）
IUPAC_COMPLEMENT: dict[str, str] = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}

_DNA_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

EndKind = Literal["blunt", "5prime", "3prime"]
Orientation = Literal["forward", "reverse"]
FragmentShape = Literal["linear", "circular"]


# ----------------------------------------------------------------------------
# 错误类型
# ----------------------------------------------------------------------------

class DigestError(ValueError):
    """本模块所有消化相关输入/拓扑错误的基类（同时是 ValueError）。"""


class InvalidSequenceError(DigestError):
    """DNA 序列非法。

    :param index: 首个非法字符的 0 基下标；若错误与具体位置无关则为 None。
    """

    def __init__(self, message: str, index: int | None = None):
        super().__init__(message)
        self.index = index


class InvalidEnzymeError(DigestError):
    """酶定义非法。``field`` 给出可定位的出错字段名。"""

    def __init__(self, message: str, field: str | None = None):
        super().__init__(message)
        self.field = field


class UnsupportedDigestError(DigestError):
    """切口拓扑在本题契约子集内无法归属为简单双链片段（明确拒绝，不做猜测配对）。"""


# ----------------------------------------------------------------------------
# 输入校验
# ----------------------------------------------------------------------------

def _require_int(name: str, value: object) -> int:
    """整数参数：拒绝 bool，拒绝一切 float（含 NaN/Infinity）与其他类型。"""
    if isinstance(value, bool):
        raise TypeError(f"{name} 必须是 int，不能是 bool（得到 {value!r}）")
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if math.isnan(value):
            raise TypeError(f"{name} 必须是 int，收到 NaN")
        if math.isinf(value):
            raise TypeError(f"{name} 必须是 int，收到 Infinity")
        raise TypeError(f"{name} 必须是 int，收到 float: {value!r}")
    raise TypeError(f"{name} 必须是 int，收到 {type(value).__name__}")


def _require_str(name: str, value: object) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{name} 必须是 str，收到 {type(value).__name__}")
    return value


# ----------------------------------------------------------------------------
# 公开数据结构
# ----------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class Enzyme:
    """限制性内切酶定义。

    :param name: 人类可读名称（仅标签，允许空串）。
    :param recognition: IUPAC 识别串，长度 1..12。
    :param top_offset: 上链切偏移（相对匹配起点，可为负/越位，类型 IIS）。
    :param bottom_offset: 下链切偏移（同一参照系，见模块坐标约定）。
    """

    name: str
    recognition: str
    top_offset: int
    bottom_offset: int

    def __post_init__(self) -> None:
        if not isinstance(self.name, str):
            raise TypeError("name 必须是 str")
        rec = _require_str("recognition", self.recognition)
        top = _require_int("top_offset", self.top_offset)
        bot = _require_int("bottom_offset", self.bottom_offset)
        if not 1 <= len(rec) <= MAX_RECOGNITION_LEN:
            raise InvalidEnzymeError(
                f"识别串长度必须在 1..{MAX_RECOGNITION_LEN} 之间，得到 {len(rec)}",
                "recognition",
            )
        for i, ch in enumerate(rec):
            if ch not in IUPAC_BASES:
                raise InvalidEnzymeError(
                    f"识别串位置 {i} 的字符 {ch!r} 不是合法 IUPAC 码",
                    "recognition",
                )
        if abs(top) > MAX_OFFSET_ABS or abs(bot) > MAX_OFFSET_ABS:
            raise InvalidEnzymeError(
                f"切偏移绝对值不得超过 {MAX_OFFSET_ABS}（接口护栏）",
                "top_offset" if abs(top) > MAX_OFFSET_ABS else "bottom_offset",
            )

    @property
    def stagger(self) -> int:
        """下链切口 - 上链切口；>0 为 5' 突出，<0 为 3' 突出，0 为平端。"""
        return self.bottom_offset - self.top_offset


@dataclass(frozen=True, slots=True)
class SiteMatch:
    """一次识别命中（只记录真正产生切口的命中；线性末端越位命中不计）。"""

    start: int
    """匹配起点（上链坐标；环状 0..n-1，线性 0..n-L）。"""
    orientation: Orientation
    """forward=模式沿上链 5'->3'；reverse=反向互补模式沿上链命中。"""


@dataclass(frozen=True, slots=True)
class CutEvent:
    """去重后的一个物理双切事件。"""

    top_cut: int
    """上链切口坐标（线性 1..n-1；环状 0..n-1）。"""
    bottom_cut: int
    """下链切口坐标（同一坐标参照系）。"""
    match_starts: tuple[tuple[Orientation, int], ...] = field(default_factory=tuple)
    """产生该物理事件的全部命中（方向, 起点），按 (起点, 方向) 排序。"""


@dataclass(frozen=True, slots=True)
class EndInfo:
    """片段的一个末端。

    :param kind: blunt / 5prime / 3prime。
    :param overhang: 突出链序列，按该链 5'->3' 书写；平端为 ''。
    :param top_boundary: 该端的上链边界（碱基前键坐标）。
    :param bottom_boundary: 该端的下链边界。
    :param ligation_key: 可连接归类键——两端 kind 相同且键相等即可连接
           （平端互连；黏端要求突出序列互为反向互补，键内取规范形式）。
    """

    kind: EndKind
    overhang: str
    top_boundary: int
    bottom_boundary: int
    ligation_key: tuple[object, ...]

    def compatible_with(self, other: "EndInfo") -> bool:
        """两个末端是否可连接（端型一致且平端/互补黏端）。"""
        return self.kind == other.kind and self.ligation_key == other.ligation_key


@dataclass(frozen=True, slots=True)
class Fragment:
    """按上链切口分段的一个双链片段。

    ``shape='circular'`` 仅用于环状 DNA 完全未切时的整环表示，此时两个末端
    为 None；所有真实消化产物均为 ``shape='linear'`` 且带有两个末端。
    """

    index: int
    shape: FragmentShape
    top_start: int
    top_end: int
    bottom_start: int
    bottom_end: int
    top_sequence: str
    """上链序列，5'->3'（环状时可能跨零拼接）。"""
    bottom_sequence_5to3: str
    """下链序列，5'->3' 方向书写（等于其覆盖位置串的反向互补）。"""
    left_end: EndInfo | None
    right_end: EndInfo | None


@dataclass(frozen=True, slots=True)
class DigestResult:
    """消化结果。"""

    sequence: str
    circular: bool
    length: int
    enzyme_name: str
    site_matches: tuple[SiteMatch, ...]
    """全部产生切口的识别命中（正反向分别计数，回文位点会有两条）。"""
    cut_events: tuple[CutEvent, ...]
    """去重后的物理双切事件，按上链切口排序。"""
    fragments: tuple[Fragment, ...]
    """按上链切口顺序分段的双链片段。"""
    positions_scanned: int
    """穷举扫描的起点数（每个起点检查正、反两个方向）。"""


# ----------------------------------------------------------------------------
# 核心算法
# ----------------------------------------------------------------------------

def reverse_complement_pattern(pattern: str) -> str:
    """IUPAC 识别串的反向互补（码集合互补后反转）。"""
    return "".join(IUPAC_COMPLEMENT[ch] for ch in reversed(pattern))


def reverse_complement_dna(sequence: str) -> str:
    """具体 ACGT 序列的反向互补。"""
    return "".join(_DNA_COMPLEMENT[ch] for ch in reversed(sequence))


def _match_at(
    code_sets: tuple[frozenset[str], ...],
    sequence: str,
    start: int,
    circular: bool,
    n: int,
) -> bool:
    if circular:
        for k, allowed in enumerate(code_sets):
            if sequence[(start + k) % n] not in allowed:
                return False
        return True
    for k, allowed in enumerate(code_sets):
        if sequence[start + k] not in allowed:
            return False
    return True


def _arc(sequence: str, start: int, length: int, circular: bool, n: int) -> str:
    """从上链坐标 start 起顺时针读取 length 个碱基（允许跨零；length<=n）。"""
    if length == 0:
        return ""
    if circular:
        start %= n
        end = start + length
        if end <= n:
            return sequence[start:end]
        return sequence[start:] + sequence[: end - n]
    return sequence[start : start + length]


def _make_end(
    top_boundary: int,
    bottom_boundary: int,
    side: Literal["left", "right"],
    stagger: int,
    sequence: str,
    circular: bool,
    n: int,
) -> EndInfo:
    """根据同一双切事件在片段左/右端的两条链边界构造末端信息。"""
    if top_boundary == bottom_boundary:
        return EndInfo("blunt", "", top_boundary, bottom_boundary, ("blunt",))

    h = abs(stagger)
    if side == "left":
        # 左末端：上链更先进入片段 -> 上链 5' 突出；否则下链 3' 突出。
        if stagger > 0:
            overhang = _arc(sequence, top_boundary, h, circular, n)
            kind: EndKind = "5prime"
        else:
            overhang = reverse_complement_dna(
                _arc(sequence, top_boundary - h, h, circular, n)
            )
            kind = "3prime"
    else:
        # 右末端：下链延伸更远 -> 下链 5' 突出；否则上链 3' 突出。
        if stagger > 0:
            overhang = reverse_complement_dna(
                _arc(sequence, top_boundary, h, circular, n)
            )
            kind = "5prime"
        else:
            overhang = _arc(sequence, top_boundary - h, h, circular, n)
            kind = "3prime"

    canonical = min(overhang, reverse_complement_dna(overhang))
    return EndInfo(kind, overhang, top_boundary, bottom_boundary, (kind, canonical))


def digest(sequence: object, enzyme: object, *, circular: object = False) -> DigestResult:
    """对线性或环状双链 DNA 执行单酶消化。

    :param sequence: 上链序列（5'->3'），仅接受大写 A/C/G/T，长度 1..5000。
    :param enzyme: :class:`Enzyme` 实例。
    :param circular: True=环状（识别与切口均可跨零）；False=线性。
    :returns: :class:`DigestResult`。
    :raises TypeError: 类型不符（含 bool 冒充 int、float/NaN/Infinity）。
    :raises InvalidSequenceError: 序列内容或长度非法（带位置）。
    :raises InvalidEnzymeError: 酶定义非法（带字段名）。
    :raises UnsupportedDigestError: 交错切口套叠导致非简单片段。
    """
    # ---- 校验 ----------------------------------------------------------------
    seq = _require_str("sequence", sequence)
    if not isinstance(enzyme, Enzyme):
        raise TypeError(f"enzyme 必须是 Enzyme，收到 {type(enzyme).__name__}")
    if not isinstance(circular, bool):
        raise TypeError("circular 必须是 bool")
    n = len(seq)
    if n == 0:
        raise InvalidSequenceError("sequence 不得为空", None)
    if n > MAX_SEQUENCE_LEN:
        raise InvalidSequenceError(
            f"序列长度 {n} 超过上限 {MAX_SEQUENCE_LEN}", None
        )
    for i, ch in enumerate(seq):
        if ch not in _DNA_COMPLEMENT:
            raise InvalidSequenceError(
                f"序列位置 {i} 的字符 {ch!r} 不是大写 A/C/G/T", i
            )

    rec = enzyme.recognition
    L = len(rec)
    ot, ob = enzyme.top_offset, enzyme.bottom_offset
    d = ob - ot  # 恒定突出方向/长度

    forward_sets = tuple(IUPAC_BASES[ch] for ch in rec)
    reverse_sets = tuple(IUPAC_BASES[ch] for ch in reverse_complement_pattern(rec))

    # ---- 穷举扫描正/反向位点并转换切口 ---------------------------------------
    starts = range(n) if circular else range(n - L + 1)
    site_matches: list[SiteMatch] = []
    # (top_cut, bottom_cut) -> 命中来源列表
    raw_events: dict[tuple[int, int], list[tuple[Orientation, int]]] = {}

    for s in starts:
        # 环状时识别窗必须落在 L 个互不相同的碱基上；L>n 物理上不可能匹配
        # （不得把同一段碱基绕环读两遍当位点）。起点仍计入扫描预算。
        site_possible = not circular or L <= n
        if site_possible and _match_at(forward_sets, seq, s, circular, n):
            t, b = s + ot, s + ob
            if circular or (0 < t < n and 0 < b < n):
                if circular:
                    t, b = t % n, b % n
                key = (t, b)
                raw_events.setdefault(key, []).append(("forward", s))
                site_matches.append(SiteMatch(s, "forward"))
        if _match_at(reverse_sets, seq, s, circular, n):
            # 反向互补命中：识别链为下链，锚点映射到 j+L。
            t, b = s + L - ob, s + L - ot
            if circular or (0 < t < n and 0 < b < n):
                if circular:
                    t, b = t % n, b % n
                key = (t, b)
                raw_events.setdefault(key, []).append(("reverse", s))
                site_matches.append(SiteMatch(s, "reverse"))

    site_matches.sort(key=lambda m: (m.start, m.orientation))
    cut_events = tuple(
        CutEvent(t, b, tuple(sorted(contrib)))
        for (t, b), contrib in sorted(raw_events.items(), key=lambda kv: (kv[0][0], kv[0][1]))
    )
    k = len(cut_events)

    # ---- 简单片段拓扑检查 -----------------------------------------------------
    # 同一酶的所有事件满足 (b - t) 的环上差恒为 d mod n；因此环上/线性上
    # 相邻上链切口间距必须严格大于 |d|，否则两条链的切口互相套叠。
    if k > 0:
        ordered = sorted(cut_events, key=lambda e: e.top_cut)
        if circular:
            for i in range(k):
                t1 = ordered[i].top_cut
                t2 = ordered[(i + 1) % k].top_cut
                gap = n if k == 1 else (t2 - t1) % n
                if gap <= abs(d):
                    raise UnsupportedDigestError(
                        f"环状切口 {t1} 与下一上链切口 {t2} 的间距 {gap} "
                        f"<= 交错长度 {abs(d)}，交错切口套叠，非简单片段不受支持"
                    )
        else:
            for i in range(k - 1):
                gap = ordered[i + 1].top_cut - ordered[i].top_cut
                if gap <= abs(d):
                    raise UnsupportedDigestError(
                        f"线性上链切口 {ordered[i].top_cut} 与 "
                        f"{ordered[i + 1].top_cut} 的间距 {gap} <= 交错长度 "
                        f"{abs(d)}，交错切口套叠，非简单片段不受支持"
                    )

    # ---- 按上链边界构造双链片段 ----------------------------------------------
    fragments: list[Fragment] = []

    if k == 0:
        if circular:
            # 整环未切：无末端的环状片段（穷举证明无位点）。
            fragments.append(
                Fragment(
                    0, "circular", 0, 0, 0, 0,
                    seq, reverse_complement_dna(seq), None, None,
                )
            )
        else:
            end0 = EndInfo("blunt", "", 0, 0, ("blunt",))
            endn = EndInfo("blunt", "", n, n, ("blunt",))
            fragments.append(
                Fragment(
                    0, "linear", 0, n, 0, n,
                    seq, reverse_complement_dna(seq), end0, endn,
                )
            )
    elif circular:
        ordered = sorted(cut_events, key=lambda e: e.top_cut)
        for i, ev in enumerate(ordered):
            nxt = ordered[(i + 1) % k]
            top_len = n if k == 1 else (nxt.top_cut - ev.top_cut) % n
            bot_len = n if k == 1 else (nxt.bottom_cut - ev.bottom_cut) % n
            top_seq = _arc(seq, ev.top_cut, top_len, True, n)
            bot_positions = _arc(seq, ev.bottom_cut, bot_len, True, n)
            fragments.append(
                Fragment(
                    i, "linear",
                    ev.top_cut, nxt.top_cut,
                    ev.bottom_cut, nxt.bottom_cut,
                    top_seq,
                    reverse_complement_dna(bot_positions),
                    _make_end(ev.top_cut, ev.bottom_cut, "left", d, seq, True, n),
                    _make_end(nxt.top_cut, nxt.bottom_cut, "right", d, seq, True, n),
                )
            )
    else:
        ordered = sorted(cut_events, key=lambda e: e.top_cut)
        pairs = [(0, 0)] + [(e.top_cut, e.bottom_cut) for e in ordered] + [(n, n)]
        for i in range(len(pairs) - 1):
            a1, b1 = pairs[i]
            a2, b2 = pairs[i + 1]
            top_seq = _arc(seq, a1, a2 - a1, False, n)
            bot_positions = _arc(seq, b1, b2 - b1, False, n)
            fragments.append(
                Fragment(
                    i, "linear", a1, a2, b1, b2,
                    top_seq,
                    reverse_complement_dna(bot_positions),
                    _make_end(a1, b1, "left", d, seq, False, n),
                    _make_end(a2, b2, "right", d, seq, False, n),
                )
            )

    return DigestResult(
        sequence=seq,
        circular=circular,
        length=n,
        enzyme_name=enzyme.name,
        site_matches=tuple(site_matches),
        cut_events=tuple(sorted(cut_events, key=lambda e: e.top_cut)),
        fragments=tuple(fragments),
        positions_scanned=len(starts),
    )
