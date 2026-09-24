"""限制性内切酶的线性/环状双链 DNA 片段归属（标准库实现）。

本模块实现题目契约的有限数学子集，不宣称行业标准（如 REBASE）全兼容：

* 底物只接受具体碱基 ``ACGT``（大写）；IUPAC 简并只允许出现在酶识别串中。
* 切偏移以上链坐标表示，相对识别串起点；同时搜索正向与反向互补位点。
* 环状 DNA 坐标对长度取模；重复物理切口去重。
* 片段按上链切口排序（环状循环排序，线性追加分子两端 0/N 的平端边界）。
* 当交错切口导致片段无法做无交叉的简单归属时，抛出
  :class:`DigestUnsupportedError`，绝不猜测配对。
* 在此之上提供候选**环状切位图**与完全消化片段长度多重集的匹配
  （:class:`CutMap` / :func:`match_candidates`）：长度精确比较、
  旋转等价去重；切位图也可由 :func:`map_from_digest` 复用上面的
  旧切割规则生成。

全部计算均为整数运算，没有浮点输出，因此不存在误差口径/容差问题。
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from typing import Final, Iterable, Union

__all__ = [
    "MAX_BASES",
    "MAX_RECOGNITION_LEN",
    "TOPOLOGIES",
    "DigestError",
    "EnzymeDefinitionError",
    "SequenceError",
    "DigestUnsupportedError",
    "CandidateError",
    "Enzyme",
    "EndInfo",
    "Fragment",
    "CutSite",
    "DigestResult",
    "CutMap",
    "MatchResult",
    "compatible",
    "digest",
    "circular_fragment_lengths",
    "map_from_digest",
    "match_candidates",
]

MAX_BASES: Final[int] = 5000
MAX_RECOGNITION_LEN: Final[int] = 12
TOPOLOGIES: Final[tuple[str, ...]] = ("linear", "circular")

# 具体碱基（底物只允许这四个）。
_BASES: Final[frozenset[str]] = frozenset("ACGT")

# IUPAC 简并碱基 -> 可匹配的具体碱基集合。
_IUPAC: Final[dict[str, frozenset[str]]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),  # puRine
    "Y": frozenset("CT"),  # pYrimidine
    "S": frozenset("GC"),  # Strong
    "W": frozenset("AT"),  # Weak
    "K": frozenset("GT"),  # Keto
    "M": frozenset("AC"),  # aMino
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

_COMPLEMENT: Final[dict[str, str]] = {
    "A": "T", "C": "G", "G": "C", "T": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}


# ---------------------------------------------------------------------------
# 错误类型
# ---------------------------------------------------------------------------


class DigestError(ValueError):
    """本模块所有可定位输入/归属错误的基类。"""


class EnzymeDefinitionError(DigestError):
    """酶定义非法（识别串、偏移等）。"""


class SequenceError(DigestError):
    """底物序列非法（含可定位的下标信息）。"""


class DigestUnsupportedError(DigestError):
    """输入合法，但交错切口产生非简单片段，本契约明确拒绝猜测配对。"""


class CandidateError(DigestError):
    """候选切位图/片段长度多重集非法（含可定位的下标信息）。"""


# ---------------------------------------------------------------------------
# 校验小工具
# ---------------------------------------------------------------------------


def _check_int(value: object, what: str) -> int:
    """整数参数：拒绝 bool（bool 是 int 子类型），拒绝任何 float（含 NaN/Inf）。"""
    if isinstance(value, bool):
        raise EnzymeDefinitionError(f"{what} 不能是 bool 类型，必须是 int")
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        # 明确履行“数值输入拒绝 NaN/Infinity”：浮点一律不接受，
        # 对 NaN/Infinity 给出专门信息。
        if math.isnan(value):
            raise EnzymeDefinitionError(f"{what} 不能是 NaN，必须是有限整数")
        if math.isinf(value):
            raise EnzymeDefinitionError(f"{what} 不能是 Infinity，必须是有限整数")
        raise EnzymeDefinitionError(f"{what} 必须是 int，收到 float: {value!r}")
    raise EnzymeDefinitionError(f"{what} 必须是 int，收到 {type(value).__name__}")


def _check_str(value: object, what: str) -> str:
    if isinstance(value, bool) or not isinstance(value, str):
        raise EnzymeDefinitionError(f"{what} 必须是 str，收到 {type(value).__name__}")
    return value


# ---------------------------------------------------------------------------
# 公开数据结构
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Enzyme:
    """限制性内切酶定义。

    ``top_cut``/``bottom_cut`` 均为相对识别串起点的上链坐标边界：

    * 边界 k 表示位于第 k-1 与第 k 个碱基之间（k=0 在识别串之前，k=L 在其后）。
    * 正向命中起点 s：上链切口 s+top_cut，下链切口 s+bottom_cut。
    * 例：EcoRI ``GAATTC`` 为 G^AATTC / CTTAA^G，即 (1, 5)。
    """

    name: str
    recognition: str
    top_cut: int
    bottom_cut: int

    def __post_init__(self) -> None:
        name = _check_str(self.name, "name")
        if not name:
            raise EnzymeDefinitionError("name 不能为空")
        recognition = _check_str(self.recognition, "recognition")
        if not 1 <= len(recognition) <= MAX_RECOGNITION_LEN:
            raise EnzymeDefinitionError(
                f"recognition 长度必须在 1..{MAX_RECOGNITION_LEN} 之间，"
                f"收到 {len(recognition)}"
            )
        for i, ch in enumerate(recognition):
            if ch not in _IUPAC:
                raise EnzymeDefinitionError(
                    f"recognition 位置 {i} 的字符 {ch!r} 不是合法 IUPAC 碱基"
                    "（仅允许大写 ACGTRYSWKMBDHVN）"
                )
        top_cut = _check_int(self.top_cut, "top_cut")
        bottom_cut = _check_int(self.bottom_cut, "bottom_cut")

        # dataclass(frozen=True) 用 object.__setattr__ 完成归一化后的写回。
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "recognition", recognition)
        object.__setattr__(self, "top_cut", top_cut)
        object.__setattr__(self, "bottom_cut", bottom_cut)

    @property
    def stagger(self) -> int:
        """固有错位量 bottom_cut - top_cut：>0 为 5' 突出，<0 为 3' 突出。"""
        return self.bottom_cut - self.top_cut


@dataclass(frozen=True)
class EndInfo:
    """片段一个末端的信息。

    kind:
        * ``"blunt"``    平端，overhang 为 ""
        * ``"5_prime"``  5' 突出端
        * ``"3_prime"``  3' 突出端
    ``overhang`` 一律按 5'→3' 方向报告（平端为 ""），因此两个同类末端
    当且仅当 overhang 相同时可互补拼接。
    """

    kind: str
    top_boundary: int
    bottom_boundary: int
    overhang: str
    overhang_length: int


@dataclass(frozen=True)
class Fragment:
    """按上链切口分段得到的双链片段。

    坐标均沿上链 5'→3' 方向编号；环状时片段可跨越坐标 0。
    以 (起点, 长度) 表示，避免半开区间在跨零时歧义。

    * top_strand / bottom_strand 都按各自 5'→3' 报告（下链方向与上链相反）。
    * paired_length 为双链配对碱基数；两条链总长之差分布在两端突出上。
    * 未被切的环状分子：circular=True 且 left_end/right_end 为 None。
    """

    index: int
    topology: str
    circular: bool
    top_start: int
    top_length: int
    bottom_start: int
    bottom_length: int
    paired_length: int
    top_strand: str
    bottom_strand: str
    left_end: Union[EndInfo, None]
    right_end: Union[EndInfo, None]


@dataclass(frozen=True)
class CutSite:
    """一次识别命中（搜索证据；物理切口可能在多条命中间去重）。"""

    enzyme: str
    start: int
    strand: str  # "+" 正向，"-" 反向互补
    top_boundary: int
    bottom_boundary: int


@dataclass(frozen=True)
class DigestResult:
    topology: str
    length: int
    cut_sites: tuple[CutSite, ...]
    fragments: tuple[Fragment, ...]

    @property
    def cut_count(self) -> int:
        """去重后的物理双链断裂数量（按 (上链边界, 下链边界) 去重）。"""
        return len({(s.top_boundary, s.bottom_boundary) for s in self.cut_sites})


# ---------------------------------------------------------------------------
# 核心算法
# ---------------------------------------------------------------------------


def _validate_sequence(sequence: object) -> str:
    if isinstance(sequence, bool) or not isinstance(sequence, str):
        raise SequenceError(
            f"sequence 必须是 str，收到 {type(sequence).__name__}"
        )
    if not 1 <= len(sequence) <= MAX_BASES:
        raise SequenceError(
            f"sequence 长度必须在 1..{MAX_BASES} 之间，收到 {len(sequence)}"
        )
    for i, ch in enumerate(sequence):
        if ch not in _BASES:
            raise SequenceError(
                f"sequence 位置 {i} 的字符 {ch!r} 非法：底物只允许大写 ACGT"
                "（不做静默大小写/空白纠正，N 等简并碱基也不接受）"
            )
    return sequence


def _pattern_sets(pattern: str) -> tuple[frozenset[str], ...]:
    return tuple(_IUPAC[ch] for ch in pattern)


def _reverse_complement_iupac(pattern: str) -> str:
    return "".join(_COMPLEMENT[pattern[i]] for i in range(len(pattern) - 1, -1, -1))


def _scan(
    top: str,
    sets: tuple[frozenset[str], ...],
    circular: bool,
) -> list[int]:
    n = len(top)
    length = len(sets)
    last = n if circular else n - length + 1
    hits: list[int] = []
    for s in range(last):
        ok = True
        for i in range(length):
            idx = (s + i) % n if circular else s + i
            if top[idx] not in sets[i]:
                ok = False
                break
        if ok:
            hits.append(s)
    return hits


def _arc(seq: str, start: int, length: int, circular: bool) -> str:
    """沿上链坐标正向取 length 个碱基；环状可跨 0，线性不跨（由构造保证）。"""
    if length == 0:
        return ""
    n = len(seq)
    if circular:
        doubled = seq + seq
        return doubled[start:start + length]
    return seq[start:start + length]


def _revcomp_concrete(s: str) -> str:
    return "".join(_COMPLEMENT[ch] for ch in reversed(s))


def _make_end(
    kind_side: str,
    t: int,
    b: int,
    delta: int,
    top: str,
    circular: bool,
) -> EndInfo:
    """根据 (上链边界 t, 下链边界 b, 固有错位 delta) 构造左/右末端。"""
    if delta == 0:
        return EndInfo("blunt", t, b, "", 0)
    kind = "5_prime" if delta > 0 else "3_prime"
    if kind_side == "left":
        if delta > 0:
            # 上链 5' 端多出来：top[t:b]
            overhang = _arc(top, t, delta, circular)
        else:
            # 下链 3' 端多出来，5'→3' 报告 = revcomp(top[b:t])
            overhang = _revcomp_concrete(_arc(top, b, -delta, circular))
    else:  # right
        if delta > 0:
            # 下链 5' 端多出来，5'→3' 报告 = revcomp(top[t:b])
            overhang = _revcomp_concrete(_arc(top, t, delta, circular))
        else:
            # 上链 3' 端多出来：top[b:t]
            overhang = _arc(top, b, -delta, circular)
    return EndInfo(kind, t, b, overhang, abs(delta))


def compatible(left: Union[EndInfo, None], right: Union[EndInfo, None]) -> bool:
    """两个末端是否可互补拼接。

    突出序列一律按各自突出链 5'→3' 报告；可配对接合处两条突出方向相反，
    因此序列互为反向互补（平端配平端）。EcoRI 的 AATT 自补，所以看起来
    “相同即可”，但非回文 4 碱基（如 BsaI 产生的端）必须按反向互补判定。

    未切环状分子的 None 末端不可与任何末端拼接。
    """
    if left is None or right is None:
        return False
    if left.kind == "blunt":
        return right.kind == "blunt"
    if right.kind != left.kind or right.overhang_length != left.overhang_length:
        return False
    return _revcomp_concrete(right.overhang) == left.overhang


# 内部用的“物理末端”：去重后的 (上链边界, 下链边界, 固有错位)。
_PhysEnd = tuple[int, int, int]


def digest(
    sequence: str,
    enzymes: Union[Enzyme, Iterable[Enzyme]],
    topology: str = "linear",
) -> DigestResult:
    """对线性/环状双链 DNA 做一个或多个限制酶的完全酶切片段归属。

    纯函数：不修改任何输入，也没有内部状态，多酶批次要么整体返回结果，
    要么在任何片段构造前抛异常，不存在部分变更。

    搜索为穷举（O(N·L)，N≤5000，L≤12）：“没有位点”是被证明的结果
    （返回完整分子的单个片段），不存在搜索预算耗尽这一失败路径。
    """
    top = _validate_sequence(sequence)
    if isinstance(topology, bool) or not isinstance(topology, str) or topology not in TOPOLOGIES:
        raise DigestError(
            f"topology 必须是 {TOPOLOGIES} 之一，收到 {topology!r}"
        )
    circular = topology == "circular"
    n = len(top)

    if isinstance(enzymes, Enzyme):
        enzyme_list = (enzymes,)
    else:
        try:
            enzyme_list = tuple(enzymes)
        except TypeError as exc:
            raise EnzymeDefinitionError("enzymes 必须是 Enzyme 或 Enzyme 可迭代对象") from exc
    if not enzyme_list:
        raise EnzymeDefinitionError("至少需要一个 Enzyme")
    for i, enz in enumerate(enzyme_list):
        if not isinstance(enz, Enzyme):
            raise EnzymeDefinitionError(
                f"enzymes[{i}] 不是 Enzyme 实例（类型 {type(enz).__name__}）"
            )

    # ---- 1. 搜索正/反向位点并转换切口坐标 ----
    raw_sites: list[CutSite] = []
    # 物理末端键 (t mod n, b mod n) -> 固有 delta（环状必须唯一，否则歧义）
    phys: dict[tuple[int, int], int] = {}
    # 上/下链切口坐标 -> 对侧坐标，检测“同一物理切口、不同伴侣”的冲突。
    top_seen: dict[int, int] = {}
    bot_seen: dict[int, int] = {}

    def register(t_raw: int, b_raw: int, delta: int, where: str) -> Union[tuple[int, int], None]:
        """登记一次双链断裂；返回去重后的 (上链边界, 下链边界)。

        线性分子的边缘语义：
        * (0,0) 或 (n,n)：切口与天然末端重合，不产生实际切割（忽略登记）。
        * 仅一条链切在 0/n 而另一条在内部：会产生单链片段，unsupported。
        * 切口落在 [0,n] 之外：酶在分子外下刀，unsupported。
        """
        if circular:
            # 错位量相对整圆不可区分时无法确定端型，明确拒绝。
            if abs(delta) >= n:
                raise DigestUnsupportedError(
                    f"{where} 的交错量 |{delta}| 不小于环长 {n}，"
                    "取模后无法确定上/下链切口的配对关系"
                )
            t = t_raw % n
            b = b_raw % n
        else:
            if not (0 <= t_raw <= n and 0 <= b_raw <= n):
                raise DigestUnsupportedError(
                    f"{where} 的切口坐标 (上链 {t_raw}, 下链 {b_raw}) 超出线性分子 "
                    f"[0, {n}]：酶切点落在分子之外，无法形成双链断裂片段"
                )
            t, b = t_raw, b_raw
            if (t, b) in ((0, 0), (n, n)):
                return None  # 平端切口与天然末端重合，无实际切割
            if t in (0, n) or b in (0, n):
                raise DigestUnsupportedError(
                    f"{where} 的一条链切口落在分子末端（上链 {t}, 下链 {b}），"
                    "另一条链切口却在分子内部：会切出单链片段，"
                    "不属于本契约支持的简单双链片段"
                )
        key = (t, b)
        if key in phys:
            if phys[key] != delta:
                raise DigestUnsupportedError(
                    f"{where} 与既有切口在坐标 ({t}, {b}) 重合但错位量不同"
                    f"（{delta} vs {phys[key]}），末端归属有歧义"
                )
        else:
            if t in top_seen and top_seen[t] != b:
                raise DigestUnsupportedError(
                    f"{where} 的上链切口 {t} 已与下链切口 {top_seen[t]} 配对，"
                    f"又出现配对下链切口 {b}，重复物理切口无法去重为同一末端"
                )
            if b in bot_seen and bot_seen[b] != t:
                raise DigestUnsupportedError(
                    f"{where} 的下链切口 {b} 已与上链切口 {bot_seen[b]} 配对，"
                    f"又出现配对上链切口 {t}，重复物理切口无法去重为同一末端"
                )
            phys[key] = delta
            top_seen[t] = b
            bot_seen[b] = t
        return t, b

    for enz in enzyme_list:
        length = len(enz.recognition)
        fwd = _pattern_sets(enz.recognition)
        rc_pattern = _reverse_complement_iupac(enz.recognition)
        # 回文识别串（如 GAATTC）在每个位置的反向命中就是同一次物理断裂，
        # 不能重复计数；只有非回文酶才需要再扫反向互补。
        scan_reverse = rc_pattern != enz.recognition
        rev = _pattern_sets(rc_pattern)
        a = enz.top_cut
        c = enz.bottom_cut
        for s in _scan(top, fwd, circular):
            reg = register(s + a, s + c, c - a,
                           f"酶 {enz.name!r} 正向位点 start={s}")
            if reg is not None:
                t, b = reg
                raw_sites.append(CutSite(enz.name, s, "+", t, b))
        if scan_reverse:
            for s in _scan(top, rev, circular):
                # 反向互补命中：酶旋转 180°，上链 = s+L-b，下链 = s+L-a。
                t_raw = s + length - c
                b_raw = s + length - a
                reg = register(t_raw, b_raw, c - a,
                               f"酶 {enz.name!r} 反向互补位点 start={s}")
                if reg is not None:
                    t, b = reg
                    raw_sites.append(CutSite(enz.name, s, "-", t, b))

    raw_sites.sort(key=lambda x: (x.start, x.strand, x.enzyme))
    cut_sites = tuple(raw_sites)
    ends: list[_PhysEnd] = [(t, b, d) for (t, b), d in phys.items()]
    ends.sort(key=lambda e: e[0])

    # ---- 2. 未切环状分子：一个无末端的环状片段 ----
    if circular and not ends:
        bottom = _revcomp_concrete(top)  # 5'→3'
        frag = Fragment(
            index=0, topology=topology, circular=True,
            top_start=0, top_length=n,
            bottom_start=0, bottom_length=n,
            paired_length=n,
            top_strand=top, bottom_strand=bottom,
            left_end=None, right_end=None,
        )
        return DigestResult(topology, n, cut_sites, (frag,))

    # ---- 3. 末端排序列表（线性追加 0/N 平端分子末端）----
    if circular:
        ordered: list[_PhysEnd] = ends
    else:
        ordered = [(0, 0, 0), *ends, (n, n, 0)]
        # 若酶切恰好落在分子末端且为平端，上面 register 已去重；
        # 这里再按 t 排序并做一次严格性检查。
        ordered.sort(key=lambda e: e[0])
        deduped: list[_PhysEnd] = []
        for e in ordered:
            if deduped and deduped[-1][0] == e[0]:
                if deduped[-1] != e:
                    raise DigestUnsupportedError(
                        f"上链切口坐标 {e[0]} 与分子末端重合但下链边界不同"
                        f"（{deduped[-1][1]} vs {e[1]}），末端归属有歧义"
                    )
                continue
            deduped.append(e)
        ordered = deduped

    k = len(ordered)

    # ---- 4. 连续上链切口之间构造片段，逐条做简单性判定 ----
    fragments: list[Fragment] = []
    for i in range(k if circular else k - 1):
        t0, b0, d0 = ordered[i]
        t1, b1, d1 = ordered[(i + 1) % k]

        top_len = (t1 - t0) % n if circular else t1 - t0
        if circular and i == k - 1 and k == 1:
            top_len = n
        if top_len <= 0:
            raise DigestUnsupportedError(
                f"上链切口 {t0} 与 {t1} 之间距离为 0，无法形成片段"
            )
        # 下链正向跨度（局部展开坐标）。
        bot_len = top_len - d0 + d1
        if circular:
            if not 0 < bot_len <= n:
                raise DigestUnsupportedError(
                    f"上链切口 {t0}→{t1} 之间下链边界交错（下链跨度 {bot_len}），"
                    "交错切口相互重叠/交叉，无法做简单片段归属"
                )
        elif bot_len <= 0:
            raise DigestUnsupportedError(
                f"上链切口 {t0}→{t1} 之间下链跨度 {bot_len}，"
                "交错切口相互重叠，无法做简单片段归属"
            )

        # 配对区 = 两端边界的局部交集：[max(0,d0), top_len+min(0,d1)]。
        paired = top_len + min(0, d1) - max(0, d0)
        if paired <= 0:
            raise DigestUnsupportedError(
                f"上链切口 {t0}→{t1}（间距 {top_len}）容不下两端的突出"
                f"（错位 {d0}, {d1}）：配对区长度 {paired}，"
                "交错切口互相侵入或使双链完全解离，无法做简单片段归属"
            )

        top_strand = _arc(top, t0, top_len, circular)
        bottom_aligned = _arc(top, b0, bot_len, circular)
        bottom_strand = _revcomp_concrete(bottom_aligned)

        left = _make_end("left", t0, b0, d0, top, circular)
        right = _make_end("right", t1, b1, d1, top, circular)

        fragments.append(Fragment(
            index=i, topology=topology,
            circular=False,
            top_start=t0, top_length=top_len,
            bottom_start=b0, bottom_length=bot_len,
            paired_length=paired,
            top_strand=top_strand,
            bottom_strand=bottom_strand,
            left_end=left, right_end=right,
        ))

    return DigestResult(topology, n, cut_sites, tuple(fragments))


# ---------------------------------------------------------------------------
# 候选环状切位图 × 完全消化片段长度多重集
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CutMap:
    """一个候选环状切位图：环长 ``length`` 与升序、不重复的切口边界。

    边界 k 位于第 k-1 与第 k 个碱基之间（与 :class:`Enzyme` 的切偏移
    约定一致），取值 ``0..length-1``；环上坐标只相差旋转时视为同一张图。

    * 0 个切口表示**未切环**：完全消化后仍是一整个环（一个长度 length
      的环状片段）。
    * 1 个切口表示**单切环**：完全消化后得到一个长度 length 的线性片段。
      两者长度多重集都是 ``(length,)``，但切口数不同，因此互不等价、
      分别作为候选保留（去重键含切口数）。
    """

    length: int
    cuts: tuple[int, ...]
    label: str = ""

    def __post_init__(self) -> None:
        length = _check_map_length(self.length)
        cuts, label = _check_map_parts(self.cuts, self.label, length)
        object.__setattr__(self, "length", length)
        object.__setattr__(self, "cuts", cuts)
        object.__setattr__(self, "label", label)

    @property
    def cut_count(self) -> int:
        return len(self.cuts)

    @property
    def is_uncut(self) -> bool:
        """是否为未切环（0 个切口）。"""
        return not self.cuts

    def fragment_lengths(self) -> tuple[int, ...]:
        """完全消化后的片段长度元组（按切口升序的循环序，非多重集排序）。"""
        return circular_fragment_lengths(self.length, self.cuts)


@dataclass(frozen=True)
class MatchResult:
    """候选匹配结果。

    * ``matches``：与片段长度多重集匹配的候选，已做**旋转等价去重**，
      保留每个等价类中先出现的候选（输入顺序）。
    * ``examined``：去重前检查过的候选总数。
    * ``unique_classes``：参与匹配的全部候选（含不匹配的）按旋转归一化
      后的等价类数量，即去重后的候选数。
    """

    matches: tuple[CutMap, ...]
    examined: int
    unique_classes: int

    @property
    def matched(self) -> bool:
        """是否至少存在一个匹配候选。"""
        return bool(self.matches)


def _check_map_length(value: object) -> int:
    """切位图环长：拒绝 bool/float（含 NaN/Inf），要求 >= 1。"""
    if isinstance(value, bool):
        raise CandidateError("length 不能是 bool 类型，必须是 int")
    if isinstance(value, int):
        if value < 1:
            raise CandidateError(f"length 必须 >= 1，收到 {value}")
        return value
    if isinstance(value, float):
        if math.isnan(value):
            raise CandidateError("length 不能是 NaN，必须是有限正整数")
        if math.isinf(value):
            raise CandidateError("length 不能是 Infinity，必须是有限正整数")
        raise CandidateError(f"length 必须是 int，收到 float: {value!r}")
    raise CandidateError(f"length 必须是 int，收到 {type(value).__name__}")


def _check_map_parts(
    cuts: object, label: object, length: int
) -> tuple[tuple[int, ...], str]:
    if isinstance(label, bool) or not isinstance(label, str):
        raise CandidateError(f"label 必须是 str，收到 {type(label).__name__}")
    if isinstance(cuts, (str, bytes)) or not isinstance(cuts, Iterable):
        raise CandidateError(
            f"cuts 必须是 int 可迭代对象，收到 {type(cuts).__name__}"
        )
    norm: list[int] = []
    for i, c in enumerate(cuts):
        if isinstance(c, bool):
            raise CandidateError(f"cuts[{i}] 不能是 bool，必须是 int")
        if not isinstance(c, int):
            if isinstance(c, float):
                raise CandidateError(f"cuts[{i}] 必须是 int，收到 float: {c!r}")
            raise CandidateError(
                f"cuts[{i}] 必须是 int，收到 {type(c).__name__}"
            )
        if not 0 <= c < length:
            raise CandidateError(
                f"cuts[{i}] = {c} 超出环长 {length} 的切口边界范围 "
                f"[0, {length - 1}]"
            )
        if norm and c <= norm[-1]:
            # 升序且不重复；显式报错便于定位（不静默排序/去重）。
            raise CandidateError(
                f"cuts[{i}] = {c} 不满足严格升序/不重复"
                f"（前一个切口为 {norm[-1]}）；切口必须显式、升序、不重复"
            )
        norm.append(c)
    return tuple(norm), label


def circular_fragment_lengths(
    length: int, cuts: Iterable[int]
) -> tuple[int, ...]:
    """环长 ``length``、切口升序 ``cuts`` 的完全消化片段长度元组。

    切口 p0 < p1 < ... < p(k-1) 把环分成 k 段，循环序长度为
    ``(p1-p0, p2-p1, ..., length-p(k-1)+p0)``；0 个切口时整环未切，
    返回 ``(length,)``。长度为精确整数，调用方做多重集比较。
    """
    n = _check_map_length(length)
    pts, _ = _check_map_parts(cuts, "", n)
    if not pts:
        return (n,)
    gaps = [pts[i + 1] - pts[i] for i in range(len(pts) - 1)]
    gaps.append(n - pts[-1] + pts[0])
    return tuple(gaps)


def _rotation_key(n: int, cuts: tuple[int, ...]) -> tuple[int, ...]:
    """旋转等价归一化键：间隙循环序列的字典序最小移位（前面补切口数）。

    环旋转只改变切口的绝对坐标，不改变相邻切口的循环间隙序列；
    补上切口数是为了让“未切环”和“单切环”（间隙同为 (n,)）不会被
    错误地归为同一个等价类。
    """
    if not cuts:
        return (0,)
    gaps = circular_fragment_lengths(n, cuts)
    k = len(gaps)
    best = gaps
    for i in range(1, k):
        rot = gaps[i:] + gaps[:i]
        if rot < best:
            best = rot
    return (k, *best)


def map_from_digest(result: DigestResult) -> CutMap:
    """从环状 :func:`digest` 结果提取候选切位图（复用旧切割规则）。

    取去重后的物理断裂上链边界（模环长）；线性 ``DigestResult`` 不适用
    （本功能只处理环状候选），分子天然末端不是酶切口。
    """
    if not isinstance(result, DigestResult):
        raise CandidateError(
            f"result 必须是 DigestResult，收到 {type(result).__name__}"
        )
    if result.topology != "circular":
        raise CandidateError(
            f"map_from_digest 只接受环状 DigestResult，收到 topology="
            f"{result.topology!r}"
        )
    cuts = tuple(sorted({s.top_boundary for s in result.cut_sites}))
    return CutMap(result.length, cuts)


def match_candidates(
    observed: Iterable[int], candidates: Iterable[CutMap]
) -> MatchResult:
    """在候选环状切位图中找出与完全消化片段长度多重集匹配者。

    :param observed: 完全消化测得的片段长度多重集（顺序无关，可重复），
        每个长度必须是正整数。
    :param candidates: 候选 :class:`CutMap` 可迭代对象；先全部校验、
        归一化，再开始匹配，任何一个非法都整体抛出 :class:`CandidateError`，
        不产生半更新结果（本函数无内部状态，失败即无返回值）。
    :returns: :class:`MatchResult`；匹配标准是长度多重集**精确相等**
        （``Counter`` 比较，不容差），并对候选做旋转等价去重：
        同一旋转等价类只保留先出现者；未切环与单切环不作为等价类合并。

    穷举小规模候选，无预算/超时路径；“无匹配”是被证明的结果。
    """
    if isinstance(observed, (str, bytes)) or not isinstance(observed, Iterable):
        raise CandidateError(
            f"observed 必须是 int 可迭代对象，收到 {type(observed).__name__}"
        )
    target_counter: Counter[int] = Counter()
    total = 0
    for i, x in enumerate(observed):
        if isinstance(x, bool):
            raise CandidateError(f"observed[{i}] 不能是 bool，必须是正整数")
        if not isinstance(x, int):
            if isinstance(x, float):
                if math.isnan(x):
                    raise CandidateError(f"observed[{i}] 不能是 NaN")
                if math.isinf(x):
                    raise CandidateError(f"observed[{i}] 不能是 Infinity")
                raise CandidateError(
                    f"observed[{i}] 必须是 int，收到 float: {x!r}"
                )
            raise CandidateError(
                f"observed[{i}] 必须是 int，收到 {type(x).__name__}"
            )
        if x <= 0:
            raise CandidateError(
                f"observed[{i}] = {x} 不是正整数：片段长度必须 >= 1"
            )
        target_counter[x] += 1
        total += x
    if not target_counter:
        raise CandidateError(
            "observed 不能为空：完全消化至少产生一个片段"
            "（未切环请给 [环长]）"
        )

    if isinstance(candidates, (str, bytes)) or not isinstance(candidates, Iterable):
        raise CandidateError(
            f"candidates 必须是 CutMap 可迭代对象，收到 {type(candidates).__name__}"
        )
    # 先物化并让每个 CutMap 完成自身校验（构造时即校验），保证“全部校验
    # 通过后才匹配”的原子语义；惰性生成器也在此刻被完整消费一次。
    try:
        cand_list = list(candidates)
    except CandidateError:
        raise
    except TypeError as exc:
        raise CandidateError(f"candidates 迭代失败：{exc}") from exc
    for i, cand in enumerate(cand_list):
        if not isinstance(cand, CutMap):
            raise CandidateError(
                f"candidates[{i}] 不是 CutMap 实例（类型 {type(cand).__name__}）"
            )

    seen: set[tuple[int, ...]] = set()
    matches: list[CutMap] = []
    examined = 0
    for cand in cand_list:
        examined += 1
        key = _rotation_key(cand.length, cand.cuts)
        if key in seen:
            continue  # 与先前候选旋转等价：去重，不重复列入匹配
        seen.add(key)
        # Counter 相等已蕴含长度总和（=环长）相等；total 仅用于早筛。
        if cand.length == total and Counter(cand.fragment_lengths()) == target_counter:
            matches.append(cand)

    return MatchResult(tuple(matches), examined, len(seen))
