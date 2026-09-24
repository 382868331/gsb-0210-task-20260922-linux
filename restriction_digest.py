"""限制性内切酶的线性/环状双链 DNA 片段归属（标准库实现）。

本模块实现题目契约的有限数学子集，不宣称行业标准（如 REBASE）全兼容：

* 底物只接受具体碱基 ``ACGT``（大写）；IUPAC 简并只允许出现在酶识别串中。
* 切偏移以上链坐标表示，相对识别串起点；同时搜索正向与反向互补位点。
* 环状 DNA 坐标对长度取模；重复物理切口去重。
* 片段按上链切口排序（环状循环排序，线性追加分子两端 0/N 的平端边界）。
* 当交错切口导致片段无法做无交叉的简单归属时，抛出
  :class:`DigestUnsupportedError`，绝不猜测配对。
* :func:`match_candidates` 在旧切割规则之上做候选环状切位图筛查：
  输入一组候选环状序列与完全消化的片段长度多重集，精确匹配并
  按旋转等价类去重后列出命中候选。

全部计算均为整数运算，没有浮点输出，因此不存在误差口径/容差问题。
"""

from __future__ import annotations

import math
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
    "FragmentLengthsError",
    "Enzyme",
    "EndInfo",
    "Fragment",
    "CutSite",
    "DigestResult",
    "CandidateMatch",
    "UnsupportedCandidate",
    "ScreenResult",
    "compatible",
    "digest",
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


class FragmentLengthsError(DigestError):
    """观测的完全消化片段长度多重集非法（空、非正整数、bool/浮点等）。"""


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
# 候选环状切位图筛查（调用上面的旧切割规则）
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CandidateMatch:
    """一个与观测片段长度多重集精确匹配的候选旋转等价类。

    * representative：输入中该类首次出现的候选序列（同一物理环状分子的
      不同旋转写法只报告一次）。
    * indices：该等价类在输入候选中的全部下标（按输入顺序）。
    * fragment_lengths：完全消化后的片段长度多重集（升序）。片段长度取
      相邻上链切口之间的碱基跨度（酶切图谱的通行约定，总和恒等于环长）。
    * cut_count：去重后的物理断裂数；0 表示未切环。注意未切环与单切环
      的长度多重集都是 ``(N,)``，只能靠 cut_count 区分。
    """

    representative: str
    indices: tuple[int, ...]
    length: int
    fragment_lengths: tuple[int, ...]
    cut_count: int


@dataclass(frozen=True)
class UnsupportedCandidate:
    """几何上无法做简单片段归属的候选类（不视为匹配，也不中断筛查）。"""

    representative: str
    indices: tuple[int, ...]
    reason: str


@dataclass(frozen=True)
class ScreenResult:
    """候选筛查结果。

    * target_lengths：校验并升序排序后的观测片段长度多重集。
    * candidate_count：输入候选总数（含旋转重复）。
    * class_count：旋转等价去重后实际评估的候选类数。
    * matches：长度多重集与观测精确相等（含重复计数）的候选类，
      按输入中首次出现的顺序排列。
    * unsupported：digest 判定为非简单几何（DigestUnsupportedError）的
      候选类，附原因；它们不参与匹配。
    """

    target_lengths: tuple[int, ...]
    candidate_count: int
    class_count: int
    matches: tuple[CandidateMatch, ...]
    unsupported: tuple[UnsupportedCandidate, ...]


def _validate_fragment_lengths(value: object) -> tuple[int, ...]:
    """校验观测片段长度多重集：非空、元素为正整数；拒绝 bool 与任何 float。"""
    if isinstance(value, bool) or isinstance(value, (str, bytes)):
        raise FragmentLengthsError(
            f"fragment_lengths 必须是正整数的可迭代对象，收到 {type(value).__name__}"
        )
    try:
        items = tuple(value)  # type: ignore[call-overload]
    except TypeError:
        raise FragmentLengthsError(
            f"fragment_lengths 必须是正整数的可迭代对象，收到 {type(value).__name__}"
        ) from None
    if not items:
        raise FragmentLengthsError(
            "fragment_lengths 不能为空：完全消化至少产生一个片段"
            "（未切环为全长单片段）"
        )
    out: list[int] = []
    for i, v in enumerate(items):
        if isinstance(v, bool):
            raise FragmentLengthsError(
                f"fragment_lengths[{i}] 不能是 bool，必须是正整数"
            )
        if isinstance(v, float):
            if math.isnan(v):
                raise FragmentLengthsError(
                    f"fragment_lengths[{i}] 不能是 NaN，长度必须是精确正整数"
                )
            if math.isinf(v):
                raise FragmentLengthsError(
                    f"fragment_lengths[{i}] 不能是 Infinity，长度必须是精确正整数"
                )
            raise FragmentLengthsError(
                f"fragment_lengths[{i}] 必须是 int，收到 float: {v!r}"
                "（本契约只做精确长度匹配，不接受测量近似值）"
            )
        if not isinstance(v, int):
            raise FragmentLengthsError(
                f"fragment_lengths[{i}] 必须是 int，收到 {type(v).__name__}"
            )
        if v < 1:
            raise FragmentLengthsError(
                f"fragment_lengths[{i}] 必须是正整数，收到 {v}"
            )
        out.append(v)
    return tuple(sorted(out))


def _validate_candidates(value: object) -> tuple[str, ...]:
    """校验候选序列集合；错误信息带候选下标，便于定位。"""
    if isinstance(value, bool) or isinstance(value, (str, bytes)):
        raise SequenceError(
            "candidates 必须是候选序列（str）的可迭代对象，不是单个字符串"
        )
    try:
        items = tuple(value)  # type: ignore[call-overload]
    except TypeError:
        raise SequenceError(
            f"candidates 必须是候选序列的可迭代对象，收到 {type(value).__name__}"
        ) from None
    out: list[str] = []
    for i, cand in enumerate(items):
        try:
            out.append(_validate_sequence(cand))
        except SequenceError as exc:
            raise SequenceError(f"candidates[{i}] 非法：{exc}") from exc
    return tuple(out)


def _canonical_rotation(s: str) -> str:
    """最小字典序旋转（Booth 算法，O(n)），作为旋转等价类的规范键。"""
    n = len(s)
    doubled = s + s
    i, j, k = 0, 1, 0
    while i < n and j < n and k < n:
        a = doubled[i + k]
        b = doubled[j + k]
        if a == b:
            k += 1
        elif a > b:
            i = i + k + 1
            if i <= j:
                i = j + 1
            k = 0
        else:
            j = j + k + 1
            if j <= i:
                j = i + 1
            k = 0
    start = min(i, j)
    return doubled[start:start + n]


def match_candidates(
    candidates: Iterable[str],
    enzymes: Union[Enzyme, Iterable[Enzyme]],
    fragment_lengths: Iterable[int],
) -> ScreenResult:
    """用旧切割规则筛查候选环状切位图，列出与观测片段长度匹配的候选。

    对每个候选环状序列调用 :func:`digest`（``topology="circular"``）做
    完全消化，取各片段的上链切口跨度作为片段长度，与观测的
    ``fragment_lengths`` 做多重集精确比较（重复片段按计数匹配）。

    * 同一物理环状分子的不同旋转写法是同一个候选：按最小字典序旋转
      归一去重，每个等价类只评估一次、只报告一次（附全部输入下标）。
    * digest 判定为非简单几何（:class:`DigestUnsupportedError`）的候选类
      不视为匹配，收入 ``ScreenResult.unsupported``，不中断其余候选。
    * 纯函数：候选与片段长度在任何评估之前完成校验，非法输入立即抛出
      （:class:`SequenceError` / :class:`FragmentLengthsError` /
      :class:`EnzymeDefinitionError`），不存在部分结果。
    * 长度精确匹配，不做测量区间/容差反演；未切环的长度多重集为 ``(N,)``，
      与单切环不可区分（用 ``CandidateMatch.cut_count`` 辨别）。
    """
    target = _validate_fragment_lengths(fragment_lengths)
    candidate_list = _validate_candidates(candidates)

    # 酶列表预先归一化（避免可迭代对象在多轮 digest 中被耗尽），
    # 并保证酶定义错误在任何候选评估之前抛出。
    if isinstance(enzymes, Enzyme):
        enzyme_tuple = (enzymes,)
    else:
        try:
            enzyme_tuple = tuple(enzymes)
        except TypeError as exc:
            raise EnzymeDefinitionError(
                "enzymes 必须是 Enzyme 或 Enzyme 可迭代对象"
            ) from exc
    if not enzyme_tuple:
        raise EnzymeDefinitionError("至少需要一个 Enzyme")
    for i, enz in enumerate(enzyme_tuple):
        if not isinstance(enz, Enzyme):
            raise EnzymeDefinitionError(
                f"enzymes[{i}] 不是 Enzyme 实例（类型 {type(enz).__name__}）"
            )

    # 旋转等价去重：规范键 -> 该类首次出现的序列与全部输入下标。
    reps: dict[str, str] = {}
    members: dict[str, list[int]] = {}
    for i, seq in enumerate(candidate_list):
        key = _canonical_rotation(seq)
        if key in reps:
            members[key].append(i)
        else:
            reps[key] = seq
            members[key] = [i]

    matches: list[CandidateMatch] = []
    unsupported: list[UnsupportedCandidate] = []
    for key, rep in reps.items():
        indices = tuple(members[key])
        try:
            res = digest(rep, enzyme_tuple, topology="circular")
        except DigestUnsupportedError as exc:
            unsupported.append(UnsupportedCandidate(rep, indices, str(exc)))
            continue
        lengths = tuple(sorted(f.top_length for f in res.fragments))
        if lengths == target:
            matches.append(CandidateMatch(
                representative=rep,
                indices=indices,
                length=len(rep),
                fragment_lengths=lengths,
                cut_count=res.cut_count,
            ))

    return ScreenResult(
        target_lengths=target,
        candidate_count=len(candidate_list),
        class_count=len(reps),
        matches=tuple(matches),
        unsupported=tuple(unsupported),
    )
