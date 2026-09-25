"""限制性内切酶的线性/环状双链 DNA 片段归属（标准库实现）。

本模块实现题目契约的有限数学子集，不宣称行业标准（如 REBASE）全兼容：

* 底物只接受具体碱基 ``ACGT``（大写）；IUPAC 简并只允许出现在酶识别串中。
* 切偏移以上链坐标表示，相对识别串起点；同时搜索正向与反向互补位点。
* :func:`recognize_sites` 先返回带稳定 ID 的物理识别位点；:func:`digest` 可只启用子集。
* 环状 DNA 坐标对长度取模；重复物理切口去重。
* 片段按上链切口排序（环状循环排序，线性追加分子两端 0/N 的平端边界）。
* 当交错切口导致片段无法做无交叉的简单归属时，抛出
  :class:`DigestUnsupportedError`，绝不猜测配对。

全部计算均为整数运算，没有浮点输出，因此不存在误差口径/容差问题。
"""

from __future__ import annotations

import math
import hashlib
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
    "UnknownSiteError",
    "Enzyme",
    "EndInfo",
    "Fragment",
    "CutSite",
    "RecognitionSite",
    "DigestResult",
    "compatible",
    "recognize_sites",
    "digest",
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


class UnknownSiteError(DigestError):
    """启用的位点 ID 不属于本次 DNA、拓扑与酶表的识别结果。"""


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
    """一个已启用识别位点产生的一条物理双链断裂证据。

    回文识别串的正/反向来源合并为一条记录，方向保存在 ``sources`` 中；
    两个不同酶即使产生完全相同切口，也仍保留为各自的证据。
    """

    enzyme: str
    start: int
    strand: str  # 主方向："+" 正向，"-" 反向互补
    top_boundary: int
    bottom_boundary: int
    site_id: str = ""
    sources: tuple[str, ...] = ("+",)


@dataclass(frozen=True)
class RecognitionSite:
    """识别阶段的一个物理识别位点，尚未因为该位点执行切口登记。"""

    site_id: str
    enzyme_index: int
    enzyme: str
    start: int
    length: int
    matched_sequence: str
    sources: tuple[str, ...]


@dataclass(frozen=True)
class DigestResult:
    topology: str
    length: int
    cut_sites: tuple[CutSite, ...]
    fragments: tuple[Fragment, ...]
    recognition_sites: tuple[RecognitionSite, ...] = ()

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


def _normalize_enzymes(
    enzymes: Union[Enzyme, Iterable[Enzyme]],
) -> tuple[Enzyme, ...]:
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
    return enzyme_list


def _input_fingerprint(
    top: str,
    topology: str,
    enzymes: tuple[Enzyme, ...],
) -> str:
    """生成绑定本次 DNA、拓扑和有序酶表的短指纹，防止复用旧结果 ID。"""
    enzyme_spec = "|".join(
        f"{i}:{enz.name}:{enz.recognition}:{enz.top_cut}:{enz.bottom_cut}"
        for i, enz in enumerate(enzymes)
    )
    payload = f"{topology}\n{len(top)}\n{top}\n{enzyme_spec}".encode("ascii")
    return hashlib.sha256(payload).hexdigest()[:12]


def _recognize(
    top: str,
    topology: str,
    enzymes: tuple[Enzyme, ...],
) -> tuple[RecognitionSite, ...]:
    circular = topology == "circular"
    fingerprint = _input_fingerprint(top, topology, enzymes)
    grouped: dict[tuple[int, int], dict[str, None]] = {}
    metadata: dict[tuple[int, int], tuple[int, int, str]] = {}

    for enzyme_index, enz in enumerate(enzymes):
        length = len(enz.recognition)
        patterns = (
            ("+", _pattern_sets(enz.recognition)),
            ("-", _pattern_sets(_reverse_complement_iupac(enz.recognition))),
        )
        for strand, sets in patterns:
            for start in _scan(top, sets, circular):
                key = (enzyme_index, start)
                grouped.setdefault(key, {})[strand] = None
                metadata[key] = (length, enzyme_index, enz.name)

    sites: list[RecognitionSite] = []
    for enzyme_index, start in sorted(
        grouped,
        key=lambda key: (key[1], key[0]),
    ):
        length, site_enzyme_index, enzyme_name = metadata[(enzyme_index, start)]
        sources = tuple(strand for strand in ("+", "-") if strand in grouped[(enzyme_index, start)])
        site_id = f"RS-{fingerprint}-E{enzyme_index + 1:03d}-S{start:04d}"
        sites.append(RecognitionSite(
            site_id=site_id,
            enzyme_index=site_enzyme_index,
            enzyme=enzyme_name,
            start=start,
            length=length,
            matched_sequence=_arc(top, start, length, circular),
            sources=sources,
        ))
    return tuple(sites)


def recognize_sites(
    sequence: str,
    enzymes: Union[Enzyme, Iterable[Enzyme]],
    topology: str = "linear",
) -> tuple[RecognitionSite, ...]:
    """只执行识别阶段，返回带稳定 ID 的物理识别位点。

    该函数不登记切口，也不判断交错切口是否可构造片段；调用者可把其中的
    ``site_id`` 子集传给 :func:`digest` 做确定性的部分酶切。
    """
    top = _validate_sequence(sequence)
    if isinstance(topology, bool) or not isinstance(topology, str) or topology not in TOPOLOGIES:
        raise DigestError(
            f"topology 必须是 {TOPOLOGIES} 之一，收到 {topology!r}"
        )
    enzyme_list = _normalize_enzymes(enzymes)
    return _recognize(top, topology, enzyme_list)


def _normalize_enabled_sites(
    enabled_site_ids: object,
    sites: tuple[RecognitionSite, ...],
) -> Union[frozenset[str], None]:
    if enabled_site_ids is None:
        return None
    if isinstance(enabled_site_ids, (str, bytes)) or not isinstance(enabled_site_ids, Iterable):
        raise DigestError("enabled_site_ids 必须是站点 ID 的可迭代对象，或 None")
    enabled: set[str] = set()
    for site_id in enabled_site_ids:
        if not isinstance(site_id, str):
            raise UnknownSiteError(
                f"启用位点 ID 必须是 str，收到 {type(site_id).__name__}"
            )
        enabled.add(site_id)
    valid = {site.site_id for site in sites}
    unknown = sorted(enabled - valid)
    if unknown:
        shown = ", ".join(repr(site_id) for site_id in unknown[:5])
        suffix = " ..." if len(unknown) > 5 else ""
        raise UnknownSiteError(
            f"{len(unknown)} 个启用位点 ID 不属于本次输入 DNA、拓扑与酶表的识别结果：{shown}{suffix}"
        )
    return frozenset(enabled)


def digest(
    sequence: str,
    enzymes: Union[Enzyme, Iterable[Enzyme]],
    topology: str = "linear",
    enabled_site_ids: Union[Iterable[str], None] = None,
) -> DigestResult:
    """对线性/环状双链 DNA 做完全酶切或显式位点子集的片段归属。

    纯函数：不修改任何输入，也没有内部状态，多酶批次要么整体返回结果，
    要么在任何片段构造前抛异常，不存在部分变更。

    ``enabled_site_ids=None`` 启用全部识别位点；传入空集合为确定性不切；
    传入其他可迭代对象时，仅其中属于本次识别结果的物理位点会登记切口。
    只要有一个旧/外来 ID，本函数在片段构造前抛出 :class:`UnknownSiteError`。

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

    enzyme_list = _normalize_enzymes(enzymes)

    # ---- 1. 先穷举识别，再只对启用的物理位点登记切口 ----
    recognition_sites = _recognize(top, topology, enzyme_list)
    selected_ids = _normalize_enabled_sites(enabled_site_ids, recognition_sites)
    if selected_ids is None:
        selected_sites = recognition_sites
    else:
        selected_sites = tuple(
            site for site in recognition_sites if site.site_id in selected_ids
        )
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

    for site in selected_sites:
        enzyme_index = site.enzyme_index
        enz = enzyme_list[enzyme_index]
        length = site.length
        a = enz.top_cut
        c = enz.bottom_cut
        delta = c - a
        abstract_palindrome = (
            _reverse_complement_iupac(enz.recognition) == enz.recognition
        )
        # 相同 (上链,下链) 切口只登记一次，但保留正/反向来源。
        cut_sources: dict[tuple[int, int], list[str]] = {}
        cut_strands = site.sources
        if abstract_palindrome:
            # 模式级回文时，正/反向是同一识别方向证据；规范切口按正向偏移。
            # 例如简并回文 GANTC 不能因具体窗 GATTC 再产生第二个平切口。
            cut_strands = ("+",)
        for strand in cut_strands:
            t_raw = site.start + a
            b_raw = site.start + c
            if strand == "-":
                t_raw = site.start + length - c
                b_raw = site.start + length - a
            reg = register(
                t_raw,
                b_raw,
                delta,
                f"酶 {enz.name!r} {strand} 向位点 start={site.start}",
            )
            if reg is None:
                continue
            cut_key = reg
            cut_sources.setdefault(cut_key, []).append(strand)
            if abstract_palindrome:
                cut_sources[cut_key] = list(site.sources)
        for (t, b), sources in cut_sources.items():
            raw_sites.append(CutSite(
                enzyme=site.enzyme,
                start=site.start,
                strand=sources[0],
                top_boundary=t,
                bottom_boundary=b,
                site_id=site.site_id,
                sources=tuple(sources),
            ))

    raw_sites.sort(key=lambda x: (x.start, x.strand, x.enzyme, x.site_id))
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
        return DigestResult(topology, n, cut_sites, (frag,), recognition_sites)

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

    return DigestResult(topology, n, cut_sites, tuple(fragments), recognition_sites)
