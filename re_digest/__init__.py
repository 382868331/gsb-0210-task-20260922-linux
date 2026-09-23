"""re_digest —— 限制性内切酶线性/环状双链 DNA 片段归属（标准库）。

公开接口：

* :class:`Enzyme`            酶定义（IUPAC 识别串 + 上/下链切偏移）
* :func:`digest`             单酶消化，返回 :class:`DigestResult`
* 数据结构：:class:`DigestResult` / :class:`Fragment` / :class:`EndInfo`
            / :class:`CutEvent` / :class:`SiteMatch`
* 错误：:class:`DigestError` / :class:`InvalidSequenceError`
        / :class:`InvalidEnzymeError` / :class:`UnsupportedDigestError`
* 常量：:data:`MAX_SEQUENCE_LEN` / :data:`MAX_RECOGNITION_LEN`
        / :data:`IUPAC_BASES`
* 工具：:func:`reverse_complement_dna` / :func:`reverse_complement_pattern`

坐标约定见 ``re_digest._core`` 模块文档串。
"""

from __future__ import annotations

from ._core import (
    MAX_OFFSET_ABS,
    MAX_RECOGNITION_LEN,
    MAX_SEQUENCE_LEN,
    CutEvent,
    DigestError,
    DigestResult,
    Enzyme,
    EndInfo,
    Fragment,
    InvalidEnzymeError,
    InvalidSequenceError,
    IUPAC_BASES,
    IUPAC_COMPLEMENT,
    SiteMatch,
    UnsupportedDigestError,
    digest,
    reverse_complement_dna,
    reverse_complement_pattern,
)

__all__ = [
    "MAX_OFFSET_ABS",
    "MAX_RECOGNITION_LEN",
    "MAX_SEQUENCE_LEN",
    "CutEvent",
    "DigestError",
    "DigestResult",
    "Enzyme",
    "EndInfo",
    "Fragment",
    "InvalidEnzymeError",
    "InvalidSequenceError",
    "IUPAC_BASES",
    "IUPAC_COMPLEMENT",
    "SiteMatch",
    "UnsupportedDigestError",
    "digest",
    "reverse_complement_dna",
    "reverse_complement_pattern",
]
