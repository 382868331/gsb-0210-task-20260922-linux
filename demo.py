#!/usr/bin/env python3
"""固定输入的离线演示：真实调用 re_digest 计算并打印。

内容：
  1. 正常结果：一条线性 dsDNA 被 5' 黏端酶切两次（含回文位点正/反向命中去重），
     打印双链片段、上下链边界、5'/3' 端与突出端，并按 ligation_key 列出
     “可拼接端型”分组（哪些末端彼此兼容可连接）。
  2. 边界拒绝：交错切口互相套叠的输入真实触发 UnsupportedDigestError；
     另有一个非法序列输入真实触发带位置的 InvalidSequenceError。

不联网、不 sleep、不含预录结果；整段运行毫秒级。
"""

from __future__ import annotations

from re_digest import (
    DigestResult,
    Enzyme,
    Fragment,
    InvalidSequenceError,
    UnsupportedDigestError,
    digest,
)

KIND_CN = {"blunt": "平端", "5prime": "5'突出", "3prime": "3'突出"}


def show_result(title: str, result: DigestResult) -> None:
    bar = "=" * 72
    print(bar)
    print(title)
    print(
        f"序列长度={result.length} 拓扑={'环状' if result.circular else '线性'} "
        f"酶={result.enzyme_name} 扫描起点数={result.positions_scanned}"
    )
    print(f"识别命中 {len(result.site_matches)} 条；去重后物理双切事件 "
          f"{len(result.cut_events)} 个：")
    for ev in result.cut_events:
        src = ", ".join(f"{d}:{s}" for d, s in ev.match_starts)
        print(f"  上链切口={ev.top_cut:>3}  下链切口={ev.bottom_cut:>3}  <- [{src}]")

    print(f"双链片段 {len(result.fragments)} 个：")
    for f in result.fragments:
        show_fragment(f)

    show_ligation_groups(result)


def show_fragment(f: Fragment) -> None:
    if f.shape == "circular":
        print(f"  [片段 {f.index}] 未切环状分子（无末端）")
        print(f"      上链 5'->3': {f.top_sequence}")
        print(f"      下链 5'->3': {f.bottom_sequence_5to3}")
        return
    le, re_ = f.left_end, f.right_end
    print(
        f"  [片段 {f.index}] 上链[{f.top_start}..{f.top_end}) "
        f"下链[{f.bottom_start}..{f.bottom_end})"
    )
    print(f"      上链 5'->3': {f.top_sequence}")
    print(f"      下链 5'->3': {f.bottom_sequence_5to3}")
    print(
        f"      左端: {KIND_CN[le.kind]:>6} 突出={le.overhang or '-'} | "
        f"右端: {KIND_CN[re_.kind]:>6} 突出={re_.overhang or '-'}"
    )


def show_ligation_groups(result: DigestResult) -> None:
    """按兼容键对所有末端分组：同组内任意两端可连接（平端全互通）。"""
    groups: dict[tuple[object, ...], list[str]] = {}
    for f in result.fragments:
        if f.shape == "circular":
            continue
        for side, end in (("左", f.left_end), ("右", f.right_end)):
            groups.setdefault(end.ligation_key, []).append(f"片段{f.index}{side}端")
    print("可拼接端型分组（同组两端可连接）：")
    for i, (key, members) in enumerate(sorted(groups.items(), key=lambda kv: kv[0])):
        kind = key[0]
        detail = "平端互连" if kind == "blunt" else f"{KIND_CN[kind]} 互补对"
        print(f"  组{i}: {detail:>10} <- {', '.join(members)}")


def main() -> None:
    # ---- 正常结果：线性，EcoRI 样 5' 黏端，两个位点 -------------------------
    eco = Enzyme("EcoRI-like G^AATTC/CTTAA^G", "GAATTC", 1, 5)
    seq = "AAGAATTCAAAGAATTCGG"   # 位点 s=2 与 s=11；回文位点正反双向命中
    result = digest(seq, eco, circular=False)
    show_result("正常结果：线性 DNA 双位点 5' 黏端消化", result)

    # 演示一个具体的连接判断（实时计算，非预录）
    f0r, f1l = result.fragments[0].right_end, result.fragments[1].left_end
    print(
        f"连接检查: 片段0右端[{f0r.overhang}] 与 片段1左端[{f1l.overhang}] "
        f"可连接 = {f0r.compatible_with(f1l)}"
    )

    # ---- 边界 1：交错切口套叠 -> 明确 unsupported ----------------------------
    print("=" * 72)
    enz = Enzyme("dense-G (0,2)", "G", 0, 2)  # 5' 突出 2nt，G 处处是位点
    try:
        digest("AGGGA", enz)
    except UnsupportedDigestError as exc:
        print(f"实际触发拒绝（非静默纠正、非任意配对）: {type(exc).__name__}: {exc}")

    # ---- 边界 2：非法序列 -> 带位置的可定位错误 ------------------------------
    try:
        digest("ACGTxACGT", eco)
    except InvalidSequenceError as exc:
        print(f"实际触发输入拒绝: {type(exc).__name__}: {exc}（定位 index={exc.index}）")


if __name__ == "__main__":
    main()
