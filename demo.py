"""固定输入演示：真实计算限制性酶切片段归属，并打印一个正常结果和一个实际触发的拒绝边界。

运行：python demo.py
数据全部在本文件内确定性给出，不读网络/文件/环境变量。
"""

from restriction_digest import (
    CandidateError,
    CutMap,
    DigestUnsupportedError,
    Enzyme,
    compatible,
    digest,
    map_from_digest,
    match_candidates,
)


def describe_end(end):
    if end is None:
        return "无末端（环状片段）"
    label = {"blunt": "平端", "5_prime": "5' 黏端", "3_prime": "3' 黏端"}[end.kind]
    if end.kind == "blunt":
        return f"{label}（上链边界 {end.top_boundary} / 下链边界 {end.bottom_boundary}）"
    return (f"{label}，突出 {end.overhang!r}（5'→3'，{end.overhang_length} nt）"
            f"，边界 上链 {end.top_boundary} / 下链 {end.bottom_boundary}")


def main():
    print("=" * 70)
    print("正常结果：EcoRI（G^AATTC，5' 黏端）单切环状 DNA")
    print("=" * 70)
    ecori = Enzyme("EcoRI", "GAATTC", 1, 5)
    seq = "GGGAATTCCCAA"  # 12 bp 环状，唯一位点起点 3
    res = digest(seq, ecori, topology="circular")
    print(f"底物（环）: {seq}  长度 {res.length}")
    site = res.cut_sites[0]
    print(f"命中: 酶 {site.enzyme}，起点 {site.start}，链 {site.strand}；"
          f"上链切口 {site.top_boundary}，下链切口 {site.bottom_boundary}（已取模）")
    print(f"物理断裂数: {res.cut_count}，片段数: {len(res.fragments)}")
    frag = res.fragments[0]
    print(f"片段0（线性化，从上链坐标 {frag.top_start} 起绕环）:")
    print(f"  上链 5'→3': {frag.top_strand}  ({frag.top_length} nt)")
    print(f"  下链 5'→3': {frag.bottom_strand}  ({frag.bottom_length} nt)")
    print(f"  配对区长度: {frag.paired_length} bp")
    print(f"  左端: {describe_end(frag.left_end)}")
    print(f"  右端: {describe_end(frag.right_end)}")
    print()
    can_recircularize = compatible(frag.right_end, frag.left_end)
    print(f"可拼接端型检查：右端 × 左端 -> {'可互补拼接（可重新环化）' if can_recircularize else '不可拼接'}")

    # 再附一个平端线性小例，直观对照端型。
    smai = Enzyme("SmaI", "CCCGGG", 3, 3)
    r2 = digest("AAACCCGGGTT", smai)
    print()
    print("对照（SmaI 平端线性酶切）:")
    for f in r2.fragments:
        print(f"  片段{f.index}: 上链 {f.top_strand} / 下链 {f.bottom_strand}；"
              f"左 {f.left_end.kind}，右 {f.right_end.kind}")
    print(f"  相邻端可拼接: {compatible(r2.fragments[0].right_end, r2.fragments[1].left_end)}")

    print()
    print("=" * 70)
    print("新能力：候选环状切位图 × 完全消化片段长度多重集（环长 12）")
    print("=" * 70)
    # 完全消化后测得四段，长度多重集 {2,2,4,4}（测量顺序无意义）。
    observed = [4, 2, 4, 2]
    # 候选池（切位边界显式给出；图A' 是图A 整体旋转 +1，图B 排列不同）：
    candidates = [
        CutMap(12, (0, 2, 4, 8), label="图A (0,2,4,8) 间隙2,2,4,4"),
        CutMap(12, (1, 3, 5, 9), label="图A旋转+1 (1,3,5,9)"),
        CutMap(12, (0, 2, 6, 8), label="图B (0,2,6,8) 间隙2,4,2,4"),
        CutMap(12, (0, 3, 6, 9), label="图C (0,3,6,9) 间隙3,3,3,3"),
        CutMap(12, (), label="未切环 () 间隙12"),
    ]
    print(f"测得片段长度: {observed}（多重集，顺序无关）")
    for c in candidates:
        print(f"  候选 {c.label}: 完全消化长度 {c.fragment_lengths()}")
    result = match_candidates(observed, candidates)
    print(f"检查候选 {result.examined} 个，旋转等价去重后 {result.unique_classes} 类")
    print(f"匹配 {len(result.matches)} 个：")
    for m in result.matches:
        print(f"  ✓ {m.label}")
    print("  （图A 与其旋转图只保留先出现者；图B 多重集相同但排列不同，")
    print("    不与图A 旋转等价，故作为第二个独立答案保留；图C/未切环不匹配）")

    # 与旧切割规则打通：真实酶切一个环，提出切位图，再用测得长度匹配。
    seq = "AAACCCGGGTTT"  # 12 bp 环，SmaI 唯一位点起点 3 -> 平端切口 6
    dr = digest(seq, smai, topology="circular")
    derived = map_from_digest(dr)
    print()
    print(f"旧规则联动：digest({seq!r}, SmaI, 环状) -> 切位图 "
          f"length={derived.length}, cuts={derived.cuts}")
    one_cut_pool = [CutMap(12, (), label="未切环"),
                    CutMap(12, (6,), label="单切@6"),
                    CutMap(12, (0,), label="单切@0（与@6旋转等价）")]
    r1 = match_candidates([derived.length], one_cut_pool)
    print(f"  对测得长度 [{derived.length}] 匹配：{[m.label for m in r1.matches]}")
    print("  -> 未切环与单切环都只产生一个全长片段，且互不等价，均列出；")
    print("     两个不同位置的单切口互为旋转，只保留一个。")

    # 边界：无匹配是穷举后被证明的结果。
    r_none = match_candidates([5, 7], candidates)
    print(f"  对找不到的长度 [5,7]：matched={r_none.matched}，matches={r_none.matches}")

    # 实际错误：非法候选（重复切口 / 越界）必须明确拒绝，不静默纠正。
    print()
    print("-- 新功能的实际错误：非法切位图候选 --")
    try:
        CutMap(12, (3, 3))
    except CandidateError as exc:
        print(f"CutMap(12, (3, 3)) -> CandidateError：{exc}")
    try:
        match_candidates([2, 0, 10], [CutMap(12, (0, 2))])
    except CandidateError as exc:
        print(f"observed 含 0 长片段 -> CandidateError：{exc}")
    else:
        raise AssertionError("应当触发 CandidateError 但未触发")

    print()
    print("=" * 70)
    print("拒绝边界：同一位点两个酶的交错切口重叠，归属非简单片段")
    print("=" * 70)
    # 同一识别串 GAATTC：E1 切 (1,5)，E2 切 (2,5)。
    # 两个酶产生的断裂争夺同一个物理下链切口 5（分别配上链切口 1 和 2），
    # 重复物理切口无法去重为同一末端 -> 非简单归属，明确拒绝。
    e1 = Enzyme("EcoRI", "GAATTC", 1, 5)
    e2 = Enzyme("EcoRI-shift", "GAATTC", 2, 5)
    print(f"底物（环）: GAATTCAAAA；酶: {e1.name}({e1.top_cut},{e1.bottom_cut})，"
          f"{e2.name}({e2.top_cut},{e2.bottom_cut})")
    try:
        digest("GAATTCAAAA", [e1, e2], topology="circular")
    except DigestUnsupportedError as exc:
        print(f"如契约预期抛出 DigestUnsupportedError：\n  {exc}")
        print("-> 模块拒绝猜测配对，没有返回乱配的片段。")
    else:
        raise AssertionError("应当触发 DigestUnsupportedError 但未触发")


if __name__ == "__main__":
    main()
