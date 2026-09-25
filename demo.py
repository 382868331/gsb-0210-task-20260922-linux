"""固定输入演示：真实计算限制性酶切片段归属，并打印一个正常结果和一个实际触发的拒绝边界。

运行：python demo.py
数据全部在本文件内确定性给出，不读网络/文件/环境变量。
"""

from restriction_digest import (
    DigestUnsupportedError,
    Enzyme,
    UnknownSiteError,
    compatible,
    digest,
    recognize_sites,
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
    print("部分酶切：两个酶识别同一物理位点、上链切口相同而下链切口不同")
    print("=" * 70)
    partial_seq = "GAATTCAAAA"
    cut_a = Enzyme("SharedTop-A", "GAATTC", 1, 5)
    cut_b = Enzyme("SharedTop-B", "GAATTC", 1, 6)
    partial_sites = recognize_sites(partial_seq, [cut_a, cut_b])
    print(f"底物（线性）: {partial_seq}；识别到 {len(partial_sites)} 个可独立选择的位点")
    for site in partial_sites:
        print(f"  ID={site.site_id}；酶={site.enzyme}；来源={site.sources}")
    selected = digest(
        partial_seq,
        [cut_a, cut_b],
        enabled_site_ids={partial_sites[0].site_id},
    )
    cut = selected.cut_sites[0]
    print(f"只启用 {cut.enzyme}: 上链切口 {cut.top_boundary}，下链切口 {cut.bottom_boundary}")
    for frag in selected.fragments:
        print(f"  片段{frag.index}: 上链 {frag.top_strand} / 下链 {frag.bottom_strand}")
    empty = digest(partial_seq, cut_a, enabled_site_ids=frozenset())
    print(f"空集合: 识别位点 {len(empty.recognition_sites)} 个，实际切口 {len(empty.cut_sites)} 个，"
          f"片段为完整分子 {empty.fragments[0].top_strand!r}")
    try:
        digest(partial_seq, cut_b, enabled_site_ids={partial_sites[0].site_id})
    except UnknownSiteError as exc:
        print(f"复用旧酶表 ID 时如契约预期拒绝：{exc}")
    else:
        raise AssertionError("旧酶表 ID 应当被拒绝但未触发")

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
