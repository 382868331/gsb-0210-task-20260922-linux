"""环状 DNA 一次部分酶切的期望片段长度分布演示。

运行：python demo_iteration.py
全部输入由代码即时构造与计算，不访问网络、文件、账号或外部服务。
"""

from restriction_digest import (
    Enzyme,
    partial_digest_expectation,
)


def main():
    print("=" * 72)
    print("例 1：两个酶识别到完全重合的切点，按‘至少一个切开’合并")
    print("=" * 72)
    eco = Enzyme("EcoRI", "GAATTC", 1, 5)
    clone = Enzyme("EcoRI-clone", "GAATTC", 1, 5)
    seq12 = "GGGAATTCCCAA"
    small = partial_digest_expectation(
        seq12,
        [eco, clone],
        (0.4, 0.5),
    )
    merged = small.merged_cut_sites[0]
    print(f"环长: {small.length}；识别证据数: {len(merged.evidence)}")
    print(
        "合并后切点: 上链边界 "
        f"{merged.top_boundary}，下链边界 {merged.bottom_boundary}，"
        f"P(至少一个切开) = {merged.probability:.6f}"
    )
    for fragment in small.expected_fragments:
        label = "未切整环" if fragment.kind == "circular" else "全长线性片段"
        print(
            f"  长度 {fragment.length:2d}，{label}："
            f"期望 {fragment.expected_count:.6f} 条"
        )
    print(f"期望总片段长度: {small.expected_total_length:.12f}")

    print()
    print("=" * 72)
    print("例 2：1200 个不同切点的部分酶切；O(m²) 公式，不枚举 2^1200 子集")
    print("=" * 72)
    one_base = Enzyme("OneBaseA", "A", 1, 1)
    size = 1200
    large_sequence = "A" * size

    def probability_for(site):
        return 0.01 + 0.0001 * (site.start % 31)

    large = partial_digest_expectation(
        large_sequence,
        one_base,
        probability_for,
    )
    nonzero_lengths = len(large.expected_fragments)
    print(f"环长/合并后的实际切点数: {large.length} / {len(large.merged_cut_sites)}")
    print(f"非零期望长度类型数: {nonzero_lengths}")
    print("前 8 个长度分布（长度, 类型, 期望条数）:")
    for fragment in large.expected_fragments[:8]:
        print(
            f"  ({fragment.length:4d}, {fragment.kind:8s}, "
            f"{fragment.expected_count:.12f})"
        )
    print(f"期望总片段长度: {large.expected_total_length:.12f}")
    assert abs(large.expected_total_length - size) < 1e-8


if __name__ == "__main__":
    main()
