# re_digest —— 限制性内切酶线性/环状片段归属

离线、零依赖的可复用 Python 库：给定线性或环状双链 DNA 与限制性内切酶定义
（IUPAC 识别串 + 上/下链切偏移），穷举搜索正向与反向互补位点，转换双链接点
坐标（环状取模、重复物理切口去重），返回**按上链切口分段的双链片段**、
每条链的边界、各末端 5'/3' 突出端序列及可连接（可拼接）端型。

## 运行环境

- Python **3.14.7**（仅用标准库 `dataclasses`、`math`、`typing`）
- Linux 原生运行；装好 Python 后**无需安装任何其他依赖**，无需网络/账号/
  密钥/数据库/中间件/Docker。
- 全部计算离线、确定性（无随机源影响库行为；测试中的随机用例使用固定种子）。

## 准确命令

```bash
# 必要测试（44 个用例）
python -m unittest discover -s tests -v

# 固定输入演示（毫秒级，远小于 8 秒；真实计算 + 真实触发的拒绝边界）
python demo.py
```

## 快速上手

```python
from re_digest import Enzyme, digest

# EcoRI: G^AATTC / CTTAA^G —— 上链偏移 1，下链偏移 5（都相对匹配起点）
eco = Enzyme("EcoRI-like", "GAATTC", top_offset=1, bottom_offset=5)
r = digest("AAGAATTCAA", eco)                 # 线性
rc = digest("AAGAATTCTT", eco, circular=True) # 环状

ev = r.cut_events[0]
ev.top_cut, ev.bottom_cut      # (3, 7)
f0, f1 = r.fragments
f0.top_sequence                # "AAG"
f0.bottom_sequence_5to3        # "AATTCTT"（下链按 5'->3' 书写）
f0.right_end.kind              # "5prime"
f0.right_end.overhang          # "AATT"
f0.right_end.compatible_with(f1.left_end)  # True —— 可拼接
```

## 公开接口

| 名称 | 说明 |
|---|---|
| `Enzyme(name, recognition, top_offset, bottom_offset)` | 不可变酶定义；`enzyme.stagger = bottom_offset - top_offset`（>0：5'突出，<0：3'突出，0：平端） |
| `digest(sequence, enzyme, *, circular=False) -> DigestResult` | 单酶消化主入口 |
| `DigestResult` | `sequence/circular/length/enzyme_name`、`site_matches`、`cut_events`、`fragments`、`positions_scanned` |
| `SiteMatch(start, orientation)` | 识别命中；`orientation` 为 `"forward"`/`"reverse"`，回文位点正反向各一条 |
| `CutEvent(top_cut, bottom_cut, match_starts)` | 去重后的物理双切事件；`match_starts` 记录全部贡献命中 |
| `Fragment` | `index/shape`、`top_start/top_end/bottom_start/bottom_end`、`top_sequence`、`bottom_sequence_5to3`、`left_end/right_end` |
| `EndInfo` | `kind`（`blunt`/`5prime`/`3prime`）、`overhang`（突出链 5'->3'，平端为 `""`）、`top_boundary/bottom_boundary`、`ligation_key`、`compatible_with(other)` |
| 错误 | `DigestError`（基类，也是 `ValueError`）、`InvalidSequenceError(index=)`、`InvalidEnzymeError(field=)`、`UnsupportedDigestError` |
| 常量 | `MAX_SEQUENCE_LEN=5000`、`MAX_RECOGNITION_LEN=12`、`IUPAC_BASES` |
| 工具 | `reverse_complement_dna`、`reverse_complement_pattern` |

所有数据结构均为 `frozen=True` 的不可变值对象；`digest` 无状态，重复调用
结果相等，失败时不产生任何部分变更（纯函数，出错只抛异常）。

## 坐标与切偏移约定（契约）

1. 碱基沿上链 5'->3' 编号 `0..n-1`；切口 `c` 表示“碱基 c 之前”的键。
2. 上/下链切偏移**都相对匹配起点 s**（模式第一碱基位置）：
   上链切口 `s + top_offset`，下链切口 `s + bottom_offset`。
3. 反向互补命中（模式在 j..j+L-1 的下链意义上命中）：
   上链切口 `j + L - bottom_offset`，下链切口 `j + L - top_offset`。
4. 线性 DNA：两个切口都必须严格位于 `(0, n)` 内该事件才切断；任一切点
   落在末端键（0/n）或分子之外 = “切到空气”，整事件不产生切口。
5. 环状 DNA：识别允许跨零，切口坐标对 n 取模，键 0 与 n 等同；同一
   `(top_cut, bottom_cut)` 的物理事件去重（回文酶正/反向命中即合并，
   贡献来源保留在 `match_starts`）。
6. 当相邻上链切口的环上间距 `<= |stagger|` 时，切口互相套叠，无法归属为
   彼此独立的常规双链片段，**明确抛出 `UnsupportedDigestError`**，不做
   任意配对。环状单切且 `n <= |stagger|`（取模吞掉交错）同样拒绝。
7. 片段按上链切口有序划分；`Fragment` 同时给出下链覆盖弧与两端边界。
   环未切时返回唯一一个 `shape="circular"`、无末端的片段；切开后所有产物
   都是 `shape="linear"`。

`EndInfo.overhang` 始终按**突出链自身 5'->3'** 书写；同一物理切口两侧的
两个突出端序列互为反向互补，`compatible_with`/`ligation_key` 据此判定
（平端全部互连；5' 与 3' 端不互通；非互补黏端不互通）。

## 输入上限与错误语义

- 序列长度 `1..5000`（含），仅接受大写 `A/C/G/T`；空串、超长、含小写或
  非法字符均抛 `InvalidSequenceError`，非法字符带 0 基 `index`。
- 识别串长度 `1..12`，仅接受 IUPAC 码（`ACGTRYSWKMBDHVN`）；非法带出错
  位置，酶字段错误带 `field` 名。切偏移支持负数/越位（类型 IIS 酶），
  绝对值护栏为 5000。
- 类型契约：`top_offset/bottom_offset` 必须是真正的 `int`——拒绝 `bool`
  冒充、拒绝一切 `float`（含 NaN / ±Infinity）；`circular` 必须是真正的
  `bool`；`sequence` 必须是 `str`。类型不符抛 `TypeError`。
- 非法输入一律抛出**可定位**异常，不静默纠正、不猜测补全。

## “无解”与“预算耗尽”

搜索是对全部合法起点的**穷举**（线性 `n-L+1` 个、环状 `n` 个，
`DigestResult.positions_scanned` 如实记录）。因此无命中是**证明无解**，
返回未切分子（线性：一个平端完整片段；环状：一个无末端环状片段）。
本算法不存在搜索预算这一状态，也就不存在“预算耗尽”，两者不可能混淆。

## 数值/浮点口径

本库不产生任何浮点结果：匹配、坐标、取模、长度比较全部是整数/字符串精确
运算，序列比较逐字符精确相等。因此**没有误差带、没有容差**，边界条件
（间距恰等于交错长度、切点恰在末端键、长度恰为 5000 等）按精确值判定，
不以宽容差掩盖。

## 设计取舍

- **单酶、单分子、经典双链模型**：不实现多酶混合消化、甲基化阻断、切刻酶、
  切口平移等；有限子集只按本 README 契约工作，不宣称行业标准（如
  REBASE/SBOL）完整兼容。
- 线性末端越位事件按物理事实整体丢弃（“切空气”），而不是把切口夹到边界。
- 简并 IUPAC 按集合匹配；自反向互补的奇数长度码（如 GANTC 偏移居中）正、
  反方向切在中央简并碱基两侧的不同键上，是两个不同平端事件，不去重。
- 环状且 `L > n` 时识别窗不可能落在 L 个不同碱基上，直接判定无位点
  （不把同一段碱基绕环读两遍）。
- 复杂度 O(n·L) 时间、O(n) 输出规模；n=5000、L=12 时为亚毫秒~毫秒级。

## 已知限制

- 不支持交错切口套叠拓扑（抛 `UnsupportedDigestError`，见约定 6）。
- 不做连接反应模拟，只提供末端兼容性判定（`compatible_with`）。
- 不提供 CLI 产品/文件读写；`demo.py` 仅是固定输入的打印演示。
- 仓库布局：`re_digest/`（核心库）、`tests/`（测试与独立预言机）、
  `demo.py`、`README.md`。测试中的期望值来自手工常量与 `tests/oracle.py`
  的独立暴力实现（**不导入被测核心来生成期望值**）；重组不变量仅作补充。
