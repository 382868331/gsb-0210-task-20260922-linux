# 限制性内切酶的环状片段归属（restriction_digest）

离线、纯标准库的限制性酶切片段归属库：给定线性或环状双链 DNA 与一个/多个
限制酶定义，搜索正向与反向互补识别位点，转换并取模切口坐标，对重复物理切口
去重，返回**按上链切口分段**的双链片段、上/下链边界以及 5′/3′ 突出端序列；
当交错切口造成无法做无交叉简单归属的片段时，明确抛出 `DigestUnsupportedError`，
绝不猜测配对。

在此之上提供**候选环状切位图筛查**：给定一组候选环状序列、显式酶定义和
完全消化的片段长度多重集，调用同一切割规则逐一评估，按旋转等价类去重后
列出长度图样精确匹配的候选（`match_candidates`）。

## 运行环境

- Python **3.14.7**（仅使用标准库 `dataclasses`、`math`、`typing`）。
- Linux 原生运行；装好 Python 后**无需安装任何第三方依赖**，无需网络、
  账号、密钥、数据库、中间件或 Docker。

## 运行命令（在仓库根目录执行）

```bash
python -m unittest discover -s tests -v   # 必要测试（77 个用例）
python demo.py                            # 固定输入演示（约 0.1 秒）
```

## 公开接口

```python
from restriction_digest import Enzyme, digest, compatible, match_candidates

ecori = Enzyme(name="EcoRI", recognition="GAATTC", top_cut=1, bottom_cut=5)
result = digest("GGGAATTCCC", ecori, topology="linear")   # 或 "circular"

# 候选环状切位图筛查（旋转等价去重，片段长度精确多重集匹配）
screen = match_candidates(
    ["GAATTCGAATTC", "TTCGAATTCGAA"], ecori, [6, 6],
)
screen.matches[0].indices          # (0, 1)：两条是同一环的不同旋转写法
```

### `Enzyme`

| 字段 | 含义 |
|---|---|
| `name` | 非空字符串 |
| `recognition` | 大写 IUPAC 识别串，长度 1..12：`ACGT R Y S W K M B D H V N` |
| `top_cut` (`a`) | 上链切偏移（边界编号，相对识别串起点） |
| `bottom_cut` (`b`) | 下链切偏移（**同样以上链坐标表示**，相对识别串起点） |

边界 `k` 位于第 `k-1` 与第 `k` 个碱基之间：`k=0` 在识别串之前，`k=L` 在其后。
- 正向命中起点 `s`：上链切口 `s+a`，下链切口 `s+b`。
- 反向互补命中起点 `s`（反向互补串在上链的起点）：
  上链切口 `s+L−b`，下链切口 `s+L−a`（酶整体旋转 180°）。
- 回文识别串（如 `GAATTC`）的正/反向命中是同一次物理断裂，只计一次。
- `stagger = b−a`：`>0` 为 5′ 突出，`<0` 为 3′ 突出，`=0` 为平端。

示例：EcoRI `G^AATTC` / `CTTAA^G` → `(1, 5)`；SmaI `CCC^GGG` → `(3, 3)`；
3′ 突出的 `ATGCA^T` / `T^ACGTA` → `(5, 1)`。

### `digest(sequence, enzymes, topology="linear") -> DigestResult`

- `sequence`：大写 `ACGT` 字符串，长度 1..**5000**。底物不接受简并碱基、
  小写或空白（不做静默纠正）。
- `enzymes`：单个 `Enzyme` 或非空可迭代对象（多酶同时完全酶切）。
- `topology`：`"linear"` 或 `"circular"`。
- 纯函数、无内部状态；多酶批次在任何片段构造前完成全部校验与登记，
  失败即整体抛出异常，不存在部分变更。

`DigestResult`：

- `topology`, `length`；
- `cut_sites: tuple[CutSite, ...]`：每条识别命中的证据（酶名、起点、
  `+`/`-` 链、取模后的上/下链切口）；`cut_count` 为按
  `(上链边界, 下链边界)` 去重后的物理断裂数；
- `fragments: tuple[Fragment, ...]`：按上链切口排序的片段（环状按循环序，
  可跨越坐标 0；线性顺序的两端为分子天然末端）。

`Fragment`：`index`、`topology`、`circular`（仅未切环状分子为 `True`）、
`top_start/top_length`（上链起点坐标与跨度，跨零无歧义）、
`bottom_start/bottom_length`、`paired_length`（双链配对碱基数）、
`top_strand`/`bottom_strand`（**均按各自 5′→3′** 报告）、
`left_end`/`right_end: EndInfo | None`。

`EndInfo`：`kind ∈ {"blunt","5_prime","3_prime"}`、
`top_boundary`/`bottom_boundary`、`overhang`（突出序列，**一律按突出链
5′→3′**，平端为 `""`）、`overhang_length`。

未切环状分子返回单个 `circular=True`、两端为 `None` 的片段；未切线性分子
返回单个两端平端的完整片段。

### `compatible(left_end, right_end) -> bool`

两个末端能否互补拼接：平端配平端；黏端要求同类且一条突出的 5′→3′ 序列
是另一条的反向互补（如 EcoRI 的 `AATT` 自补）。`None` 末端不可拼接。

### `match_candidates(candidates, enzymes, fragment_lengths) -> ScreenResult`

候选环状切位图筛查：对每个候选环状序列调用 `digest(..., topology="circular")`
做完全消化，取各片段的**上链切口跨度**作为片段长度（酶切图谱通行约定，
总和恒等于环长），与观测的 `fragment_lengths` 做**多重集精确匹配**
（重复片段按计数匹配，如 `(6,6)` 不等于 `(6,)`）。

- `candidates`：候选环状序列的可迭代对象（大写 `ACGT`，校验规则同
  `digest` 的底物，错误信息带 `candidates[i]` 下标）；空列表合法，
  结果为零匹配。
- `enzymes`：单个 `Enzyme` 或非空可迭代对象（同 `digest`）。
- `fragment_lengths`：观测片段长度多重集，非空正整数可迭代对象；
  拒绝 `bool`、任何 `float`（`NaN`/`±Infinity` 有专门报错）与非正整数，
  抛 `FragmentLengthsError`。顺序无关，内部升序归一。
- **旋转等价去重**：同一物理环状分子的不同旋转写法是同一个候选，
  按最小字典序旋转（Booth 算法，O(n)）归一去重，每个等价类只评估一次、
  只报告一次，`representative` 取输入中该类首次出现的写法，`indices`
  记录该类全部输入下标。
- 纯函数：候选与片段长度在任何评估之前完成整体校验，非法输入立即抛出，
  不存在部分结果。

`ScreenResult`：

- `target_lengths`：升序归一后的观测多重集；
- `candidate_count` / `class_count`：输入候选总数 / 旋转去重后的类数；
- `matches: tuple[CandidateMatch, ...]`：精确匹配的候选类，按输入首次
  出现顺序。`CandidateMatch` 含 `representative`、`indices`、`length`、
  `fragment_lengths`（升序）、`cut_count`（物理断裂数，0 = 未切环）；
- `unsupported: tuple[UnsupportedCandidate, ...]`：`digest` 判定为非简单
  几何（`DigestUnsupportedError`）的候选类，附 `reason`；它们不参与匹配，
  也不中断其余候选的评估。

注意：未切环与单切环的长度多重集都是 `(N,)`，凝胶长度图样无法区分，
用 `CandidateMatch.cut_count` 辨别。

## 坐标与归属约定（设计取舍）

1. **上链统一坐标系**：所有切口都是上链 5′→3′ 的边界编号；下链切口也用
   上链坐标表达，避免双坐标制换算错误。环状坐标对 `N` 取模。
2. **物理切口去重**：`(链, 坐标)` 相同的切口只算一次；若同一物理切口被
   两个不同的对侧切口争夺（无法去重为同一末端），报 `DigestUnsupportedError`。
3. **片段 = 相邻上链切口之间的区间**，下链边界由该物理断裂的固有错位
   `δ = b−a` 决定；简单性判据是：下链跨度为正、配对区长度 `> 0`
   （两端突出互不侵入、双链不解离）且环状交错量 `< N`。任何一条不满足
   都明确报 unsupported，而不是输出交叉配对的片段。
4. **线性分子边缘**：平端切口恰在 0/N 与天然末端重合，按无切割处理；
   一条链切在末端、另一条切在内部会切出单链片段，报 unsupported；
   切口超出 `[0, N]` 同样报 unsupported。
5. **突出序列方向**：统一按突出链自身 5′→3′ 报告，因此可直接用于
   拼接判定（配对接合处两条链方向相反，比较反向互补）。

## 错误语义

- `DigestError(ValueError)`：本模块错误基类。
  - `EnzymeDefinitionError`：酶定义非法（识别串字符/长度、偏移类型等）。
  - `SequenceError`：底物非法，错误信息含**从 0 开始的位置下标**；
    候选筛查时带 `candidates[i]` 下标。
  - `DigestUnsupportedError`：输入合法，但几何上不属于支持的简单片段。
  - `FragmentLengthsError`：观测片段长度多重集非法（空、非正整数、
    bool、float/NaN/Infinity 等）。
- 整数参数（`top_cut`/`bottom_cut`）拒绝 `bool`；任何 `float` 都拒绝，
  其中 `NaN`/`±Infinity` 有专门报错信息。
- “没有识别位点”是穷举搜索后**被证明的无解**（返回未切分子），与任何
  预算耗尽无关；搜索复杂度 O(N·L)（N≤5000，L≤12），无预算/超时失败路径。
- 全部结果为整数与字符串，**没有浮点输出**，故无误差口径/容差问题。

## 已知限制（本契约的有限子集，非行业标准全兼容）

- 底物只接受具体碱基 `ACGT`；IUPAC 简并仅用于酶识别串；不支持 U、
  小写自动转换、RNA、修饰碱基或缺口/泡状底物。
- 不模拟酶切反应动力学、星号活性、甲基化阻断、温度/缓冲液条件，
  也不输出片段浓度；只做确定性的完全酶切几何归属。
- 交错切口互相侵入（配对区 ≤ 0）、同一物理切口配对冲突、
  `|错位| ≥ 环长`、线性边缘单链切割等情形一律 unsupported。
- 单条序列上限 5000 bp、识别串上限 12 bp，均为接口硬边界。
- 候选筛查只做**精确长度多重集匹配**：不做测量区间/容差反演、
  不做部分消化或混合样本建模；片段长度取上链切口跨度（黏端片段的
  上/下链跨度可能不同，本接口固定取上链口径）。
- 长度多重集无法区分未切环与单切环（都是 `(N,)`），也无法区分
  产生相同切口间距图样的不同序列；这些都如实作为不同/相同候选类报告，
  不做额外猜测。

## 测试与独立参考说明

- 测试用例的期望值全部手算，或由测试文件内独立的小函数（`ref_revcomp`、
  `ref_circular_lengths`、显式列窗等）给出，**不调用被测核心生成期望值**；
  链长守恒、沃森-克里克配对等往返/不变量检查只作为补充
  （`TestRoundTripInvariants`）。
- 覆盖：平端、5′/3′ 黏端、单切环状（平端/黏端再环化）、跨零识别、
  反向酶坐标转换、回文单次计数、重复位点、重复物理切口去重、不切、
  IUPAC（含 N）、各类 unsupported 与全部输入校验。
- 候选筛查覆盖：精确匹配、重复片段多重性、旋转等价去重（含全部 12 个
  旋转写法）、未切环/单切环同图样、无匹配、多酶图样、几何不可归属候选
  收入 unsupported 不中断、以及片段长度与候选的全部输入校验。
