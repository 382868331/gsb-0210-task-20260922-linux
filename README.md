# 限制性内切酶的环状片段归属（restriction_digest）

离线、纯标准库的限制性酶切片段归属库：给定线性或环状双链 DNA 与一个/多个
限制酶定义，搜索正向与反向互补识别位点，转换并取模切口坐标，对重复物理切口
去重，返回**按上链切口分段**的双链片段、上/下链边界以及 5′/3′ 突出端序列；
当交错切口造成无法做无交叉简单归属的片段时，明确抛出 `DigestUnsupportedError`，
绝不猜测配对。

在此之上提供第二轮新增能力——**候选环状切位图 × 完全消化片段长度多重集
匹配**：给定测得的完全消化片段长度（多重集，顺序无关、允许重复）和一组
候选环状切位图，用上面的旧切割规则（切位边界显式给出）列出匹配候选，
长度精确比较（无容差），并按环旋转等价去重。切位图既可手工构造，也可由
环状 `digest` 结果直接导出。

## 运行环境

- Python **3.14.7**（仅使用标准库 `dataclasses`、`math`、`collections`、
  `typing`）。
- Linux 原生运行；装好 Python 后**无需安装任何第三方依赖**，无需网络、
  账号、密钥、数据库、中间件或 Docker。

## 运行命令（在仓库根目录执行）

```bash
python -m unittest discover -s tests -v   # 必要测试（75 个用例）
python demo.py                            # 固定输入演示（约 0.05 秒）
```

## 公开接口

```python
from restriction_digest import (
    Enzyme, digest, compatible,
    CutMap, MatchResult, CandidateError,
    circular_fragment_lengths, map_from_digest, match_candidates,
)

ecori = Enzyme(name="EcoRI", recognition="GAATTC", top_cut=1, bottom_cut=5)
result = digest("GGGAATTCCC", ecori, topology="linear")   # 或 "circular"

# 第二轮：候选环状切位图匹配完全消化片段长度多重集
m = match_candidates(
    observed=[4, 2, 4, 2],                          # 测得长度，顺序无关
    candidates=[CutMap(12, (0, 2, 4, 8), "A"),     # 环长 12，切口显式
                CutMap(12, (0, 2, 6, 8), "B")],
)
print([c.label for c in m.matches])                 # ['A', 'B']
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

### 候选环状切位图匹配（第二轮新增）

```python
CutMap(length: int, cuts: tuple[int, ...], label: str = "")
circular_fragment_lengths(length, cuts) -> tuple[int, ...]
map_from_digest(result: DigestResult) -> CutMap
match_candidates(observed, candidates) -> MatchResult
```

- **`CutMap`**：候选环状切位图。`length` 为环长（≥1 的 int，拒绝
  bool/float，含 NaN/Inf 专门报错）；`cuts` 为显式给出的上链切位边界，
  必须是 `0..length-1` 内**严格升序、不重复**的整数（不做静默排序/去重，
  非法时下标可定位）；`label` 为可选非空用途的任意字符串标签（可空）。
  0 个切口表示未切环。
- **`circular_fragment_lengths`**：纯函数。切口 `p0<p1<…<p(k-1)` 的完全
  消化循环间隙为 `(p1-p0, …, length-p(k-1)+p0)`；0 切口返回 `(length,)`。
  全部为精确整数，无浮点/容差。
- **`map_from_digest`**：从**环状** `DigestResult` 提取去重物理断裂的
  上链边界，构造 `CutMap`（复用旧切割规则）；线性结果抛 `CandidateError`。
- **`match_candidates(observed, candidates)`**：
  - `observed`：测得的完全消化片段长度可迭代对象（多重集：顺序无关、
    允许重复），每个为正整数，且不能为空（未切环也要给 `[环长]`）。
  - `candidates`：`CutMap` 可迭代对象；先全部物化并完成校验，再开始
    匹配，任一非法整体抛 `CandidateError`。纯函数、无内部状态，失败即
    异常返回，不存在半更新。
  - 匹配判据：候选长度多重集与 `observed` 的 `Counter` **精确相等**，
    且环长等于长度总和（后者被前者蕴含，作为早筛）。
  - **旋转等价去重**：环旋转只改绝对坐标、不改循环间隙序列；以“间隙
    循环序列的字典序最小移位”（键首附切口数）为归一化键，同一等价类
    只保留输入中先出现的候选。未切环（0 切口）与单切环（1 切口）间隙
    同为 `(length,)` 但切口数不同，**不**合并，二者都会匹配 `[length]`。
- **`MatchResult`**：`matches`（去重后匹配候选，保持输入顺序）、
  `examined`（去重前检查总数）、`unique_classes`（全部候选去重后的
  等价类数）；属性 `matched` 表示是否非空。“无匹配”是穷举后被证明的
  结果，无预算/超时路径。

手算示例（n=12，observed=`[4,2,4,2]`）：

| 候选 cuts | 循环间隙 | 结论 |
|---|---|---|
| `(0,2,4,8)` | `(2,2,4,4)` | 匹配 |
| `(1,3,5,9)` | `(2,2,4,4)` | 与上者旋转等价，去重 |
| `(0,2,6,8)` | `(2,4,2,4)` | 多重集同、排列不同，**不**等价，独立匹配 |
| `(0,3,6,9)` | `(3,3,3,3)` | 不匹配 |

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
  - `SequenceError`：底物非法，错误信息含**从 0 开始的位置下标**。
  - `DigestUnsupportedError`：输入合法，但几何上不属于支持的简单片段。
  - `CandidateError`：候选切位图/片段长度非法（环长、切口越界或非升序、
    非正片段长度、成员类型错误、线性 digest 导出等），信息含可定位下标。
- 整数参数（`top_cut`/`bottom_cut`/切位图环长与切口/测得长度）拒绝
  `bool`；任何 `float` 都拒绝，其中 `NaN`/`±Infinity` 有专门报错信息。
- “没有识别位点”“没有匹配候选”都是穷举后**被证明的无解**（前者返回未切
  分子，后者返回空 `matches`），与任何预算耗尽无关；酶切搜索复杂度
  O(N·L)（N≤5000，L≤12），匹配为候选数的线性扫描加小间隙序列归一化，
  无预算/超时失败路径。
- 全部结果为整数与字符串，**没有浮点输出**，故无误差口径/容差问题。

## 已知限制（本契约的有限子集，非行业标准全兼容）

- 底物只接受具体碱基 `ACGT`；IUPAC 简并仅用于酶识别串；不支持 U、
  小写自动转换、RNA、修饰碱基或缺口/泡状底物。
- 不模拟酶切反应动力学、星号活性、甲基化阻断、温度/缓冲液条件，
  也不输出片段浓度；只做确定性的完全酶切几何归属。
- 交错切口互相侵入（配对区 ≤ 0）、同一物理切口配对冲突、
  `|错位| ≥ 环长`、线性边缘单链切割等情形一律 unsupported。
- 单条序列上限 5000 bp、识别串上限 12 bp，均为接口硬边界。
- 候选匹配只处理**环状**切位图与**完全消化**长度：不做部分消化、
  混合样本、测量区间/误差反演；切位边界必须显式给出（不按识别串
  反推未知位点）；匹配只按长度多重集，不使用片段序列或端型信息。

## 测试与独立参考说明

- 测试用例的期望值全部手算，或由测试文件内独立的小函数（`ref_revcomp`、
  显式列窗、独立间隙公式 `indep_gaps`、n=6 全部 64 张切位图的旋转轨道
  全枚举 `indep_orbit` 等）给出，**不调用被测核心生成期望值**；
  链长守恒、沃森-克里克配对等往返/不变量检查只作为补充
  （`TestRoundTripInvariants`）。
- 覆盖：平端、5′/3′ 黏端、单切环状（平端/黏端再环化）、跨零识别、
  反向酶坐标转换、回文单次计数、重复位点、重复物理切口去重、不切、
  IUPAC（含 N）、各类 unsupported 与全部输入校验；候选匹配覆盖正常
  匹配、重复片段、排列不同但多重集相同的非等价候选、旋转去重保留顺序、
  未切环/单切环边界、`digest` 联动、无匹配、全枚举交叉核对以及各类
  非法输入（含 bool/float/NaN/Inf、原子失败与输入不被修改）。
