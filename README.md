# 限制性内切酶的环状片段归属（restriction_digest）

离线、纯标准库的限制性酶切片段归属库：给定线性或环状双链 DNA 与一个/多个
限制酶定义，搜索正向与反向互补识别位点，转换并取模切口坐标，对重复物理切口
去重，返回**按上链切口分段**的双链片段、上/下链边界以及 5′/3′ 突出端序列；
当交错切口造成无法做无交叉简单归属的片段时，明确抛出 `DigestUnsupportedError`，
绝不猜测配对。此外可对环状分子计算一次独立随机部分酶切后，各片段长度的
期望条数；不枚举所有切点子集。

## 运行环境与隔离 venv

- Python **3.14.7**（仅使用标准库 `dataclasses`、`math`、`typing`）。
- Linux 原生运行；装好 Python 后**无需安装任何第三方依赖**，无需网络、
  账号、密钥、数据库、中间件或 Docker。

如需完全隔离环境，可在仓库根目录执行：

```bash
python3.14 -m venv .venv
. .venv/bin/activate
python -m pip --version
```

不创建 venv 也可直接运行；本项目不需要安装任何第三方包。

## 运行命令（在仓库根目录执行）

```bash
python -m unittest discover -s tests -v   # 基线回归与部分酶切新增测试
python demo.py                            # 基线完全酶切几何演示
python demo_iteration.py                  # 部分酶切期望演示（远低于 90 秒）
```

## 公开接口

```python
from restriction_digest import Enzyme, digest, compatible

ecori = Enzyme(name="EcoRI", recognition="GAATTC", top_cut=1, bottom_cut=5)
result = digest("GGGAATTCCC", ecori, topology="linear")   # 或 "circular"
```

部分酶切：

```python
from restriction_digest import partial_digest_expectation

expectation = partial_digest_expectation(
    "GGGAATTCCCAA",
    ecori,
    lambda site: 0.3,
)
for row in expectation.expected_fragments:
    print(row.length, row.kind, row.expected_count)
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

### `partial_digest_expectation(sequence, enzymes, cut_probabilities, topology="circular")`

计算**环状分子一次部分酶切**后的片段长度期望分布，当前仅支持环状拓扑。
概率赋予对象是 `digest` 识别出的每条 `CutSite` 证据：

- 可传统一概率：`0.3`；
- 可传函数：`lambda site: 0.2 + 0.1 * (site.start % 3)`；
- 可传等长序列，顺序与 `DigestResult.cut_sites` 一致；
- 概率必须是有限 `int`/`float`，范围 `[0,1]`，拒绝 `bool`、`NaN`、无穷。

位置完全重合（`(top_boundary, bottom_boundary)` 相同）的多条识别证据先合并：

```text
P(合并后切点切开) = 1 - ∏(1 - P_i)
```

即同一物理位置上“至少一个识别事件切开”的概率；合并后的切点之间再视为相互
独立。合并点只产生一次事件，因此不会产生零长片段。

`PartialDigestExpectation`：

- `length`：环长；
- `merged_cut_sites`：按上链边界排序的合并切点，含 `probability` 与 `evidence`；
- `expected_fragments`：按 `(length, kind)` 聚合；`"circular"` 仅表示零切开后
  保留的一条整环，`"linear"` 表示至少一个切点切开后的线性片段；
- `expected_total_length`：`Σ 长度 × 期望条数`，应在浮点舍入误差内等于环长。

只切一次时得到一条全长线性片段；它与零切开的全长环类型分开报告，即使二者
长度都等于环长。两个实际切点 `i→j` 沿环形成一段片段的概率为：

```text
P_i × P_j × ∏(中间实际切点 k 未切)
```

环外其他切点是否切开不影响该弧段是否存在，因此算法只做有序切点对的 O(m²)
累计；零切开用 `∏(1-P_i)`，单切开用前缀/后缀乘积 O(m) 计算，不枚举 2^m。

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
- 整数参数（`top_cut`/`bottom_cut`）拒绝 `bool`；任何 `float` 都拒绝，
  其中 `NaN`/`±Infinity` 有专门报错信息。
- “没有识别位点”是穷举搜索后**被证明的无解**（返回未切分子），与任何
  预算耗尽无关；搜索复杂度 O(N·L)（N≤5000，L≤12），无预算/超时失败路径。
- `digest` 的完全酶切结果为整数与字符串，**没有浮点输出**；部分酶切
  概率接口的浮点误差口径见下节。

## 误差口径与已知限制（本契约的有限子集，非行业标准全兼容）

- 底物只接受具体碱基 `ACGT`；IUPAC 简并仅用于酶识别串；不支持 U、
  小写自动转换、RNA、修饰碱基或缺口/泡状底物。
- 不模拟酶切反应动力学、星号活性、甲基化阻断、温度/缓冲液条件，
  也不输出片段浓度；完全酶切接口只做确定性的几何归属。
- 部分酶切接口只输出给定独立切开概率下的期望条数，不模拟真实实验、不采样，
  也不从概率推导动力学参数；输入概率的相关性和条件依赖由调用方负责。
- 部分酶切期望为 IEEE 双精度浮点，可能有极小舍入误差；不做人为容差修正，
  测试中的近似比较仅用于处理正常浮点舍入。
- 概率部分酶切当前仅支持环状分子；线性分子的天然末端归属不属于本次增量。
- 交错切口互相侵入（配对区 ≤ 0）、同一物理切口配对冲突、
  `|错位| ≥ 环长`、线性边缘单链切割等情形一律 unsupported。
- 单条序列上限 5000 bp、识别串上限 12 bp，均为接口硬边界。

## 测试与独立参考说明

- 测试用例的期望值全部手算，或由测试文件内独立的小函数（`ref_revcomp`、
  显式列窗等）给出，**不调用被测核心生成期望值**；链长守恒、沃森-克里克
  配对等往返/不变量检查只作为补充（`TestRoundTripInvariants`）。
- 覆盖：平端、5′/3′ 黏端、单切环状（平端/黏端再环化）、跨零识别、
  反向酶坐标转换、回文单次计数、重复位点、重复物理切口去重、不切、
  IUPAC（含 N）、各类 unsupported 与全部输入校验。
- `tests/test_partial_digest.py` 对少量切点显式枚举全部 0/1 切点子集，
  独立核对 O(m²) 公式；覆盖概率全 0、全 1、两个相邻位点、三个相邻位点、
  重合切点合并、总长度守恒与概率校验。
