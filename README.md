# 限制性内切酶的环状片段归属（restriction_digest）

离线、纯标准库的限制性酶切片段归属库：给定线性或环状双链 DNA 与一个/多个
限制酶定义，先穷举全部物理识别位点并给出稳定 ID，再可选择只启用其中一部分。
被启用的位点才会转换并取模上/下链切口、对重复物理切口去重，并返回**按上链
切口分段**的双链片段、边界以及 5′/3′ 突出端序列；当交错切口造成无法做无交叉
简单归属的片段时，明确抛出 `DigestUnsupportedError`，绝不猜测配对。

## 运行环境

- Python **3.14.7**（仅使用标准库 `dataclasses`、`hashlib`、`math`、`typing`）。
- Linux 原生运行；装好 Python 后**无需安装任何第三方依赖**，无需网络、
  账号、密钥、数据库、中间件或 Docker。

## 运行命令（在仓库根目录执行）

```bash
./.venv/bin/python -m unittest discover -s tests -v
./.venv/bin/python demo.py
```

## 公开接口

```python
from restriction_digest import Enzyme, digest, recognize_sites, compatible

ecori = Enzyme(name="EcoRI", recognition="GAATTC", top_cut=1, bottom_cut=5)
result = digest("GGGAATTCCC", ecori, topology="linear")   # 或 "circular"

recognized = recognize_sites("GGGAATTCCC", [ecori])
partial = digest(
    "GGGAATTCCC", [ecori], topology="linear",
    enabled_site_ids=[recognized.sites[0].site_id],
)
no_cut = digest("GGGAATTCCC", [ecori], enabled_site_ids=())
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

### `recognize_sites(sequence, enzymes, topology="linear") -> RecognitionResult`

识别阶段只做确定性穷举搜索，不执行切口几何计算、不模拟反应动力学，也不会因
启用集合之外的交错切口而抛出 `DigestUnsupportedError`。每个 `RecognitionSite`
表示一个物理识别区间：

- `site_id`：稳定、不可手工构造的 ID，绑定本次 DNA、拓扑和按调用顺序编号的
  酶表；它不是跨不同实验复用的名称。
- `start`/`length`/`matched_sequence`：当前输入上链坐标中的起点、长度和实际
  匹配碱基。
- `strand_sources`：保留 `("+",)`、`("-",)` 或 `("+", "-")` 全部方向来源。
  回文模式两个方向得到相同切口时只切一次；若简并模式两个方向得到不同切口，
  仍按原几何冲突规则拒绝猜测配对。

环状身份约定：ID 以整环字典序最小旋转为基准；旋转同一环状 DNA 后，物理位点
ID 不变，`start` 随当前坐标原点变化。跨零识别仍以循环窗口完成；线性 DNA 不
做首尾拼接。

### `digest(sequence, enzymes, topology="linear", enabled_site_ids=None) -> DigestResult`

- `sequence`：大写 `ACGT` 字符串，长度 1..**5000**。底物不接受简并碱基、
  小写或空白（不做静默纠正）。
- `enzymes`：单个 `Enzyme` 或非空可迭代对象（多酶同时完全酶切）。
- `topology`：`"linear"` 或 `"circular"`。
- `enabled_site_ids=None`：默认启用本次识别到的全部位点，旧公开调用不变。
- `enabled_site_ids=[]` 或 `set()`：显式空集合，返回未切分子。
- 传入 ID 可迭代对象：仅这些物理位点参与原上下链切口计算与片段构建。
- ID 必须属于本次 DNA、拓扑和酶表的识别结果；旧结果或其他输入的 ID 混入时，
  在片段构造前整体抛出 `DigestError`。
- 纯函数、无内部状态；ID 校验和几何登记先于片段构造，失败即整体失败。

`DigestResult`：

- `topology`, `length`；
- `recognition_sites`：本次识别到的全部物理位点，包含未启用位点。
- `cut_sites: tuple[CutSite, ...]`：已启用并实际断裂的位点证据（酶名、起点、
  代表链、取模后的上/下链切口、`site_id`）；`cut_count` 为按
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

## 坐标与归属约定（设计取舍）

1. **上链统一坐标系**：所有切口都是上链 5′→3′ 的边界编号；下链切口也用
   上链坐标表达，避免双坐标制换算错误。环状坐标对 `N` 取模。
2. **物理切口去重**：`(链, 坐标)` 相同的切口只算一次；若同一物理切口被
   两个不同的对侧切口争夺（无法去重为同一末端），报 `DigestUnsupportedError`。
   部分酶切只登记已启用位点；只启用两个冲突位点中的一个时可得到该位点片段，
   全启用仍与旧入口一致地整体拒绝。
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
- 全部结果为整数与字符串，**没有浮点输出**，故无误差口径/容差问题。

## 已知限制（本契约的有限子集，非行业标准全兼容）

- 底物只接受具体碱基 `ACGT`；IUPAC 简并仅用于酶识别串；不支持 U、
  小写自动转换、RNA、修饰碱基或缺口/泡状底物。
- 不模拟酶切反应动力学、星号活性、甲基化阻断、温度/缓冲液条件，
  也不输出片段浓度；只做确定性的识别集合选择与几何归属，不加入随机概率。
- 交错切口互相侵入（配对区 ≤ 0）、同一物理切口配对冲突、
  `|错位| ≥ 环长`、线性边缘单链切割等情形一律 unsupported。
- 单条序列上限 5000 bp、识别串上限 12 bp，均为接口硬边界。

## 测试与独立参考说明

- 测试用例的期望值全部手算，或由测试文件内独立的小函数（`ref_revcomp`、
  显式列窗等）给出，**不调用被测核心生成期望值**；链长守恒、沃森-克里克
  配对等往返/不变量检查只作为补充（`TestRoundTripInvariants`）。
- 覆盖：平端、5′/3′ 黏端、单切环状（平端/黏端再环化）、跨零识别、
  反向酶坐标转换、回文单次计数、重复位点、重复物理切口去重、不切、
  IUPAC（含 N）、部分酶切、空集合、旧/外来 ID 拒绝、环状旋转身份、
  跨原点片段、各类 unsupported 与全部输入校验。
