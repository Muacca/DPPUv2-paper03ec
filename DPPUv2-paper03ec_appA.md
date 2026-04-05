## Appendix A. 記号・規約

本付録では、本論で用いる記号規約、3 トポロジーの幾何学的パラメータ、EC 接続の定義、および計算パイプラインをまとめる。

### A.1 基本記号・規約

**指数規約**:
- 4 次元時空指数: $\mu, \nu = 0, 1, 2, 3$（0 = 内部 $S^1$, 1,2,3 = 内部 3 次元空間）
- 空間テンソル指数: $i, j = 1, 2, 3$（内部 3 次元空間）
- spin-2 五重項ラベル: $A = 1, \ldots, 5$（traceless symmetric）
- spin-1 三重項ラベル: $k = 0, 1, 2$（ベクトル）
- トポロジー指数: $S^3$, $T^3$, $Nil^3$

**幾何学的量**:
- $r$: 内部空間のスケールパラメータ（ $S^3$, $T^3$ 共通）
- $R$: Nil³ のスケールパラメータ
- $L$: $S^1$ の半径
- $\kappa$: 重力定数（ $\kappa^2 = 8\pi G$ ）

**体積因子**:

$$
\text{Vol}(S^3\times S^1) = 2\pi^2 L r^3, \quad \text{Vol}(T^3\times S^1) = (2\pi)^4 L r^3, \quad \text{Vol}(Nil^3\times S^1) = (2\pi)^4 L R^3.
$$

**有効ポテンシャル**: minisuperspace 還元後の有効ポテンシャルを

$$
V_{\rm eff} = -\left(\frac{R_{\rm EC}}{2\kappa^2} + \theta\,\mathcal{L}_{\rm NY} + \alpha\,C^2_{\rm EC}\right)\cdot\text{Vol}
$$

で定義する（符号規約: $V_{\rm eff} > 0$ が安定化方向）。

**Pontryagin 密度**: $P = \langle R, \star R\rangle = \frac{1}{2}\varepsilon^{abcd}R_{abef}R_{cd}{}^{ef}$。
**auxiliary MIXING cubic coefficient**:

$$
C_\delta := \left.\frac{\partial^3 V_{\rm eff}}{\partial\delta_0\,\partial\delta_1\,\partial\delta_2}\right|_{\delta=0}.
$$

### A.2 3 トポロジーの構造定数と接続

#### A.2.1 $S^3\times S^1$

$S^3$ の構造定数（Maurer-Cartan 形式 $\sigma^i$ に関する正規直交基底 $e^i = r\sigma^i/2$）:

$$
C^i_{jk}(S^3) = \frac{4}{r}\varepsilon_{ijk} \quad (i, j, k = 1, 2, 3).
$$

**重要**: 係数は $4/r$（ $2/r$ ではない）。これは $e^i = r\sigma^i/2$ の規格化から来る（付録 B の参照）。

Koszul 公式による接続係数:

$$
\omega_{abc} = \frac{1}{2}\bigl(C_{abc} + C_{cba} - C_{bac}\bigr).
$$

体積保存 squash+shear ansatz（ $\det = 1$ ）:

$$
e^0 = r(1+\varepsilon)^{2/3}(1+s)^2\,\sigma^0/2, \quad e^1 = r(1+\varepsilon)^{2/3}(1+s)^{-2}\,\sigma^1/2, \quad e^2 = r(1+\varepsilon)^{-4/3}\,\sigma^2/2.
$$

#### A.2.2 $T^3\times S^1$

$T^3 = \mathbb{R}^3/\mathbb{Z}^3$ の構造定数:

$$
C^a_{bc}(T^3) = 0 \quad \forall\, a, b, c.
$$

体積非保存 squash+shear ansatz:

$$
e^1 = r(1+\varepsilon)\,dx^1, \quad e^2 = r(1+s)\,dx^2, \quad e^3 = r\,dx^3.
$$

体積因子: $\text{Vol}(\varepsilon, s) = (2\pi)^4 L r^3(1+\varepsilon)(1+s)$（ $\varepsilon, s$ 依存）。

#### A.2.3 $Nil^3\times S^1$

$Nil^3$（3 次元 Heisenberg 多様体）の構造定数:

$$
C^2_{01}(Nil^3) = \frac{1}{R}, \quad \text{（その他はすべてゼロ）.}
$$

本文では Lie 括弧 $[E_0,E_1]=(1/R)E_2$ 規約でこの符号を採る。一方、実装は Maurer-Cartan 形 $d\sigma^2=-\sigma^0\wedge\sigma^1$ を基準にしているため、frame 係数は $C[2,0,1]=-(1+\varepsilon)^{-4/3}(1+s)^{-2}/R$ と書かれる。両者は同じ Heisenberg 構造を表す規約差である。

体積保存 squash+shear ansatz:

$$
e^0 = R\,\sigma^0, \quad e^1 = R(1+\varepsilon)^{2/3}(1+s)\,\sigma^1, \quad e^2 = R(1+\varepsilon)^{-2/3}(1+s)^{-1}\,\sigma^2.
$$

体積保存性の確認:

$$
\det(e^0 \wedge e^1 \wedge e^2) \propto (1+\varepsilon)^{2/3-2/3}(1+s)^{1-1} = 1 \; ✓
$$

LC 背景 Weyl scalar（ $\eta=V=0$ ）: $C^2_{\rm LC}(Nil^3) = 4/(3R^4) \neq 0$（非共形平坦性）。

### A.3 EC 接続と contortion

Einstein-Cartan 接続は LC 接続 $\Gamma^{a}{}\_{bc}$ と contortion $K^{a}{}\_{bc}$ の和：

$$
\omega_{\rm EC}^a{}_{bc} = \Gamma^a{}_{bc} + K^a{}_{bc}.
$$

Contortion の定義（torsion $T^{a}{}\_{bc} = 2\omega^{a}{}\_{[bc]}$ から）:

$$
K_{abc} = \frac{1}{2}\bigl(T_{abc} + T_{bca} - T_{cab}\bigr), \quad K_{abc} = -K_{acb}.
$$

**Torsion モードの定義**（均質 ansatz）:

- **AX（軸性スカラー）** $\eta$: torsion のトレースレス反対称部分。 $S^3$ では $T^{i}{}\_{jk} \propto \eta\,\varepsilon_{ijk}$ 。
- **VT（ベクトル-トレース混合）** $V$: torsion のトレース部分。 $S^3$ では $T^i{}_{0j} \propto V\,\delta^i_j$。
- **MX モード**: AX と VT の共存（ $\eta \neq 0$ かつ $V \neq 0$ ）。

EC Weyl scalar の MX モードでの追加項（全トポロジー）:

$$
C^2_{\rm EC} - C^2_{\rm LC} = \frac{16V^2\eta^2}{3r^2} \quad (S^3, T^3), \qquad = \frac{16V^2\eta^2}{3R^2} \quad (Nil^3).
$$

### A.4 Pontryagin 保護の機構

**Plücker 保護（AX モード）**: $V=0$ のとき、Riemann テンソルの行列表現が Plücker 関係式を満たし、 $P \equiv 0$。
**Pfaffian 保護（VT モード）**: $\eta=0$ のとき、Pfaffian 型の恒等式により $P \equiv 0$。
### A.5 計算パイプライン

本論文の計算は以下のパイプラインで実行される。

1. **幾何構造の設定**: 3 トポロジーの $e^i$（vielbein）、構造定数 $C^a_{bc}$、体積因子を設定（`script/dppu/topology/`）。
2. **EC 接続の計算**: contortion $K$ を torsion ansatz から構成し、EC 接続 $\omega_{\rm EC}$ を算出（`script/dppu/connection/`）。
3. **Riemann・Weyl テンソルの計算**: EC 接続から曲率テンソル、Ricci scalar $R_{\rm EC}$、Weyl scalar $C^2_{\rm EC}$、Pontryagin 密度 $P$ を SymPy で計算（`script/dppu/curvature/`）。
4. **有効ポテンシャルの構成**: $V_{\rm eff}(r, \eta, V; \alpha, \theta)$ を陽に展開し、Hessian 行列を取得（`script/dppu/action/`）。
5. **モード辞書の確定**: spin-2 quintet（shear, off-diagonal）, spin-1 triplet, spin-0 の質量スペクトルを一般化固有値問題や停留点探索で点検する（`script/scripts/paper03ec/`, `script/scripts/proofs/`）。

**主要スクリプト対応表**:

| 内容 | スクリプト |
|---|---|
| $C^2_{\rm EC}$ 記号計算（AX/VT/MX dropout） | `script/scripts/proofs/ax_vt_dropout_proof.py` |
| T³ mode dictionary（squash/shear null） | `script/scripts/paper03ec/t3_mode_dictionary.py` |
| Nil³ EC mode dictionary | `script/scripts/paper03ec/nil3_mode_dictionary.py` |
| T³ MX 真空探索 | `script/scripts/paper03ec/t3_mx_vacuum_search.py` |
| S³ MX 停留点探索 | `script/scripts/paper03ec/s3_mx_vacuum_search.py` |
| Nil³ MX 停留点探索 | `script/scripts/paper03ec/nil3_mx_vacuum_search.py` |
| Nil³ EFT at EC false vacuum | `script/scripts/paper03ec/nil3_false_vacuum_eft.py` |
| S³ VT spin-2/1 質量 | `script/scripts/paper03ec/s3_vt_spin_masses.py` |
| Nil³ quintet 分裂 | `script/scripts/paper03ec/nil3_spin2_quintet_splitting.py` |
| ε-s cross-term 3 トポロジー比較 | `script/scripts/paper03ec/squash_shear_cross_term.py` |
| Palatini 保護（G_ηη=G_VV=0） | `script/scripts/proofs/palatini_protection_proof.py` |
| auxiliary $C_\delta$ null theorem | `script/scripts/proofs/c_delta_null_proof.py` |
| Cubic channel 構造 | `script/scripts/paper03ec/ec_cubic_vertex.py` |
| γ=1/2 解析的証明 | `script/scripts/proofs/gamma_scaling_proof.py` |
| T³ flatness null test | `script/scripts/proofs/t3_flatness_null_test.py` |

### A.6 Sanity checks

本論文の主要命題は以下の sanity check を経て採用している。

1. AX/VT dropout: $\Delta C^2 = 0$（6 ケース, SymPy exact zero）。
2. Palatini 保護: $G_{\eta\eta} = G_{VV} = 0$（速度法, 3 トポロジー）。
3. auxiliary $\delta$-sector: $C_\delta = 0$ かつ $\partial_\alpha C_\delta = 0$（3 トポロジー × 3 背景, SymPy exact）。
4. Palatini 普遍性: VT $-$ AX spin-2/1 質量 $= 0$（SymPy exact）。
5. γ=1/2: $r_0^{\rm 解析}$ vs $r_0^{\rm 数値}$ の誤差 $< 0.001$%（全 6 点）。
6. quintet zero mode: $m^2(q_3) = 0$（SymPy, cross terms も $= 0$）。
7. 体積保存性: $S^3, Nil^3$ の squash+shear で $\det = 1$（SymPy）。
8. cross-term 起源確認: $\partial^2 V/\partial\varepsilon\partial s$ の $\alpha$ 依存性検証（ $T^3$ : 独立, $Nil^3$ : $\alpha$ 依存）。
9. S³ MX 停留点探索: 辞書対象範囲では Hessian 正定値の full stationary vacuum が存在しない。
10. Nil³ MX 停留点探索: $\alpha<0$ 分枝では $V\eta\neq 0$ を満たす実数停留点が存在せず、残る局所極小は $\eta=V=0$ false vacuum のみ。
