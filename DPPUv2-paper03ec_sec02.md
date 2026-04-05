## 2. EC-Weyl 結合の構造

本節では、Einstein-Cartan 接続の基本構造を導入し、EC-Weyl scalar $C^2_{\rm EC}$ が torsion モードに対してどのような依存性を示すかを確定する。主要結果は **定理 1 (AX/VT dropout)**、**定理 2 (Palatini 保護 EC 版)**、**定理 3 (δ-sector null theorem)** であり、前二者が EC-Weyl 結合の選択則と torsion の非伝播性を与え、後者が auxiliary な MIXING sector に新たな homogeneous cubic channel が生じないことを与える。

### 2.1 EC 接続と contortion

Einstein-Cartan 理論では、接続 $\omega^{a}{}\_{b}$ を計量と独立な変数として扱い、Levi-Civita (LC) 接続 $\Gamma^{a}{}\_{bc}$ との差を contortion $K^{a}{}_{bc}$ で表す：

$$
\omega_{\rm EC}^a{}_{bc} = \Gamma^a{}_{bc} + K^a{}_{bc}.
$$

Contortion は反対称な torsion テンソル $T^{a}{}\_{bc} = 2\omega^{a}{}\_{[bc]}$ から

$$
K_{abc} = \frac{1}{2}\bigl(T_{abc} + T_{bca} - T_{cab}\bigr)
$$

として定義され、 $K_{abc} = -K_{acb}$ を満たす。

均質 minisuperspace の torsion ansatz として、本論文では 3 つのモードを用いる：

- **AX モード**: 軸性スカラー $\eta$（torsion のトレースレス反対称部分）。 $V = 0$ を前提とする。
- **VT モード**: ベクトル-トレース混合スカラー $V$（torsion のトレース部分）。 $\eta = 0$ を前提とする。
- **MX モード**: AX と VT の共存、すなわち $V \neq 0, \eta \neq 0$。

作用は

$$
S = \int \left[\frac{R_{\rm EC}}{2\kappa^2} + \theta \cdot \mathcal{L}_{\rm NY} + \alpha \cdot C^2_{\rm EC}\right]\sqrt{g}\,d^4x
$$

と書かれ、 $R_{\rm EC}$ は EC Ricci scalar、 $\mathcal{L}\_{\rm NY}$ は Nieh-Yan 密度、 $C^2\_{\rm EC}$ は EC Weyl scalar である。パラメータ $\alpha \in \mathbb{R}$（Weyl coupling）, $\theta \in \mathbb{R}$（NY coupling）を持つ。

均質背景での有効ポテンシャル $V_{\rm eff}(r, \eta, V; \alpha, \theta)$ は $S^3\times S^1$, $T^3\times S^1$, $Nil^3\times S^1$ それぞれについて陽に計算される（スクリプト: `proofs/ax_vt_dropout_proof.py`）。

### 2.2 定理 1 (AX/VT dropout)

EC Weyl scalar の torsion モード依存性を確定する基本定理を述べる。

**定理 1 (AX/VT dropout)**: $S^3\times S^1$, $T^3\times S^1$, $Nil^3\times S^1$ の全 3 トポロジーにおいて、以下が SymPy による厳密な記号恒等式として成立する：

1. **AX モード** ($V = 0$, $\eta \neq 0$):
   $$C^2_{\rm EC} = C^2_{\rm LC} \quad \text{（EC 補正なし）.}$$

2. **VT モード** ($\eta = 0$, $V \neq 0$):
   $$C^2_{\rm EC} = C^2_{\rm LC} \quad \text{（EC 補正なし）.}$$

すなわち、AX および VT の単純背景では EC-Weyl coupling $\alpha C^2_{\rm EC}$ の効果は Levi-Civita の場合と完全に一致し、 $\alpha$ に依存する新寄与は現れない。

**証明の概略**: 各トポロジーで $C^2_{\rm EC} - C^2_{\rm LC}$ を記号的に展開して因数分解すると、EC 補正項は一様に $V\eta$ を因子にもつ。したがって AX モード ($V=0$) では補正は厳密にゼロとなり、VT モード ($\eta=0$) でも同様に消える（SymPy exact zero による記号的確認による）。

**各トポロジーの具体形**:

$S^3\times S^1$ および $T^3\times S^1$ では：
$$C^2_{\rm LC}(S^3, T^3) = 0 \quad (\text{等方点, AX/VT 背景}),$$
AX dropout により $C^2_{\rm EC} = 0$ も確認される。

$Nil^3\times S^1$ では：
$$C^2_{\rm LC}(Nil^3) = \frac{4}{3R^4} \neq 0 \quad (\text{非共形平坦性}),$$
AX dropout により $C^2_{\rm EC} = 4/(3R^4) = C^2_{\rm LC}$ が保たれる（EC による補正は無い）。

**T³ squash+shear null test**: $T^3$ においては、squash パラメータ $\varepsilon$ と shear パラメータ $s$ を両方加えた場合も $C^2_{\rm EC} = 0$ が $\forall \varepsilon, s$ で成立する（ $C^a_{bc} = 0$ による形式的消滅）。

以上の 3 トポロジー × 2 背景（AX, VT）= 6 ケースは、いずれも SymPy による厳密なゼロ恒等式として確認した。詳細な記号計算は [スクリプト: `proofs/ax_vt_dropout_proof.py`, `proofs/t3_flatness_null_test.py`] に対応する。

### 2.3 MX モード: EC 固有の Weyl 結合

定理 1 の裏面として、MX 背景（ $V \neq 0, \eta \neq 0$ ）においては EC 固有の Weyl 結合が現れる。

$S^3\times S^1$ および $T^3\times S^1$ の MX モード：
$$C^2_{\rm EC} = C^2_{\rm LC} + \frac{16V^2\eta^2}{3r^2}.$$

$Nil^3\times S^1$ の MX モード：
$$C^2_{\rm EC} = C^2_{\rm LC} + \frac{16V^2\eta^2}{3R^2} = \frac{4}{3R^4} + \frac{16V^2\eta^2}{3R^2}.$$

いずれも追加項は $V\cdot\eta$ の 2 乗に比例し、torsion の混合積が EC-Weyl 結合を通じて有効ポテンシャルを変形する。具体的には

$$
\Delta V_{\rm eff}^{\rm EC} = \alpha \cdot \frac{16V^2\eta^2}{3r^2} \cdot \text{Vol}(M^3\times S^1)
$$

が LC ポテンシャルに加わる（ $\text{Vol}$ は内部空間の体積因子）。

この EC 固有項は MX 背景における有効ポテンシャル変形の起源となる。一方、§5 で扱う $Nil^3$ の EC false vacuum は、 $\eta = V = 0$ 断面での $C^2_{\rm LC}(Nil^3)\neq 0$ という非共形平坦性に由来する。 $S^3, T^3$ では等方点周りで $V = \eta = 0$ が保たれるため、MX 項はその断面では実効的に不活性である（§3.4, §4.3 参照）。

**物理的解釈**: AX/VT dropout は、torsion 単体が Weyl テンソルの新たな寄与を生まないことを意味する。EC-Weyl 相互作用の本質は torsion の「混合積」 $V\times\eta$ にあり、これは幾何学的には contortion の asymmetric な組み合わせが Weyl 成分を活性化する機構として理解できる。3 トポロジーで共通したこの選択則は、EC 接続の一般的な代数構造から来る普遍的な性質である。

### 2.4 定理 2 (Palatini 保護 EC 版)

**定理 2 (Palatini 保護 EC 版)**: Einstein-Cartan 接続のもとでも

$$
G_{\eta\eta}^{\rm EC} = G_{VV}^{\rm EC} = 0
$$

が厳密に成立する。ここで $G_{AB}^{\rm EC} = \partial^2 \mathcal{L}_{\rm EC}/\partial(\partial_z\Phi^A)\partial(\partial_z\Phi^B)$ は minisuperspace の場の空間計量（kinetic metric）を表す。

**証明の概略（速度法）**: §2.3 で得た MX 背景の EC Weyl scalar

$$
C^2_{\rm EC}(S^3, T^3, \text{MX}) = \frac{16V^2\eta^2}{3r^2}
$$

は $(r,\eta,V)$ の代数関数であり、速度場 $\partial_z\eta$, $\partial_z V$ を含まない。したがって速度法で $G_{\eta\eta}, G_{VV}$ を計算しても、 $\Delta z \to 0$ 極限で消える項しか現れず、torsion の kinetic term は生成されない。詳細な式展開と各トポロジーでの確認は [スクリプト: `proofs/palatini_protection_proof.py`] に対応する。

**r-stiffness の不変性**:

| トポロジー | $G_{rr}$ (LC) | $G_{rr}$ (EC) | 差分 |
|---|---|---|---|
| $S^3$ | 59.22 | 59.22 | 0 |
| $T^3$ | 4675.6 | 4675.6 | $\approx 0$ |
| $Nil^3$ | 5195.2 | 5195.2 | $\approx 0$ |

spin-2 重力子の kinetic 構造も EC で不変である。

### 2.5 定理 3 (δ-sector null theorem)

EC 拡張では、physical な $(\eta, V)$ channel に

$$
\frac{\partial^3 V_{\rm eff}}{\partial\alpha\,\partial\eta\,\partial V} \neq 0
$$

という新しい homogeneous cubic coefficient が現れる（§3.4）。これに対して、auxiliary な MIXING 変数 $\delta_k$ にも独立の cubic channel が生じるかどうかは、別途確認しておく必要がある。

そこで unified engine の formal MIXING sector に対して

$$
C_\delta(M^3;{\rm mode}) \;:=\;
\left.
\frac{\partial^3 V_{\rm eff}}{\partial\delta_0\,\partial\delta_1\,\partial\delta_2}
\right|_{\delta_0=\delta_1=\delta_2=0}
$$

を定義する。ここで $M^3 = S^3, T^3, Nil^3$、mode $=$ AX/VT/MX である。

**定理 3 (δ-sector null theorem)**: 全 3 トポロジー、かつ AX/VT/MX の全背景に対して

$$
C_\delta(M^3;{\rm mode}) = 0
$$

が SymPy の厳密計算として成立する。したがって

$$
\frac{\partial C_\delta}{\partial \alpha} = 0
$$

も自動的に従う。すなわち、EC 拡張は auxiliary な $\delta$-sector に新たな homogeneous cubic channel を生成しない。

**証明の概略**: unified engine で `fiber_mode = MIXING` を有効化し、各トポロジー・各 torsion 背景に対して $V_{\rm eff}$ を構成し、 

$$
\partial_{\delta_0}\partial_{\delta_1}\partial_{\delta_2}V_{\rm eff}|_{\delta=0}
$$

を直接記号微分する。9 ケースすべてで係数は exact zero となり、したがって $\partial\_\alpha C_\delta$ も exact zero である。ここで重要なのは、この結果が「EC 固有効果が存在しない」ことではなく、EC 拡張で新たに活性化する cubic 構造が auxiliary な $\delta$ -sector ではなく physical な $\eta$ - $V$ channel に局在していることを意味する点である。

**3 トポロジーでの帰結**:

- $S^3$ : formal MIXING sector は非自明に導入できるが、 $C_\delta(S^3)=0$ 。したがって §3.4 の $\alpha$ -cubic は $\delta$ -sector ではなく $\eta$ - $V$ sector にのみ現れる。
- $T^3$ : 平坦性のもとで formal MIXING test を行っても $C_\delta(T^3)=0$ 。auxiliary cubic channel は現れない。
- $Nil^3$ : 非共形平坦性と EC false vacuum の存在にもかかわらず、formal MIXING sector では $C_\delta(Nil^3)=0$ 。 $Nil^3$ の EC 固有物理は background Weyl と $(r,\eta,V)$ sector に現れ、 $\delta$-sector には現れない。

詳細な式展開と 9 ケースの exact check は [Appendix D](paper03ec_appD.md) および [スクリプト: `proofs/c_delta_null_proof.py`] に譲る。
