## 4. Mode dictionary: $T^3\times S^1$ under EC-Weyl

本節では $T^3\times S^1$ の均質幾何自由度に対する EC-Weyl 結合の効果を解析する。 $T^3$ の平坦性（構造定数 $C^a_{bc} = 0$）が引き起こす三重の triviality と、squash-shear cross-term の kinematic 起源（**定理 5**）を確立する。

### 4.1 T³ 平坦性定理

$T^3 = \mathbb{R}^3/\mathbb{Z}^3$ は平坦なトーラスであり、構造定数

$$
C^a_{bc}(T^3) = 0 \quad \forall\, a, b, c
$$

が恒等的にゼロとなる。この平坦性は Weyl テンソルの消滅を通じて spin-2/spin-1 の全質量をゼロにする。

**spin-2 質量（squash 方向）**: SymPy 厳密計算により

$$
\frac{\partial^2 V}{\partial\varepsilon^2}\bigg|_{\varepsilon=0} = 0 \quad \forall\, \eta, V, \alpha
$$

が確認される。 $C^a_{bc} = 0$ より Weyl テンソルが消えるため、spin-2 に対する幾何的質量項が生じない。[スクリプト: `paper03ec/t3_mode_dictionary.py`]

**spin-2 質量（shear 方向）**: 同様に

$$
\frac{\partial^2 V}{\partial s^2}\bigg|_{s=0} = 0 \quad \forall\, \eta, V, \alpha
$$

が成立する（SymPy exact zero）。[スクリプト: `proofs/t3_flatness_null_test.py`]

**spin-1 質量**: $C^a_{bc} = 0$ により、全 spin-1 成分は質量ゼロとなる（ $m^2_\omega = 0$ ）。

### 4.2 Mode dictionary table (Table 2: T³×S¹)

**Table 2: $T^3\times S^1$ mode dictionary under EC-Weyl**

| 物理量 | AX ($V=0, \eta\neq 0$) | VT ($\eta=0, V\neq 0$) | MX ($V\neq 0, \eta\neq 0$) |
|---|---|---|---|
| $C^2_{\rm EC}$ | $= C^2_{\rm LC} = 0$ | $= C^2_{\rm LC} = 0$ | $16V^2\eta^2/(3r^2)$ |
| AX/VT dropout | ✓ | ✓ | ✗ |
| spin-0: $r$ | 力学的 | 力学的 | 力学的 |
| spin-0: $\eta, V$ | $\eta$ 非力学的, $V=0$ | $\eta=0$, $V$ 非力学的 | 両方非力学的 |
| spin-2 質量 $\mu^2$ | $0$（ $T^3$ 平坦性） | $0$（ $T^3$ 平坦性） | $0$（ $T^3$ 平坦性） |
| spin-1 質量 $m^2_\omega$ | $0$（質量ゼロ） | $0$ | $0$ |

### 4.3 T³ の三重 triviality

$T^3\times S^1$ における EC-Weyl 結合の効果は 3 つの意味で trivial となる。

**① AX dropout**: $C^a_{bc} = 0$ より $C^2_{\rm LC} = 0$、定理 1 の AX dropout により $C^2_{\rm EC} = C^2_{\rm LC} = 0$。 $\alpha$ は有効ポテンシャルを変えない（AX 背景）。

**② VT dropout**: 同様に $\eta = 0$ より定理 1 の VT dropout が適用され、 $C^2_{\rm EC} = 0$。VT 背景でも $\alpha$ は不活性。

**③ MX 補正**: AX/VT dropout が外れる MX 背景（ $V\eta \neq 0$ ）では $C^2_{\rm EC} = 16V^2\eta^2/(3r^2) \neq 0$ となる。 $\alpha < 0$ ではこの補正は $\Delta V_{\rm eff}^{\rm EC} > 0$ を与えてポテンシャルを押し上げる。有限な MX 臨界点の有無は §4.5 で代数的に判定する。

formal MIXING sector に対する auxiliary EFT 係数も

$$
C_\delta(T^3) = 0,
\qquad
\frac{\partial C_\delta(T^3)}{\partial\alpha} = 0
$$

を満たす（定理 3, App. D）。したがって $T^3$ では、physical sector でも auxiliary sector でも EC 拡張による新しい homogeneous cubic channel は現れない。

### 4.4 定理 5 (torsion-volume coupling)

**定理 5 (torsion-volume coupling)**: $T^3\times S^1$ において、squash パラメータ $\varepsilon$ と shear パラメータ $s$ の混合 2 次微分は

$$
\frac{\partial^2 V_{\rm eff}}{\partial\varepsilon\,\partial s}\bigg|_{\varepsilon=s=0} = \frac{48\pi^4 L\eta^2 r}{\kappa^2} \neq 0 \quad (\eta \neq 0)
$$

を満たす。この cross-term の起源は体積素 $(1+\varepsilon)(1+s)$ の kinematic mixing であり、Weyl 質量とは無関係である。

**証明の概略**: $T^3$ の有効ポテンシャルは

$$
V_{\rm eff}(T^3) = \frac{16\pi^4 r}{3}\left[V^2 r^2 + 6V\eta r\theta + 9\eta^2 - 16\alpha V^2\eta^2\right]
$$

として書かれる（AX 背景では $V=0$）。squash と shear を加えると体積因子として $(1+\varepsilon)(1+s)$ が積として現れ、この因子が $\partial^2/\partial\varepsilon\partial s$ に対して非零寄与を与える。具体的には

$$
V_{\rm eff}^{T^3}(\varepsilon, s) \supset \frac{16\pi^4 r}{\kappa^2}\cdot 9L\eta^2 \cdot (1+\varepsilon)(1+s) + \cdots
$$

から $\partial^2 V/\partial\varepsilon\partial s = 48\pi^4 L\eta^2 r/\kappa^2$ が導かれる。

**Weyl 質量でないことの確認**: AX 背景（ $V=0$ ）では $C^a_{bc}=0$ により $C^2_{\rm LC}=0$、さらに定理 1（AX dropout）により $C^2_{\rm EC} = C^2_{\rm LC} = 0$ が成立するため、この cross-term は $\alpha$ に独立であり Weyl 質量への寄与はゼロ。数値確認（AX 背景）： $\alpha=0$ と $\alpha=-1$ で cross-term 値が同一（SymPy exact）。なお MX 背景（ $V\neq 0$ ）では $C^2_{\rm EC} = 16V^2\eta^2/(3r^2) \neq 0$ となり、cross-term に $\alpha$ 依存項が加わるが、spin-2 null test（ $\partial^2 V/\partial\varepsilon^2 = 0$, $\partial^2 V/\partial s^2 = 0$ ）は依然成立する（ $C^2_{\rm EC}$ が $\varepsilon, s$ に非依存なため）。[スクリプト: `proofs/t3_flatness_null_test.py`, `paper03ec/squash_shear_cross_term.py`]

### 4.5 MX 臨界点の不存在

$\alpha < 0$ での $T^3$ MX sector の臨界点を解析する。

$\partial V_{\rm eff}/\partial\eta = \partial V_{\rm eff}/\partial V = 0$ の拘束面を SymPy で解くと、 $T^3$ の MX 解は

$$
V^* = \pm\frac{3}{4}\sqrt{\frac{1}{\alpha}}, \quad \eta^* = \pm\frac{r}{4}\sqrt{\frac{1}{\alpha}}
$$

であり、 $\alpha < 0$ では純虚数になる（実数解なし）。したがって **$T^3$ には $\alpha < 0$ のもとで有限な実 MX 臨界点が存在しない**。

これは §4.3 の `MX 補正はポテンシャルを押し上げる` という解析と整合しており、 $T^3$ の MX sector が新たな均質局所極小を与えないことを示す。[スクリプト: `paper03ec/t3_mx_vacuum_search.py`]
