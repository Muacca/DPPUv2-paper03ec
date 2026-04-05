## Appendix C. ε-s cross-term の詳細導出

本付録では、定理 9（ε-s cross-term 分類）の詳細な証明を与える。squash パラメータ $\varepsilon$ と shear パラメータ $s$ の混合二次微分 $\partial^2 V/\partial\varepsilon\partial s|_{\varepsilon=s=0}$ が 3 トポロジーで本質的に異なる起源を持つことを、体積保存性の代数的確認と SymPy 計算から示す。

### C.1 体積保存性の代数的証明

#### C.1.1 $S^3$ の squash+shear 体積保存

$S^3$ の volume-preserving squash+shear ansatz（付録 A.2.1）では、vielbein の決定因子が

$$
\det(e^1, e^2, e^3) \propto (1+\varepsilon)^{2/3}(1+s)^2 \cdot (1+\varepsilon)^{2/3}(1+s)^{-2} \cdot (1+\varepsilon)^{-4/3}
$$

となる（ $r$ を除くスケール因子の積）。指数を整理すると：

$$
(1+\varepsilon)^{2/3+2/3-4/3}\cdot(1+s)^{2-2} = (1+\varepsilon)^0\cdot(1+s)^0 = 1.
$$

すなわち $S^3$ の体積は $\varepsilon, s$ に完全に非依存である（SymPy で確認）：

$$
\text{Vol}(S^3; \varepsilon, s) = 2\pi^2 L r^3 = \text{const w.r.t. }\varepsilon, s.
$$

#### C.1.2 $Nil^3$ の squash+shear 体積保存

$Nil^3$ の volume-preserving squash+shear ansatz（付録 A.2.3）では：

$$
e^0 = R\,\sigma^0, \quad e^1 = R(1+\varepsilon)^{2/3}(1+s)\,\sigma^1, \quad e^2 = R(1+\varepsilon)^{-2/3}(1+s)^{-1}\,\sigma^2.
$$

体積因子の $\varepsilon, s$ 依存部分：

$$
(1+\varepsilon)^{2/3}\cdot(1+s)\cdot(1+\varepsilon)^{-2/3}\cdot(1+s)^{-1} = 1.
$$

したがって（ $e^0 = R$ は $\varepsilon, s$ 非依存）:

$$
\text{Vol}(Nil^3; \varepsilon, s) = (2\pi)^4 L R^3 = \text{const w.r.t. }\varepsilon, s.
$$

$Nil^3$ でも体積保存が成立する（SymPy 確認）。

#### C.1.3 $T^3$ の体積非保存

$T^3$ の squash+shear ansatz では：

$$
e^1 = r(1+\varepsilon)\,dx^1, \quad e^2 = r(1+s)\,dx^2, \quad e^3 = r\,dx^3.
$$

体積因子は積の形：

$$
\text{Vol}(T^3; \varepsilon, s) = (2\pi)^4 L r^3(1+\varepsilon)(1+s),
$$

これは $\varepsilon, s$ に依存する。特に $\partial^2\text{Vol}/\partial\varepsilon\partial s = (2\pi)^4 L r^3 \neq 0$。

### C.2 T³ cross-term の導出

$T^3$ の AX 背景（ $V=0$ ）での有効ポテンシャルは、squash+shear を含めると

$$
V_{\rm eff}^{T^3}(\varepsilon, s) = \frac{16\pi^4 L r}{3\kappa^2}\cdot 9\eta^2\cdot(1+\varepsilon)(1+s) + (\text{他の項})
$$

と書かれる（ $C^a_{bc}=0$ より Weyl 寄与はなく、torsion kinetic 項 $\propto \eta^2$ が体積因子 $(1+\varepsilon)(1+s)$ と積になる）。

混合二次微分：

$$
\frac{\partial^2 V_{\rm eff}^{T^3}}{\partial\varepsilon\,\partial s}\bigg|_{\varepsilon=s=0} = \frac{16\pi^4 L r}{\kappa^2}\cdot 9\eta^2\cdot\frac{\partial^2}{\partial\varepsilon\partial s}[(1+\varepsilon)(1+s)]\bigg|_{\varepsilon=s=0} = \frac{48\pi^4 L\eta^2 r}{\kappa^2}.
$$

**Weyl 質量でないことの確認**: AX 背景（ $V=0$ ）では $C^a_{bc}=0$ により $C^2_{\rm LC}=0$、定理 1（AX dropout）により $C^2_{\rm EC}=0$ が成立するため、 $\alpha$ 依存項はこの cross-term に現れない。数値検証（AX 背景）： $\alpha=0$ と $\alpha=-1$ で同一値を確認（SymPy exact）。MX 背景では $C^2_{\rm EC} = 16V^2\eta^2/(3r^2)\neq 0$ により cross-term に $\alpha$ 依存項が加わるが（ §4.4 参照）、spin-2 null test は保たれる。

数値例（ $r=3, L=\kappa=\eta=1$ ）:

$$
\partial^2 V/\partial\varepsilon\partial s = 48\pi^4\times 3 \approx 1.403\times 10^4.
$$

[スクリプト: `proofs/t3_flatness_null_test.py`]

### C.3 S³ cross-term の消滅

$S^3$ での cross-term の消滅を 2 段階で示す。

**Step 1（体積保存）**: §C.1.1 より $\text{Vol}(S^3; \varepsilon, s) = \text{const}$。体積因子から来る kinematic cross-term はゼロ。

**Step 2（ $s \to -s$ 対称性）**: $S^3$ の squash+shear 有効ポテンシャルは反射対称性

$$
V_{\rm eff}^{S^3}(\varepsilon, s) = V_{\rm eff}^{S^3}(\varepsilon, -s)
$$

を持つ（ $S^3$ の Weyl テンソルと LC 接続が $s \to -s$ の対称性から来る偶数べき展開のみを含む）。したがって

$$
\frac{\partial^2 V^{S^3}}{\partial\varepsilon\partial s}\bigg|_{\varepsilon=s=0} = 0 \quad (\text{SymPy exact}).
$$

数値確認（ $r=3, L=\kappa=\eta=1, \alpha=0$ ）: $\partial^2 V/\partial\varepsilon\partial s = 0.000000\times 10^{0}$（機械精度）。

### C.4 Nil³ cross-term の Weyl 起源

$Nil^3$ の cross-term は体積保存（§C.1.2）にもかかわらず非零となる。

$Nil^3$ の AX 背景では、体積保存により kinematic 項は cross-term を生まない。SymPy 計算で得られる混合二次微分の最終式は

$$
\frac{\partial^2 V^{Nil^3}}{\partial\varepsilon\,\partial s}\bigg|_{\varepsilon=s=0} = \frac{384\pi^4 LR^2 - 8192\pi^4 L\alpha\kappa^2}{9R\kappa^2}.
$$

**α 依存性による Weyl 起源の確認**： $T^3$ とは異なり $\alpha$ 依存項が現れる。本文では $[E_0,E_1]=(1/R)E_2$ 規約により $C^2_{01}(Nil^3)=1/R$ と書くが、以下の式はコードと同じく Maurer-Cartan 基底 $d\sigma^2=-\sigma^0\wedge\sigma^1$ に対応する frame 係数であるため負号を持つ：

$$
C^2_{01}(\varepsilon, s) = -\frac{(1+\varepsilon)^{-4/3}(1+s)^{-2}}{R}
$$

が $\varepsilon, s$ に関して非対称（特に $(1+s)^{-2}$ の $s$ 依存性）であるため、体積保存にもかかわらず $\partial^2 C^2_{01}/\partial\varepsilon\partial s \neq 0$ となり、Weyl 結合 $\alpha \cdot C^2_{\rm EC}$ を通じた curvature 起源の cross-term が生じる。

数値例（ $R=3, L=\kappa=1$ ）:

| $\alpha$ | $\partial^2 V/\partial\varepsilon\partial s$ |
|---|---|
| $0$（LC） | $128\pi^4 \approx 1.247\times 10^4$ |
| $-1$（EC） | $11648\pi^4/27 \approx 4.202\times 10^4$ |

[スクリプト: `paper03ec/squash_shear_cross_term.py`]

### C.5 3 トポロジー cross-term 比較まとめ

| トポロジー | $\partial^2 V/\partial\varepsilon\partial s\|_{\varepsilon=s=0}$ | 体積保存 | $\alpha$ 依存 | 起源 |
|---|---|---|---|---|
| $T^3$ | $48\pi^4 L\eta^2 r/\kappa^2 \neq 0$ | ✗ | なし | kinematic（体積因子積） |
| $S^3$ | $0$（SymPy exact） | ✓ | なし | 消滅（ $s\to -s$ 対称） |
| $Nil^3$ | $(384\pi^4 LR^2 - 8192\pi^4 L\alpha\kappa^2)/(9R\kappa^2) \neq 0$ | ✓ | あり | curvature/Weyl（ $C^2_{01}$ 非対称） |

この表が定理 9 の証明の核心を構成する： $T^3$ は体積非保存による kinematic 起源、 $S^3$ は体積保存 + 対称性による消滅、 $Nil^3$ は体積保存にもかかわらず $C^2_{01}$ の方向非対称性による Weyl 起源の cross-term を持つ。
