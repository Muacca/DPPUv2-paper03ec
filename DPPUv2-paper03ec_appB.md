## Appendix B. EFT 係数完全リスト

本付録では、3 トポロジーそれぞれの EC-Weyl 下での EFT 係数（spin-0, spin-2, spin-1 の質量スペクトルおよびポテンシャル係数）をまとめる。

### B.1 S³×S¹ — AX EFT（等方点 $r = r_0$, paper03 App H 参照 + EC 補正）

**設定**: $S^3\times S^1$, AX 背景（ $V=0, \eta \neq 0$ ）, 等方点（ $\varepsilon = s = 0$ ）。

**spin-0 スペクトル**:

| DOF | 質量 $m^2$（解析式）| 数値（ $r=3, L=\kappa=1, \alpha=0$ ）|
|---|---|---|
| $r$（径方向） | $\partial^2 V/\partial r^2\|_{\rm on-shell}$ | paper03 §4.3 参照 |
| $\eta$（AX torsion） | 非伝播（ $G_{\eta\eta}=0$ ） | — |

**spin-2 スペクトル（5 重縮退）**:

$$
\mu^2_{\rm spin2}(S^3, {\rm AX}) = \frac{128\pi^2 Lr}{3\kappa^2} + \frac{-16384\pi^2 L\alpha}{3r}
$$

- LC 部分（ $\alpha=0$ ）: $128\pi^2 Lr/(3\kappa^2)$
- EC 補正: $-16384\pi^2 L\alpha/(3r)$（ $\alpha$ 依存, SymPy exact）
- 数値（ $r=3, L=\kappa=1, \alpha=0$ ）: $128\pi^2 \approx 1263.3$（5 重縮退）
- 注: 本稿では $r=3$ を用いている。LC 版 paper03 の本文では $r=1$ ベースの数値表示を採っており、その場合の本式は $128\pi^2/3 \approx 421.2$ となる。paper03 の $48\pi^2$ との差は、半径の取り方に加えて spin-2 基底/Hessian の規約差による。
- VT 背景での値は AX と完全に一致（定理 4）

**spin-1 スペクトル（3 重縮退）**:

$$
m^2_{\omega_k}(S^3, {\rm AX}) = \frac{16\pi^2 L^3}{\kappa^2 r}, \quad k = 0, 1, 2.
$$

- 数値（ $r=3, L=\kappa=1$）: $16\pi^2/3 \approx 52.6$（3 重縮退）
- $\eta, V$ に非依存（Palatini 普遍性, 定理 4）; VT でも同値

---

### B.2 Nil³×S¹ — EC false vacuum EFT

**設定**: $Nil^3\times S^1$, EC false vacuum（ $\alpha = -a^2 < 0, \eta = V = 0$, $r_0 = 4\kappa a/\sqrt{3}$ ）。

ここで false vacuum は、homogeneous effective potential 上の局所極小を指す。

**False vacuum の同定**:

$$
r_0 = \frac{4\kappa}{\sqrt{3}}\sqrt{|\alpha|} = \frac{4\kappa a}{\sqrt{3}}, \qquad V_0 = V_{\rm eff}(r_0) = \frac{32\sqrt{3}\,\pi^4 La}{3\kappa} > 0.
$$

**全モード二次係数スペクトル**（ $\delta\Phi = \Phi - \Phi^*$ 周りの二次作用）:

| DOF | 展開変数 | 二次係数（解析式） | 数値（ $a=\kappa=L=1$ ） | $\alpha$ 依存性 |
|---|---|---|---|---|
| $r$ | $\delta r$ | $2\sqrt{3}\pi^4 L/(a\kappa^3)$ | 337.4 | $\propto 1/a$ |
| $\eta$（非伝播） | $\delta\eta$ | $128\sqrt{3}\pi^4 La/\kappa$ | 21,596 | $\propto a$（ $\eta$ 方向曲率） |
| spin-2 squash $\varepsilon$ | $\delta\varepsilon$ | $6272\sqrt{3}\pi^4 La/(27\kappa)$ | 39,192 | $\propto a$ |
| spin-1 $\omega_2$ | $\delta\omega_2$ | $6\sqrt{3}\pi^4 L^3/(a\kappa^3)$ | 1,012 | $\propto 1/a$ |
| spin-1 $\omega_0, \omega_1$ | $\delta\omega_{0,1}$ | $0$ | 0 | — |
| spin-2 zero mode $q_3$ | $\delta q_3$ | $0$ | 0 | — |
| spin-2 $q_4, q_5$ | $\delta q_{4,5}$ | $m^2(q_4) = m^2(q_5)$ | 10,798 | — |

ここで $\eta$ 行は propagating mass ではなく、Palatini 保護（ $G_{\eta\eta}=0$ ）の下での $\eta$ 方向ポテンシャル曲率を表す。

**二次有効作用（EFT）**:

$$
{V_{\rm eff}}^{Nil^3}_{\rm EFT}(\delta\Phi) \approx V_0 + 
\frac{1}{2}\frac{2\sqrt{3}\pi^4 L}{a\kappa^3}(\delta r)^2 + 
\frac{1}{2}\frac{128\sqrt{3}\pi^4 La}{\kappa}(\delta\eta)^2 + 
\frac{1}{2}\frac{6272\sqrt{3}\pi^4 La}{27\kappa}(\delta\varepsilon)^2 + 
\frac{1}{2}\frac{6\sqrt{3}\pi^4 L^3}{a\kappa^3}(\delta\omega_2)^2 + 
O(\delta^3).
$$

（ $\delta r$ と $\delta\eta$ の交差項はゼロ: $\partial^2 V/\partial R\partial\eta|_{r_0,0} = 0$）

[スクリプト: `paper03ec/nil3_false_vacuum_eft.py`, `paper03ec/nil3_spin2_quintet_splitting.py`]

---

### B.3 T³×S¹ — Trivial EFT

**設定**: $T^3\times S^1$, $C^a_{bc} = 0$。

**有効ポテンシャルの閉じた式**（AX 背景, $V=0$ ）:

$$
V_{\rm eff}(T^3,\, r,\, \eta;\, \alpha,\, \theta)\big|_{V=0} = \frac{48\pi^4 Lr}{\kappa^2}\eta^2.
$$

MX 背景での全式（ $V\neq 0, \eta\neq 0$ ）:

$$
V_{\rm eff}(T^3,\, r,\, \eta,\, V;\, \alpha,\, \theta) = \frac{16\pi^4 Lr}{3}\left[\frac{V^2 r^2 + 9\eta^2}{\kappa^2} + 6V\eta r\theta - 16\alpha V^2\eta^2\right].
$$

（注： $R_{\rm EC}/(2\kappa^2)$ 起源の項（ $V^2 r^2$, $9\eta^2$ ）のみ $1/\kappa^2$ を持ち、 $\theta\cdot N_{\rm NY}$ 起源の項と $\alpha\cdot C^2_{\rm EC}$ 起源の項は $\kappa^2$ に非依存である。 $\kappa = L = 1$ では sec. 4.4 の式と一致する。）

**全モード質量**:

| DOF | 質量 $m^2$ | 備考 |
|---|---|---|
| $r$ | on-shell 依存 | LC 部分のみが寄与 |
| $\eta, V$ | $G_{\eta\eta} = G_{VV} = 0$ | 非力学的（Palatini 保護） |
| spin-2（squash, shear） | $= 0$ | $C^a_{bc}=0$ による平坦性 |
| spin-1（ $\omega_k$ ） | $= 0$ | 質量ゼロ |
| ε-s cross-term | $48\pi^4 L\eta^2 r/\kappa^2$ | kinematic 起源（定理 5） |

[スクリプト: `paper03ec/t3_mode_dictionary.py`, `paper03ec/t3_mx_vacuum_search.py`]
