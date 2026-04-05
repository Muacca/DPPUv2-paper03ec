## 3. Mode dictionary: $S^3\times S^1$ under EC-Weyl

本節では $S^3\times S^1$ の均質幾何自由度に対する EC-Weyl 結合の効果を、 $S^3$ 固有の mode dictionary として確定する。spin-0 セクターの非伝播性は §2.4 の **定理 2 (Palatini 保護 EC 版)** により既に確定しているので、ここではその $S^3$ での具体的帰結、すなわち **定理 4 (Palatini 普遍性)** とモード辞書への反映に集中する。

### 3.1 定理 4 (Palatini 普遍性)

**定理 4 (Palatini 普遍性)**: $S^3\times S^1$ において、AX 背景と VT 背景の spin-2/spin-1 質量は SymPy 厳密計算で完全に一致する：

**spin-2 質量（5 重縮退）**:
$$\mu^2_{\rm spin2}(S^3) = \frac{128\pi^2 Lr}{3\kappa^2} - \frac{16384\pi^2 L\alpha}{3r}$$

AX 背景での $\partial^2 V/\partial\varepsilon^2|_{\varepsilon=0}$ と VT 背景での同量は完全に等しい（VT $-$ AX $= 0$, SymPy exact）。また、この式は $\eta, V$ を含まない：torsion 背景の種類に非依存である。

LC 部分（ $\alpha=0$ ）の寄与は $128\pi^2 Lr/(3\kappa^2)$、EC 補正は $-16384\pi^2 L\alpha/(3r)$ であり、後者は $\eta, V$ に非依存でトポロジー（ $S^3$ 幾何）のみに依存する。

*注*: 本稿のログと数値例は $r=3, L=\kappa=1, \alpha=0$ を用いているため、上式は $\mu^2_{\rm spin2}=128\pi^2 \approx 1263.3$ と評価される。LC 版 paper03 の本文は $r=1$ ベースの数値表示を用いており、その場合の本式は $128\pi^2/3 \approx 421.2$ となる。paper03 の「 $\mu^2_q = 48\pi^2$ 」との差は、半径の取り方に加えて spin-2 基底/Hessian の規約差に由来する。

**spin-1 質量（3 重縮退）**:
$$m^2_{\omega_k}(S^3) = \frac{16\pi^2 L^3}{\kappa^2 r}, \quad k = 0, 1, 2.$$

AX 背景での $m^2_{\omega}$ と VT 背景での $m^2_{\omega}$ は完全に一致する（比率 $= 1.000000$）。

数値例（ $r=3, L=\kappa=1, \alpha=0$ ）:
- spin-2: $m^2 = 128\pi^2 \approx 1263.3$（VT は AX と同値）
- spin-1: $m^2 = 16\pi^2/3 \approx 52.6$（VT は AX と同値）

**物理的解釈**: Palatini 普遍性は、torsion 背景（ $\eta \neq 0$ の AX vs $V \neq 0$ の VT）の違いが spin-2/spin-1 の質量スペクトルに影響しないことを意味する。 $S^3$ 幾何（ $C^i_{jk} = 4\varepsilon_{ijk}/r$ ）から来る kinematic 部分が質量を決定し、torsion は potential の形のみを変形する（Palatini 保護の帰結）。

[スクリプト: `paper03ec/s3_vt_spin_masses.py`]

### 3.2 Mode dictionary table (Table 1: S³×S¹)

$S^3\times S^1$ の EC-Weyl 下でのモード辞書をまとめる。

**Table 1: $S^3\times S^1$ mode dictionary under EC-Weyl**

| 物理量 | AX ($V=0, \eta\neq 0$) | VT ($\eta=0, V\neq 0$) | MX ($V\neq 0, \eta\neq 0$) |
|---|---|---|---|
| $C^2_{\rm EC}$ | $= C^2_{\rm LC} = 0$ | $= C^2_{\rm LC} = 0$ | $16V^2\eta^2/(3r^2)$ |
| AX/VT dropout | ✓ | ✓ | ✗ |
| spin-0: $r$ | 力学的（伝播） | 力学的（伝播） | 力学的（伝播） |
| spin-0: $\eta, V$ | $\eta$ 非力学的, $V=0$（定理 2） | $\eta=0$, $V$ 非力学的（定理 2） | 両方非力学的（定理 2） |
| spin-2 質量 $\mu^2$ | $128\pi^2 Lr/(3\kappa^2) - 16384\pi^2 L\alpha/(3r)$（5 重縮退） | AX と同値（定理 4） | — |
| spin-1 質量 $m^2_\omega$ | $16\pi^2 L^3/(\kappa^2 r)$（3 重縮退） | AX と同値（定理 4） | — |

$S^3\times S^1$ の MX セクターについては、補助スクリプト `paper03ec/s3_mx_vacuum_search.py` で停留条件と Hessian を点検した。その結果、paper03-EC の辞書対象範囲（ $\|\eta\|\le 3$ ）では $\nabla_{(r,\eta,V)}V_{\rm eff}=0$ と Hessian 正定値を同時に満たす安定な full stationary vacuum は見つからない。 $\alpha>0$ 側に実数 MX 停留点そのものは現れ得るが、いずれも負固有値を 1 つ以上持つ鞍点である。固定 $V$ スライスに現れる有限半径の井戸は $r$ 方向の 1 次元極小にとどまり、full stationary vacuum を意味しない。

### 3.3 P 保護: AX = Plücker, VT = Pfaffian, MX = P ≠ 0

Pontryagin 密度の保護構造は paper03 から継承される。

**AX モード（Plücker 保護）**: $V = 0$ のとき、Pontryagin 密度 $P = \langle R, \star R\rangle$ は Plücker 関係式により恒等的にゼロとなる。これは AX torsion が Riemann テンソルの反対称部分に寄与しないことによる幾何的消去である。

**VT モード（Pfaffian 保護）**: $\eta = 0$ のとき、 $P$ は Pfaffian 型の恒等式により消滅する。VT torsion は Riemann テンソルの対称部分のみを修正し、Pontryagin 密度の反対称な組み合わせへの寄与がゼロとなる。

**MX モード**: $V \neq 0, \eta \neq 0$ のとき、等方点での Pontryagin 密度は

$$
P_0 = \frac{2V\eta(V^2 r^2 + 9\eta^2 - 36)}{9r^3}
$$

として非零になり得る（paper03 §3.2 から継承）。

**EC の効果**: Weyl coupling $\alpha$ は作用の coupling 定数であり、Pontryagin 密度 $P$ 自体の定義には入らない。したがって AX / VT / MX の Pontryagin 分類は paper03 と同一であり、EC 拡張で新たに変わるのは有効ポテンシャル $V_{\rm eff}$ 側である。

### 3.4 EC-Weyl 下での変化

**AX/VT 安定性の保持**: EC 補正 $\alpha < 0$ のもとで AX/VT 安定真空の構造は LC と変わらない。定理 1 の AX/VT dropout により $C^2_{\rm EC} = C^2_{\rm LC}$ であるため、Weyl coupling は有効ポテンシャルを変形しない。

**EFT レベルの EC 固有修正（α-cubic）**: LC からの新たな非線形結合として

$$
\frac{\partial^3 V_{\rm EC}}{\partial\alpha\,\partial\eta\,\partial V} = -\frac{128\pi^2 V\eta r}{3} \neq 0
$$

が出現する（ $S^3$ ）。これは $\alpha$ に関する homogeneous EFT の cubic coefficient を与える EC 固有の cubic チャンネルであり、LC では存在しない（[スクリプト: `paper03ec/ec_cubic_vertex.py`]）。この項は Table 1 の二次モード辞書には現れず、均質 EFT の cubic order にのみ現れる。ここで係数は Appendix A の $S^3\times S^1$ 体積規約 $\mathrm{Vol}(S^3\times S^1)=2\pi^2 Lr^3$ に対応している。一方、NY-cubic チャンネル $\partial^3 V/\partial\theta\partial\eta\partial V = 4\pi^2 r^2$ は EC で不変（Palatini 保護の帰結）。

これに対して auxiliary な formal MIXING sector では

$$
C_\delta(S^3) \;=\;
\left.
\frac{\partial^3 V_{\rm eff}}{\partial\delta_0\,\partial\delta_1\,\partial\delta_2}
\right|_{\delta=0}
= 0,
\qquad
\frac{\partial C_\delta(S^3)}{\partial\alpha} = 0
$$

が成立する（定理 3, App. D）。したがって EC 拡張で新たに活性化する homogeneous cubic 構造は physical な $\eta$ - $V$ channel に限られ、 $\delta$ -sector は null のまま保たれる。
