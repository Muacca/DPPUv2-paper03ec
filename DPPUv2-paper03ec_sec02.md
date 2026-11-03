# 2. 幾何・作用・変分条件

## 2.1 coframeとtorsion

向き付けたEuclidean正規直交frameを使う。空間添字は $`i,j,k=0,1,2`$ 、円方向は3とする。円座標の周期は1、 $`L`$ は物理円周長である。twistとmixingが零のcoframeを次で定める。

式(G-01)。

$$
e^i=r M^i{}_j\sigma^j,\qquad e^3=L\thinspace \mathrm{d}y,\qquad
\det M=1,\quad y\sim y+1,\quad r,L>0.
$$

球面では巡回的に $`\mathrm{d}\sigma^i=4\sigma^j\wedge\sigma^k`$ とし、round球面の物理半径は $`r/2`$ になる。トーラスの形式は閉じ、三空間角の周期を各 $`2\pi`$ とする。Nilでは $`\sigma^0=\mathrm{d}x`$ 、 $`\sigma^1=\mathrm{d}w`$ 、 $`\sigma^2=\mathrm{d}z-x\mathrm{d}w`$ を用い、周期を $`(2\pi,2\pi,4\pi^2)`$ 、 $`x`$ 面の貼合せを $`z\mapsto z+2\pi w`$ とする。変分中もこの同一視を固定する。

式(G-02)。

$$
\mathrm{Vol}=c_M Lr^3,\qquad
c_{S^3}=\frac{\pi^2}{4},\quad c_{T^3}=(2\pi)^3,\quad
c_{\mathrm{Nil}^3}=(2\pi)^4.
$$

torsionには第一Cartan式で定義する幾何学的torsionを使い、軸性成分とベクトルトレース成分を円方向へ整列させる。

式(G-03)。

$$
\begin{aligned}
T^a&=\mathrm{d}e^a+\omega^a{}_b\wedge e^b
=\frac12 T^a{}_{bc}e^b\wedge e^c,\qquad h=\frac{\eta}{r},\\
T^{\mathrm{AX}}_{ijk}&=2h\varepsilon_{ijk},\qquad
T^{\mathrm{VT}}_{abc}=\frac{V}{3}
(\delta_{ac}\delta_{b3}-\delta_{ab}\delta_{c3}).
\end{aligned}
$$

AXは $`V=0`$ 、VTは $`h=0`$ 、MXは両方を許す族とする。補足コードは幾何学的torsionを直接使う。接続の符号と添字順序を次で固定する。

式(G-04)。

$$
\begin{aligned}
\mathrm{d}e^a&=\frac12 C^a{}_{bc}e^b\wedge e^c,\qquad
\omega_{ab}=\Gamma_{abc}e^c,\\
\Gamma^{\mathrm{LC}}_{abc}&=\frac12(C_{abc}+C_{cba}-C_{bac}),\\
K_{abc}&=-\frac12(T_{abc}+T_{bca}-T_{cab}),\qquad
\Gamma=\Gamma^{\mathrm{LC}}+K.
\end{aligned}
$$

## 2.2 三つの曲率二乗項

既約Weylテンソルを $`W_{\mathrm{EC}}`$ 、Ricciトレース除去テンソルを $`\widehat R_{\mathrm{EC}}`$ と書く。付録Aでは、対交換対称性やtorsion-freeの代数Bianchi恒等式を仮定せずに射影を定義する。四次元Euclidean内積では次が成り立つ。

式(C-01)。

$$
\begin{aligned}
\widehat R_{\mathrm{EC}}&=W_{\mathrm{EC}}+R^{(2)}+R^{(3)},\\
I_{\mathrm{tf}}&=I_{\mathrm W}+I_2+I_3,\qquad
I_{\mathrm W}=\lVert W_{\mathrm{EC}}\rVert^2,\quad
I_a=\lVert R^{(a)}\rVert^2\ (a=2,3).
\end{aligned}
$$

成分数は $`10+9+1`$ である。 $`R^{(2)}`$ は対交換反対称かつRicci-free、 $`R^{(3)}`$ は4形式である。LC比較値は $`I_{\mathrm{LC}}=\lVert W_{\mathrm{LC}}\rVert^2`$ とする。次の作用で $`X\in\lbrace\mathrm{LC},\mathrm W,\mathrm{tf}\rbrace`$ を比較する。

式(C-02)。

$$
\begin{aligned}
S_X&=-\int\mathrm{vol}\left[\frac{R_{\mathrm{EC}}}{2\kappa^2}
+\alpha I_X+\theta\mathcal N\right],\qquad \kappa^2>0,\\
\mathcal N\thinspace \mathrm{vol}&=T^a\wedge T_a-e^a\wedge e^b\wedge R_{ab}
=\mathrm{d}(e^a\wedge T_a).
\end{aligned}
$$

負号はPaper03の縮約ポテンシャル規約に揃える。宇宙項やsourceは加えない。 $`\alpha`$ と $`\theta`$ は固定定数である。静的ポテンシャル $`U_X`$ は一定profileに対する作用値である。

## 2.3 固定量と変分量

| 族 | 固定量と変分量 |
|---|---|
| 円方向torsion | $`r,L,M`$ を固定しtwist・mixingを零とする。滑らかで周期的な $`\eta(y),V(y)`$ を変分する。 |
| 静的shape | $`r,L,\eta,V`$ を固定し、行列式1の指数経路で $`M`$ を変分する。 |
| 静的twist | 空間幾何とtorsion成分を固定し、 $`e^3=L(\mathrm{d}y+\omega_i\sigma^i)`$ とする。 |
| 静的mixing | twistを零とし、coframe回転後のframeに同じtorsion成分を再指定する。 |
| 連成周期shape–twist | 球面上で $`r,L,\eta,V`$ とtorsion軸を固定し、第6節の $`M=\mathrm{diag}(f,f^{-1},1)`$ 、 $`\omega=(w,0,0)`$ を変化させる。 |

静的mixingの処方は、metricの幾何に対するtorsion軸を変える。全場を同時に回すframeのgauge変換とは区別する。以下で半径を微分するときは $`\eta,V`$ を固定し、 $`h=\eta/r`$ は半径とともに変化する。半径臨界点を求める際に固定した $`L`$ についても、停留性を判定するため一次微分を示す。
