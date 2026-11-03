# 付録A. 曲率射影とLC規格化

計量適合曲率は $`R_{abcd}=-R_{bacd}=-R_{abdc}`$ を満たす。 $`Q_{bd}=\sum_aR_{abad}`$ 、 $`R=\sum_bQ_{bb}`$ と定義する。対交換やRicciの対称性は仮定しない。生のRicci除去を次の線形写像で定める。

式(A-01)。

$$
\begin{aligned}
(\Pi R)_{abcd}={}&R_{abcd}
-\frac12(\delta_{ac}Q_{bd}-\delta_{ad}Q_{bc}
-\delta_{bc}Q_{ad}+\delta_{bd}Q_{ac})\\
&+\frac{R}{6}(\delta_{ac}\delta_{bd}-\delta_{ad}\delta_{bc}).
\end{aligned}
$$

対交換対称化を $`S`$ 、その4形式部分を $`F`$ とすると、必要な三射影は次である。

式(A-02)。

$$
\begin{aligned}
S_{abcd}&=\tfrac12(R_{abcd}+R_{cdab}),\qquad
F_{abcd}=\tfrac13(S_{abcd}+S_{acdb}+S_{adbc}),\\
W&=\Pi(S-F),\qquad R^{(3)}=F,\qquad
R^{(2)}=\Pi R-W-F.
\end{aligned}
$$

$`F`$ のトレースは零である。 $`S-F`$ は代数的Riemannテンソルの対交換対称性と巡回恒等式を持つため、そのRicci部分を除くと10成分のWeylが残る。反対称添字対の $`6\times6`$ 行列では、対交換反対称部分が15次元であり、 $`\Pi`$ が反対称Ricciの6成分を除くと9成分になる。4形式部分は1次元である。これらの互いに直交する部分空間が20次元のRicci-free像を構成する。全添字を足すEuclidean内積で射影は自己随伴かつ冪等となり、ノルムは式(C-01)のように加算される。補足の再現計算でも添字対行列の36基底に対する線形写像を構成する。

一様三次元幾何と平坦な円の積では、対角構造係数 $`n_i`$ に対する空間Ricci固有値は次となる。

式(A-03)。

$$
\begin{aligned}
\lambda_i&=\frac12\big[n_i^2-(n_j-n_k)^2\big],\qquad
R_{\mathrm{LC}}=\sum_i\lambda_i,\\
I_{\mathrm{LC}}&=2\sum_i\lambda_i^2-\frac23R_{\mathrm{LC}}^2.
\end{aligned}
$$

最後の等式は、三次元の恒等式 $`\lVert\mathrm{Riem}_3\rVert^2=4\lVert\mathrm{Ric}_3\rVert^2-R_3^2`$ を四次元Weylノルムへ代入すると従う。round球面では $`R_{\mathrm{LC}}=24/r^2`$ 、 $`I_{\mathrm{LC}}=0`$ 、トーラスでは両方零、未変形Nilでは $`R_{\mathrm{LC}}=-1/(2r^2)`$ 、 $`I_{\mathrm{LC}}=4/(3r^4)`$ となる。各表はこの規格化を使う。
