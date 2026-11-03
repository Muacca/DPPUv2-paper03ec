# 5. torsion軸の向きと相互作用の次数

## 5.1 異方性に応答する追加9成分

積metricと $`n_i`$ を固定し、twistを零とする。回転前の四次元frameで一定の単位ベクトル $`u`$ を定める。一定coframe回転後に整列torsionを再指定する処方は、回転前frameで次を使うことと等価である。

式(O-01)。

$$
T^{\mathrm{AX}}_{abc}=2h\varepsilon_{abcd}u_d,\qquad
T^{\mathrm{VT}}_{abc}=\frac{V}{3}(\delta_{ac}u_b-\delta_{ab}u_c),\qquad
u_au_a=1.
$$

静的な厳密式は次となる。

式(O-02)。

$$
\begin{aligned}
W_{\mathrm{EC}}&=W_{\mathrm{LC}},\qquad
R_{\mathrm{EC}}=R_{\mathrm{LC}}-6h^2-\frac23V^2,\\
I_2&=\frac83h^2V^2+4h^2\sum_{i=0}^2(n_j-n_k)^2u_i^2,\\
I_3&=\frac83h^2V^2,\qquad \mathcal N=\mathcal P=0.
\end{aligned}
$$

和の $`(i,j,k)`$ は巡回添字とする。この向き依存を生むのは追加9成分だけである。round球面と平坦トーラスでは零だが、異方的shapeやNilでは一般に零でない。スカラー、LC-Weyl、既約EC-Weylの各作用はこの一定の向きに依存しない。式(O-02)は一定の向き $`u`$ についての恒等式である。

順序付きmixing回転 $`O=R_{03}(\delta_0)R_{13}(\delta_1)R_{23}(\delta_2)`$ 、 $`\tan\vartheta_i=\delta_i`$ を使い、 $`c_i=(1+\delta_i^2)^{-1/2}`$ 、 $`s_i=\delta_i c_i`$ と置く。torsion軸は $`u=O^{\mathsf T}e_3`$ であり、成分は次になる。

式(O-03)。

$$
u=(-s_0,-c_0s_1,-c_0c_1s_2,c_0c_1c_2).
$$

各 $`u_i^2`$ はすべての $`\delta_j`$ について偶関数である。したがって、この族の $`\delta=0`$ では純mixingの三次単項式がすべて消える。トレース除去モデルの向きに関する二次係数は次となる。

式(O-04)。

$$
H^{\delta}_{\mathrm{tf}}
=-8\alpha c_M Lr^3 h^2\thinspace 
\mathrm{diag}\big((n_1-n_2)^2,(n_0-n_2)^2,(n_0-n_1)^2\big).
$$

未変形Nilでは $`\delta_0,\delta_1`$ の係数が等しく、 $`\delta_2`$ の係数が零になる。非零の $`\bar h`$ の周りで $`h=\bar h+\varphi`$ と展開すると、次の混合三次項が現れる。

式(O-05)。

$$
\left.\Delta U_{\mathrm{tf}}\right|_{\varphi\delta^2,\mathrm{Nil}^3}
=-\frac{8\alpha c_{\mathrm{Nil}^3}Lr^3\bar h}{r^2}
\thinspace \varphi(\delta_0^2+\delta_1^2).
$$

この相互作用は、軸性振幅の揺らぎと二つの向きの揺らぎを結び、追加の対交換反対称曲率成分だけから生じる。純粋な向きの三次項は偶奇性によって消える一方、非零の軸性背景振幅はこの混合三次項を許す。

## 5.2 整列混合項の場の次数

Paper02の整列した静的項 $`-16\alpha h^2V^2/3`$ は、torsion零では場の振幅について四次となる。その展開を、上で得た向きの相互作用と比較する。 $`h=\bar h+\varphi`$ 、 $`V=\bar V+v`$ とすると、 $`h^2V^2`$ の二次・三次・四次部分はそれぞれ次となる。

式(O-06)。

$$
\begin{aligned}
[h^2V^2]_2&=\bar V^2\varphi^2+4\bar h\bar V\varphi v+\bar h^2v^2,\\
[h^2V^2]_3&=2\bar V\varphi^2v+2\bar h\varphi v^2,\\
[h^2V^2]_4&=\varphi^2v^2.
\end{aligned}
$$

非停留背景では定数項と一次項も展開に含まれる。式(O-06)は $`\alpha`$ を固定し、揺らぎ $`\varphi,v`$ を用いてPaper02の静的ポテンシャルの場の次数を表したものである。三次項は選んだ背景に依存し、場の次数と微分次数は相互作用の異なる性質を表す。

## 5.3 有限角度の極値と縮退

式(O-02)は、幾何と振幅を固定した向きの大域的極値も決定する。単位3球面上で次のように書ける。

式(O-07)。

$$
\begin{aligned}
\Delta U_{\mathrm{tf}}&=K_u u^{\mathsf T}D u,\qquad K_u=-4\alpha\mathcal V h^2,\\
D&=\mathrm{diag}\big((n_1-n_2)^2,(n_0-n_2)^2,(n_0-n_1)^2,0\big),
\qquad\mathcal V=c_M Lr^3.
\end{aligned}
$$

$`K_u\ne0`$ のときLagrange未定乗数法から $`Du=\lambda u`$ を得る。全停留集合はDの固有空間 $`E_\lambda`$ 内の単位球面であり、接空間のHessianは次となる。

式(O-08)。

$$
u\in S(E_\lambda),\qquad
H_{\mathrm{tan}}=\left.2K_u(D-\lambda I)\right|_{u^\perp}.
$$

負のαでは $`S(\ker D)`$ が最小、最大固有空間の単位球面が最大となり、正のαでは逆転する。中間固有空間の単位球面はsaddleである。縮退は次で尽くされる。相異なる3係数を順に並べ、その隣接差を $`a,b>0`$ とすると、固有値は $`0,a^2,b^2,(a+b)^2`$ となる。不等間隔では4組の対蹠点、等間隔では中間固有値にsaddleの停留円が現れる。ちょうど2係数が等しい場合は固有値が $`0,0,d^2,d^2`$ となり、2つの停留円を得る。全係数が等しい場合、または $`\alpha=0`$ 、 $`h=0`$ では全向き球面で一定となる。

Nilでは $`D=r^{-2}\mathrm{diag}(1,1,0,0)`$ である。負のαの最小円は $`\mathrm{span}(e_2,e_3)`$ 、最大円は $`\mathrm{span}(e_0,e_1)`$ にあり、各々が接線の零方向1つと等しい横断固有値2つを持つ。非正のαでは全shapeにおいて向き依存の追加項が非負で、整列軸 $`u=e_3`$ で零を達成する。したがって同じshape領域・固定制御変数で一定軸の自由度を加えても、Paper02の整列した静的infimumは下がらない。これは式(O-02)の系である。
