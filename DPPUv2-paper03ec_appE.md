# 付録E. 固定コンパクト商上のNil平坦経路

$`T_2,T_3`$ の経路は、次の正定値・単位行列式の水平計量に含まれる。座標周期と貼り合わせは第2節のものに固定する。shape変数 $`(M^{\mathsf T}M)_{22}`$ は1を保つので、式(H-03)は両経路に沿って一定である。

式(N-01)。

$$
\begin{aligned}
g_Q&=r^2\left[(\mathrm{d}x,\mathrm{d}w)Q(\mathrm{d}x,\mathrm{d}w)^{\mathsf T}
+(\mathrm{d}z-x\mathrm{d}w)^2\right]+L^2\mathrm{d}y^2,\\
Q&>0,\qquad\det Q=1,\qquad
G_{H^1}=\frac{\mathcal V}{4\pi^2r^2}Q^{-1}.
\end{aligned}
$$

最後の式は整数周期を持つ1形式 $`\mathrm{d}x/(2\pi),\mathrm{d}w/(2\pi)`$ のL² Gram行列であり、 $`\mathcal V=(2\pi)^4Lr^3`$ とする。これらは閉かつ余閉である。逆座標metricの水平blockは $`r^{-2}Q^{-1}`$ 、体積密度は一定なので、それぞれの発散が零となる。Nil latticeの可換化はrank 2であり、この二つの1形式と積空間の $`\mathrm{d}y`$ が調和1形式を張る。LC Ricci固有値は水平で $`-1/(2r^2)`$ 、中心で $`+1/(2r^2)`$ 、円で零となるため、この積metric間の等長写像は水平調和部分空間とその整数格子を保つ。固定r,Lでの正規化した最短ノルム二乗は等長不変量となる。

式(N-02)。

$$
\begin{aligned}
\sigma(Q)&=\min_{m\in\mathbb Z^2\setminus\lbrace0\rbrace}m^{\mathsf T}Q^{-1}m,
\qquad a=\sqrt2q,\\
Q_{T_2}&=\mathrm{diag}(e^a,e^{-a}),
&\sigma(Q_{T_2})&=e^{-|a|},\\
Q_{T_3}&=\begin{pmatrix}\cosh a&\sinh a\\
\sinh a&\cosh a\end{pmatrix},
&\sigma(Q_{T_3})&=\min(\cosh a,2e^{-|a|}).
\end{aligned}
$$

T2では小さい方の対角係数が全非零整数ベクトルの下限を与え、座標ベクトルで達成される。T3ではEuclidean長さ二乗1の整数ベクトルが $`\cosh a`$ を与える。それ以外は長さ二乗が2以上なのでノルムは $`2e^{-|a|}`$ 以上となり、適切な符号の $`(1,\pm1)`$ が達成する。これで式(N-02)を無限整数格子上で証明できる。

異なる $`|q|`$ はT2の全経路、およびT3の $`0\le\sqrt2|q|<\log(3)/2`$ の区間で区別される。この区間の非零点ではT2で $`\sigma<1`$ 、T3で $`\sigma>1`$ となり、両経路も区別できる。一方、符号反転はlatticeを保つ具体的写像によって同定される。正規化座標 $`(X,Y,Z)=(x/(2\pi),w/(2\pi),z/(4\pi^2))`$ では次である。

式(N-03)。

$$
(X,Y,Z)\mapsto(Y,-X,Z-XY)\quad(T_2),\qquad
(X,Y,Z)\mapsto(-X,Y,-Z)\quad(T_3).
$$

これらは整数Heisenberg群の積と $`\mathrm{d}Z-X\mathrm{d}Y`$ の二乗を保ち、 $`Q(q)`$ を $`Q(-q)`$ に写す。水平の単位行列式変換は局所Heisenberg自己同型を与えるが、latticeを保つとは限らない。この局所と大域の区別は[固定Heisenberg商の幾何](https://arxiv.org/html/2008.05328v2)とも整合する。ここでは調和形式の整数ノルムを等長性への障害として用いる。平坦経路は連続的に異なるコンパクト計量を含み、その上で式(H-05)の横断微分は非零のままである。
