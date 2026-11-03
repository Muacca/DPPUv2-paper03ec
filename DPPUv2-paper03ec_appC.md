# 付録C. 座標写像・勾配・変分経路

非停留点の二階微分には座標写像の二階微分も寄与する。 $`q=q(x)`$ の場合の変換則は次である。

式(P-01)。

$$
\frac{\partial^2U}{\partial x^a\partial x^b}
=\frac{\partial q^i}{\partial x^a}H_{ij}\frac{\partial q^j}{\partial x^b}
+\frac{\partial U}{\partial q^i}\frac{\partial^2q^i}{\partial x^a\partial x^b}.
$$

## C.1 Nilの座標重複

Nilの対角パラメータ表示 $`M=\mathrm{diag}(1,f,f^{-1})`$ 、 $`f=(1+\varepsilon)^{2/3}(1+s)`$ を考える。二つのパラメータは同じ対角shapeを変化させる。正規化基底(H-01)では次となる。

式(P-02)。

$$
q_1=\frac{\sqrt6}{2}\log f,\qquad q_2=-\frac{1}{\sqrt2}\log f,
\qquad q_3=q_4=q_5=0,\qquad
\mathrm{rank}\frac{\partial(q_1,q_2)}{\partial(\varepsilon,s)}=1.
$$

対応する不変量は $`z=f^{-2}`$ である。 $`Az^2+Bz^4`$ を直接微分しても、勾配項を含めて式(P-01)を使っても次を得る。

式(P-03)。

$$
H_{\varepsilon,s}
=\begin{pmatrix}
\frac{88A+304B}{9}&\frac{32A+128B}{3}\\
\frac{32A+128B}{3}&20A+72B
\end{pmatrix}.
$$

勾配項によって、接写像のrankが1でもこの行列のrankは2になり得る。その固有値はこのパラメータ表示に属するものであり、独立な5つのshape方向は式(H-01)で与えられる。特に交差係数は次となる。

式(P-04)。

$$
\left.\partial_\varepsilon\partial_sU_X\right|_0
=\frac{128\pi^4Lr}{3\kappa^2}
-\frac{8192\pi^4L\alpha}{9r}.
$$

第一項はECスカラー曲率、第二項はこの族で既約EC-Weylに等しいLC-Weyl不変量から生じる。追加二つのRicci-free成分は、この整列族のshapeに依存しない。式(P-04)は、この経路の交差係数に含まれる二つの幾何学的寄与を分離している。

## C.2 球面・トーラスの経路辞書

球面の対角座標は $`M=\mathrm{diag}(a(1+s),a/(1+s),a^{-2})`$ 、 $`a=(1+\varepsilon)^{1/3}`$ であり、次の写像を持つ。

式(P-05)。

$$
q_1=\frac{\sqrt6}{3}\log(1+\varepsilon),\qquad
q_2=\sqrt2\log(1+s).
$$

round点ではshape勾配が零になる。式(H-02)の係数を $`H_0`$ とすると、対角座標での成分は $`H_{\varepsilon\varepsilon}=2H_0/3`$ 、 $`H_{ss}=2H_0`$ 、 $`H_{\varepsilon s}=0`$ となる。この写像と大域的体積因子により、各座標での正規化が定まる。

トーラスでは、 $`r,L,\eta,V`$ を固定し、体積を $`(1+\varepsilon)(1+s)`$ 倍する経路を比較する。円座標の周期を1とした規約では、この経路について次を得る。

式(P-06)。

$$
\begin{aligned}
\left.\partial_\varepsilon\partial_sU_X^{\mathrm{volume\ path}}\right|_0
&=c_{T^3}Lr^3\left[\frac{3\eta^2}{\kappa^2r^2}+\frac{V^2}{3\kappa^2}
-\chi_X\frac{16\alpha\eta^2V^2}{3r^2}\right],\\
\left.\partial_\varepsilon\partial_sU_X^{\mathrm{volume\ path}}\right|_{V=0,0}
&=\frac{24\pi^3Lr\eta^2}{\kappa^2}.
\end{aligned}
$$

体積保存の指数shape族では、対応する交差項は零となる。この違いは体積の変化による。トーラスのLC曲率は零でも、ECスカラー $`R_{\mathrm{EC}}=-6h^2-2V^2/3`$ と静的なトレース除去補正 $`16h^2V^2/3`$ は非零になり得るため、体積が変わる経路には寄与する。
