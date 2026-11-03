# 付録B. 微分計算と成分別の整理

巡回的な空間構造係数と整列torsionに対して、第一・第二添字の反対称性とともに必要な非零接続係数は次である。

式(B-01)。

$$
\Gamma_{ijk}=\left(\frac{n_i+n_j-n_k}{2}-h\right)\varepsilon_{ijk},
\qquad \Gamma_{i3j}=-\frac V3\delta_{ij},\qquad
\Gamma_{3ij}=\frac V3\delta_{ij}.
$$

非座標frameにおける曲率の完全な式は次となる。

式(B-02)。

$$
\begin{aligned}
R_{abcd}={}&D_c\Gamma_{abd}-D_d\Gamma_{abc}
+\Gamma_{aec}\Gamma_{ebd}-\Gamma_{aed}\Gamma_{ebc}
+\Gamma_{abe}C^e{}_{cd}.
\end{aligned}
$$

最初の二項は、contortionを含む接続の微分を表す。円方向profileについて式(B-02)を射影すると、LC-Weylテンソルと、次の追加成分が得られる。

式(B-03)。

$$
\begin{aligned}
a_+&=\tfrac12\left(\dot h+\tfrac23hV\right),\qquad
a_-=\tfrac12\left(\dot h-\tfrac23hV\right),\\
R^{(2)}_{ijk3}&=a_+\varepsilon_{ijk},\qquad
R^{(2)}_{k3ij}=-a_+\varepsilon_{ijk},\\
R^{(3)}_{abcd}&=a_-\varepsilon_{abcd}.
\end{aligned}
$$

残りの成分は各添字対の反対称性から定まるか零となる。各成分に現れる24個の符号付き非零要素を二乗して足すと $`24a_+^2`$ と $`24a_-^2`$ となり、式(D-01)を得る。射影前曲率の縮約から式(D-01)のスカラー項、4形式部分の縮約から式(D-06)の $`12a_-=6\dot h-4hV`$ を得る。独立な外積代数の構成でも、第一Cartan式、二つの微分Bianchi恒等式、NY恒等式を照合する。

式(D-06)のPontryagin密度は、曲率2形式の外積からも接続のChern–Simons形式からも得られる。行列トレース規約では $`\mathrm{Tr}(R\wedge R)=-\sum_{a,b}R_{ab}\wedge R_{ab}`$ であるため、 $`\mathrm{d}\mathrm{CS}(\omega)=-\mathcal P\,\mathrm{vol}`$ となり、表示した原始関数の符号も固定される。

二成分間の打消しを明示するため、成分を識別する記号としてだけ別の重み $`b_2,b_3`$ を導入すると、対応するbulk項は次のように書ける。

式(B-04)。

$$
-b_2I_2-b_3I_3
=-6(b_2+b_3)\dot h^2-8(b_2-b_3)hV\dot h
-\frac83(b_2+b_3)h^2V^2.
$$

本稿で比較する作用では $`b_2=b_3=\alpha`$ であり、独立結合の物理モデルは採用していない。また $`I_2-I_3=16hV\dot h`$ は、 $`V`$ が変化する場合には積分しても零とは限らない。物理円座標 $`z`$ では次となる。

式(B-05)。

$$
\int_0^L(I_2-I_3)\thinspace \mathrm{d}z
=-8\int_0^L\dot V h^2\thinspace \mathrm{d}z.
$$

したがって式(D-02)の打消しは等しい重みの結果であり、二つの成分汎関数が一般に同一という意味ではない。 $`\alpha=0`$ でスカラー作用に戻ることと、 $`\dot h=\dot V=0`$ でPaper02の静的結果に戻ることは別々の照合である。
