## Appendix D. auxiliary $\delta$-sector null theorem の詳細

本付録では、定理 3（ $\delta$-sector null theorem）の詳細を与える。目的は、EC 拡張で physical な $(\eta, V)$ channel に新しい $\alpha$-cubic が現れる一方で、auxiliary な MIXING 変数 $\delta_k$ に対しては新しい homogeneous cubic channel が現れないことを厳密に示すことである。

### D.1 定義

unified engine の formal MIXING sector に対して

$$
C_\delta(M^3;{\rm mode}) \;:=\;
\left.
\frac{\partial^3 V_{\rm eff}}{\partial\delta_0\,\partial\delta_1\,\partial\delta_2}
\right|_{\delta_0=\delta_1=\delta_2=0}
$$

を定義する。ここで $M^3 = S^3, T^3, Nil^3$、mode $=$ AX/VT/MX である。

本稿で確認する命題は

$$
C_\delta = 0,
\qquad
\frac{\partial C_\delta}{\partial\alpha} = 0
$$

の 2 つである。後者は前者の直ちの帰結だが、script では独立に出力して確認する。

### D.2 $S^3\times S^1$ の exact zero

$S^3$ では MIXING sector 自体は非自明に導入できる。したがって、もし auxiliary cubic channel が存在するなら、最も現れやすいのはこのトポロジーである。

SymPy による直接計算では、AX/VT/MX の各背景に対して

$$
\left.
\frac{\partial^3 V_{\rm eff}^{S^3}}{\partial\delta_0\,\partial\delta_1\,\partial\delta_2}
\right|_{\delta=0}
= 0
$$

が厳密に成立する。したがって

$$
\frac{\partial C_\delta(S^3)}{\partial\alpha}=0
$$

も厳密に成立する。

この結果は、§3.4 の

$$
\frac{\partial^3 V_{\rm eff}}{\partial\alpha\,\partial\eta\,\partial V}
= -\frac{128\pi^2 V\eta r}{3} \neq 0
$$

と対照的である。すなわち、EC 拡張で新たに活性化する cubic 構造は physical な $\eta$ - $V$ channel に限られ、auxiliary な $\delta$ -sector には現れない。

### D.3 $T^3\times S^1$ の exact zero

$T^3$ では平坦性 $C^a{}_{bc}=0$ が mode dictionary の大半を trivialize する。formal MIXING sector に対しても同様に、

$$
\left.
\frac{\partial^3 V_{\rm eff}^{T^3}}{\partial\delta_0\,\partial\delta_1\,\partial\delta_2}
\right|_{\delta=0}
= 0
$$

が AX/VT/MX の全背景で exact zero となる。したがって

$$
\frac{\partial C_\delta(T^3)}{\partial\alpha}=0.
$$

$T^3$ では physical sector にも auxiliary sector にも、EC 由来の新しい homogeneous cubic channel は現れない。

### D.4 $Nil^3\times S^1$ の exact zero

$Nil^3$ は background Weyl が非零であり、EC false vacuum と spin-2/spin-1 の方向依存分裂を生む唯一のトポロジーである。それにもかかわらず、formal MIXING sector では

$$
\left.
\frac{\partial^3 V_{\rm eff}^{Nil^3}}{\partial\delta_0\,\partial\delta_1\,\partial\delta_2}
\right|_{\delta=0}
= 0
$$

が AX/VT/MX の全背景で exact zero となる。よって

$$
\frac{\partial C_\delta(Nil^3)}{\partial\alpha}=0.
$$

このことは、 $Nil^3$ の EC 固有物理が background Weyl と $(r,\eta,V)$ sector に現れるのであって、新しい auxiliary $\delta$ -sector cubic channel に現れるのではないことを意味する。

### D.5 9 ケースのまとめ

| トポロジー | AX | VT | MX | $\partial_\alpha C_\delta$ |
|---|---|---|---|---|
| $S^3\times S^1$ | $0$ | $0$ | $0$ | $0$ |
| $T^3\times S^1$ | $0$ | $0$ | $0$ | $0$ |
| $Nil^3\times S^1$ | $0$ | $0$ | $0$ | $0$ |

以上より、formal MIXING sector に対して定義される auxiliary cubic coefficient は全 3 トポロジーで厳密に消え、EC 拡張の新しい homogeneous cubic 構造は physical な $\eta$ - $V$ channel に限られる。

[スクリプト: `proofs/c_delta_null_proof.py`]
