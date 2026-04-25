## 7. Physical Interpretation and Outlook

### 7.1 3 トポロジーの物理的階層

本論文の解析を通じて、 $S^3\times S^1$、 $T^3\times S^1$、 $Nil^3\times S^1$ の 3 トポロジーが EC-Weyl 結合に対して鮮明に異なる物理的役割を持つことが明らかになった。

**$T^3\times S^1$**: 平坦多様体（ $C^a_{bc} = 0$ ）であるため、spin-2/spin-1 の全モード質量がゼロとなり、EC-Weyl 結合は AX/VT/MX のいずれの背景においても新安定点を生まない（三重の triviality, §4.3）。 $T^3$ は理論の「trivial な基準点」として機能し、 $V_{\rm eff}\to 0$ の漸近枝が比較の基準を与える。

**$S^3\times S^1$**: 共形平坦多様体（ $C^2_{\rm LC}(\eta=V=0) = 0$ ）であるため、背景 Weyl テンソルの EC-Weyl 結合は $\eta=V=0$ 断面で不活性となる。しかし $S^3$ は非自明な幾何（ $C^i_{jk} = 4\varepsilon_{ijk}/r$ ）を持ち、spin-2 五重縮退 $\mu^2 = 128\pi^2 Lr/(3\kappa^2) - 16384\pi^2 L\alpha/(3r)$ と spin-1 三重縮退 $m^2_\omega = 16\pi^2 L^3/(\kappa^2 r)$ が定理 4（Palatini 普遍性）として確立される。また Pontryagin 構造（ $P \neq 0$ in MX, paper03 継承）を通じた Pontryagin 生成の基盤を提供する。

**$Nil^3\times S^1$**: 非共形平坦多様体（ $C^2_{\rm LC} = 4/(3R^4) \neq 0$ ）であるため、EC-Weyl 結合が $\alpha < 0$ で $\eta=V=0$ 断面の EC slice minimum を生成する唯一のトポロジーである。この slice minimum は LC（ $\alpha=0$ ）には見えず、EC 拡張により初めて相構造に現れる。さらに full homogeneous potential では $|\kappa^2\theta_{\rm NY}|<1$ の範囲で genuine local minimum となる。定理 6–8 が示す EC 固有の新物理（γ=1/2 スケーリング, quintet 分裂, 1 軸質量分裂）は、 $Nil^3$ の「中間的な非可換性」（ $T^3$ の全零と $S^3$ の三方向非零の中間）が、モード辞書の方向依存性として物理化した結果と理解できる。

この 3 段階の階層は「幾何的複雑さ」の単調増加とは異なる。 $T^3$（trivial） $\to$ $S^3$（LC が活躍, EC は共形保護） $\to$ $Nil^3$（EC 固有物理の舞台）という役割分担が、均質幾何の比較枠組みを与える。

### 7.2 EC 固有の新物理

本論文で確立した EC 固有（LC では不在）の物理的効果を総括する。

**① Nil³ EC slice minimum（定理 6）**: $\alpha < 0$ での $Nil^3$ の $\eta=V=0$ slice minimum は、LC 世界では単調ポテンシャルしか持たない理論に EC が新しい stationary branch を付け加えることを意味する。これは $|\kappa^2\theta_{\rm NY}|<1$ の benchmark 範囲で full homogeneous local minimum となる。 $r_0 \propto \sqrt{|\alpha|}$ のスケーリングは、Weyl coupling の強さが幾何的サイズを規定するという、一般相対論的には非自明な機構を示す。

**② quintet 分裂（定理 7）**: $Nil^3$ の 1 次元 Lie 代数構造（ $C^2_{01}$ のみ）が、EC slice-minimum branch 周りで spin-2 の五重縮退を破り $5\to 0+2+1+1$ に分裂させる。特に zero mode $q_3$ の存在は Nil³ 固有の幾何的対称性の「残存」を示し、LC 等方点とは本質的に異なる物理的内容を持つ。

**③ 1 軸質量分裂（定理 8）**: spin-1 triplet の中で（ $C^2_{01}$ の）非可換方向のみが質量を持つという 1 軸選択性は、 $Nil^3$ の非対称な内部構造が kinematic matrix に写し込まれた結果である。これは $Nil^3$ における方向依存的なモード辞書の最も明瞭な現れである。

**④ α-cubic 結合**: $S^3$ では $\partial^3 V/\partial\alpha\partial\eta\partial V = -128\pi^2 V\eta r/3$ という EC 固有の cubic coefficient が出現する（§3.4）。これは LC の NY-cubic チャンネル（ $= 4\pi^2 r^2$ , EC で不変）に加わる新たな非線形結合であり、homogeneous EFT における cubic channel を与える。一方、Pontryagin 密度 $P$ は $\alpha$ に依存しない（ $\partial P/\partial\alpha = 0$ ）ため、この α-cubic は Pontryagin 構造の外部に位置する結合として理解できる。

**⑤ δ-sector null result（定理 3）**: これと補完的に、formal auxiliary MIXING sector に対して定義される $C_\delta$ は $S^3, T^3, Nil^3$ の全てでゼロであり、 $\partial_\alpha C_\delta = 0$ である。したがって EC 拡張で新たに活性化する homogeneous cubic 構造は physical な $\eta$ - $V$ channel に限られ、auxiliary な $\delta$ -sector は null のまま保たれる。


### 7.3 What is not claimed here

本論は、以下を保留している。

- minisuperspace を超えた非一様モードを含む完全安定性
- full 4D deformed geometry が新しい Lie group manifold をなすという主張
- arbitrary non-invariant polynomial まで含めた完全な作用生成
- irreducible tensor torsion を含む fully general homogeneous torsion ansatz （本論は paper01–03 で採用した axial + vector-trace torsion ansatz に限定する）

これらは本論の結果と矛盾するものではなく、本論の言及範囲外として意識的に線引きしたものである。

### 7.4 Toward Future Extensions

本論文で確立した homogeneous three-topology bulk sector は、それ自体で閉じた比較枠組みを与えると同時に、各トポロジーで異なる mode dictionary と保護構造がどのように現れるかを明示した。以下の課題は、この比較枠組みを保ったまま、局所自由度・境界応答・物理的解釈へと拡張していく自然な次段階である。

1. 局所自由度を含む非一様ゆらぎの導入。
2. 境界面や interface を含む問題設定。
3. homogeneous bulk mode と localized mode の接続。
4. thermal / Lorentzian 解釈を含む物理的 observables への橋渡し。
