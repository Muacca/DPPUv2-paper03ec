# Einstein–Cartan minisuperspaceにおける曲率成分と制限付きtorsion変分

初版：v1、2026-04-06。過去の改訂：v2、2026-04-25；v3、2026-06-20。
現行版： **v4**、改訂日 **2026-11-02**。 [English](README.md) · [改訂履歴](CHANGELOG.md) · [記録DOI](https://doi.org/10.5281/zenodo.19425147)

ECスカラー曲率項を共通に保ち、LC-Weyl、既約EC-Weyl、Ricciトレース除去後のEC曲率を比較する。S³、T³、Nil³について円方向torsionの微分効果を導出し、指定した静的shape・twist・torsion軸の経路を解析する。円方向torsionの微分項は追加9成分と4形式部分が担う。独立した5つのshape方向からNilの勾配と厳密に平坦な二つの経路を示し、torsion軸の回転から追加9成分の向き依存性を明らかにする。 Ricciトレース除去作用の正の結合の周期停留枝と、有限角度の向きの極値を分類する。連成shape–twist族では、追加9成分の二次周波数項と三作用共通の正の曲率四次係数を分離する。Nilの平坦経路は固定コンパクト商の不変量でも調べる。

- [英語PDF](DPPUv2-paper03ec.pdf)
- 日本語本文：[導入](DPPUv2-paper03ec_sec01.md)、[幾何と作用](DPPUv2-paper03ec_sec02.md)、[円方向torsion](DPPUv2-paper03ec_sec03.md)、[静的変形](DPPUv2-paper03ec_sec04.md)、[向きと相互作用](DPPUv2-paper03ec_sec05.md)、[連成周期profile](DPPUv2-paper03ec_sec06.md)、[結論](DPPUv2-paper03ec_sec07.md)
- 付録：[射影](DPPUv2-paper03ec_appA.md)、[微分導出](DPPUv2-paper03ec_appB.md)、[座標写像と経路](DPPUv2-paper03ec_appC.md)、[再現方法と改訂注記](DPPUv2-paper03ec_appD.md)、[Nilの固定商](DPPUv2-paper03ec_appE.md)
- [補足コードの説明](script/README.md)、[保存データ](data/README.md)、[図のnotebook](script/DPPUv2_Paper_Figures_EC.ipynb)

[Paper02](https://github.com/Muacca/DPPUv2-paper02)の静的恒等式とポテンシャル分類は参照する。規約は改訂[Paper01](https://github.com/Muacca/DPPUv2-paper01)に揃え、スカラー作用は改訂[Paper03](https://github.com/Muacca/DPPUv2-paper03)の適用可能な制限へ接続する。結果は指定したEuclidean変分族についてのものであり、一般の連成profileや全理論の場の方程式と拘束への拡張は今後の課題である。

## 再現方法

Python 3.11以降を使い、 `script` ディレクトリで実行する。

```bash
python -m pip install -r requirements.txt
python -B paper03ec.py verify
python -B paper03ec.py recompute --output ../recomputed
python -B paper03ec.py render --data ../data/paper03ec_v4 --output ../rendered
```

`verify` はファイルを書かず厳密恒等式を導出する。 `recompute` は式と固定データを再計算し、 `render` は記号計算を実行せず保存データを描画する。出力先は新しいフォルダに限る。旧データと旧汎用engineをv4結果の再現に使わない。

## PDFの作成

```bash
cd LaTeX
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

出力は `LaTeX/main.pdf` となる。リポジトリ直下のPDFは閲覧用原稿である。リポジトリまたは記録DOIからファイルを選ぶ際は、原稿とデータの版情報を照合する。

著者：Muacca。コードライセンス：[MIT](LICENSE.txt)。
