# 付録D. 再現方法と改訂注記

## D.1 コードとデータ

公開補足資料には、現行のCartan式・射影、外積代数計算、厳密な記号結果、3図の保存データを収録する。実行にはPython 3.11以降と `script/requirements.txt` の依存ライブラリが必要であり、notebookにはさらにJupyter環境が必要である。再計算と保存データの描画を分離する。 `script` ディレクトリから次を実行する。

```bash
python -m pip install -r requirements.txt
python -B paper03ec.py verify
python -B paper03ec.py recompute --output ../recomputed
python -B paper03ec.py render --data ../data/paper03ec_v4 --output ../rendered
```

出力先には新しいフォルダを指定する。 `DPPUv2_Paper_Figures_EC.ipynb` は隣接するviewerから保存データを読み、再計算は明示的に選択した場合に行い、実行ごとに日時付きの新しいフォルダへ出力する。manifestには指定profile、座標、規約、ファイルハッシュを記録する。記号計算はテンソル恒等式、スカラー縮約、静的横断変形、Euler式、停留枝の係数、連成profileの周波数項、Nilの調和形式格子の式を再現し、図の描画はその計算とは独立に保存数値データを使う。

## D.2 改訂注記

v4では完全NYの曲率側係数を $`1/2`$ とし、大域的規格化をPaper01〜02に揃えた。以前の計算入力との対応は $`T_{\mathrm{geom}}=-\tau`$ であり、現行補足コードは幾何学的torsionを直接使う。既約EC-WeylとRicciトレース除去後のEC曲率を区別し、Paper02の静的恒等式を本稿で指定した変分へ展開した。

従来の普遍的なPalatini保護とspin-2運動構造の不変性、全一様変数に関するNil極小と $`\lvert\kappa^2\theta\rvert=1`$ の閾値、粒子質量およびspin-1 tachyonとしての解釈は撤回する。第3〜5節には、それらに対応する制限された方程式と変形係数を示した。内部双対縮約によるPontryagin解釈と、それに基づくトポロジカル・chiralな結論も撤回する。付録Cに座標と経路の比較を示し、相互作用の次数や規格化を含む旧版からの変更内容は、リポジトリの `CHANGELOG.md` に記載する。
