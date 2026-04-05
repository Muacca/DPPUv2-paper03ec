# DPPUv2 計算エンジン v5 — スクリプトディレクトリ

⇒ [English](README.md)

**論文**: "Einstein-Cartan + Nieh-Yan Theory on Homogeneous 3-Manifolds: Unified Analysis of S³×S¹, T³×S¹, and Nil³×S¹"（paper03-EC）

Einstein-Cartan + Nieh-Yan フレームワークにおける3つの一様トポロジー（S³×S¹・T³×S¹・Nil³×S¹）を対象とした数値・記号計算のための Python パッケージ群と実行スクリプト。

---

## ディレクトリ構成

```
script/
├── docs/                      # 技術ドキュメントと規約
├── dppu/                      # メイン Python パッケージ（DPPUv2 Engine v5）
│   ├── geometry/              # 計量・体積形式・構造定数
│   ├── connection/            # Levi-Civita 接続・Contortion・EC 接続
│   ├── curvature/             # Riemann・Ricci・Hodge 双対・Pontryagin・Weyl
│   ├── torsion/               # トーションモード・Ansatz・Nieh-Yan 密度
│   ├── action/                # ラグランジアン・有効ポテンシャル・安定性分類
│   ├── topology/              # 統一エンジン（S³×S¹・T³×S¹・Nil³×S¹）
│   ├── engine/                # 計算パイプライン・ロギング・チェックポイント
│   ├── kk/                    # Kaluza-Klein モード抽出と検証
│   └── utils/                 # 共通ユーティリティ（Levi-Civita 記号・記号計算・可視化）
│
└── scripts/                   # 実行スクリプト
    ├── paper03ec/             # paper03-EC 固有の解析スクリプト
    ├── proofs/                # 解析的・記号的証明スクリプト
    └── visualize/             # インタラクティブビューアとノートブック
```

### `docs/` — ドキュメント

技術ドキュメントと規約：
- [DPPUv2 Engine CONVENTIONS](docs/CONVENTIONS_ja.md) — エンジンコアの規約と仕様
- [DPPUv2 SymPy guideline](docs/SymPy_guideline_ja.md) — SymPy 使用ガイドラインとベストプラクティス

---

## パッケージ概要（dppu/）

| モジュール | 役割 | 主要クラス・関数 |
|-----------|------|----------------|
| [`geometry`](dppu/geometry/README_ja.md) | 計量・フレーム場定義 | `build_metric`, `frame_field` |
| [`connection`](dppu/connection/README_ja.md) | EC 接続の構築 | `levi_civita`, `contortion`, `ec_connection` |
| [`curvature`](dppu/curvature/README_ja.md) | 曲率テンソル群・Pontryagin・Weyl | `RiemannTensor`, `compute_pontryagin_inner_product`, `WeylTensor` |
| [`torsion`](dppu/torsion/README_ja.md) | トーション構造 | `Mode`, `NyVariant`, `build_torsion_tensor` |
| [`action`](dppu/action/README_ja.md) | 作用・安定性解析 | `build_lagrangian`, `classify_stability` |
| [`topology`](dppu/topology/README_ja.md) | 全トポロジー・全 DOF 対応統一エンジン | `UnifiedEngine`, `DOFConfig`, `TopologyType`, `FiberMode` |
| [`engine`](dppu/engine/README_ja.md) | 15 ステップ計算パイプライン | `BaseFrameEngine`, `ComputationLogger`, `CheckpointManager` |
| [`kk`](dppu/kk/README_ja.md) | Kaluza-Klein モード抽出 | `extract_maxwell`, `extract_mass`, `extract_cs` |
| [`utils`](dppu/utils/README_ja.md) | 共通ユーティリティ | `epsilon_symbol`, `prove_zero`, `set_style` |

---

## 実行スクリプト概要（scripts/）

### paper03ec/ — paper03-EC 固有解析

#### モード辞書とプロパティ

| スクリプト | 説明 |
|-----------|------|
| `nil3_mode_dictionary.py` | Nil³ モードプロパティ：AX/VT ドロップアウト・spin-2 五重項分裂・spin-1 質量分裂・MX 背景 |
| `t3_mode_dictionary.py` | T³ モードプロパティ：平坦幾何・α 独立性・torsion-volume coupling |

#### MX 真空探索

| スクリプト | 説明 |
|-----------|------|
| `nil3_mx_vacuum_search.py` | Nil³×S¹ パラメータ空間での MX 真空探索 |
| `s3_mx_vacuum_search.py` | S³×S¹ パラメータ空間での MX 真空探索 |
| `t3_mx_vacuum_search.py` | T³×S¹ パラメータ空間での MX 真空探索 |

#### スペクトル解析と定理

| スクリプト | 説明 |
|-----------|------|
| `nil3_spin2_quintet_splitting.py` | Nil³ spin-2 五重項分裂（q₃ = ゼロモード）|
| `s3_vt_spin_masses.py` | Palatini 普遍性検証：AX/VT で spin-2/1 スペクトルが一致（V\*=0 保護）|
| `squash_shear_cross_term.py` | ε-s 交差項の分解（T³=kinematic / S³=0 / Nil³=curvature）|

#### EC 頂点と EFT

| スクリプト | 説明 |
|-----------|------|
| `ec_cubic_vertex.py` | EC 3点頂点：θ-cubic の保存と α-cubic（EC 固有）チャンネルの比較 |
| `nil3_false_vacuum_eft.py` | Nil³ 偽真空の EFT 構造：η=V=0 での α 誘起真空 |

---

### proofs/ — 解析的・記号的証明

| スクリプト | 定理・内容 |
|-----------|----------|
| `ax_vt_dropout_proof.py` | **定理 1（AX/VT ドロップアウト）**：C²_EC = C²_LC — 単成分トーションでは Weyl 二乗が不活性 |
| `c_delta_null_proof.py` | **定理 3（δ セクターヌル）**：全トポロジー・全モードで MIXING セクターの C_δ ≡ 0 |
| `palatini_protection_proof.py` | **Palatini 保護**：VT 背景でオンシェル V\* = 0（dV/dV = 0）|
| `t3_flatness_null_test.py` | **T³ 平坦性**：C²_LC ≡ 0 — 平坦トポロジーでは EC 偽真空なし |
| `gamma_scaling_proof.py` | **定理 6（γ=1/2 スケーリング）**：Nil³ 偽真空での r₀(α) ∝ √\|α\|・V_eff^c(α) ∝ √\|α\| |

---

### visualize/ — インタラクティブビューアとノートブック

| ファイル | 内容 |
|---------|------|
| `DPPUv2_Paper_Figures_EC.ipynb` | 論文図表生成用 Jupyter ノートブック |
| `DPPUv2_visualize_notebook_v4ec.ipynb` | 可視化ノートブック |
| `DPPUv2_interactive_viewer_v4ec.py` | インタラクティブ Jupyter ウィジェット：トポロジー・モード・NY バリアント選択・位相図・ポテンシャル V_eff(r)・軸スケール切替 |

---

## クイックスタート

### 依存関係のインストール

```bash
pip install sympy numpy scipy pandas tqdm matplotlib
```

### UnifiedEngine でエンジンを実行

```python
from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

# S³×S¹ squash・AX トーション・完全 Nieh-Yan：
cfg = DOFConfig(
    topology=TopologyType.S3,
    enable_squash=True,
    torsion_mode=Mode.AX,
    ny_variant=NyVariant.FULL,
)
engine = UnifiedEngine(cfg)
engine.run()

Veff = engine.data['potential']
fp   = engine.get_free_params()   # 有効な SymPy Symbol の dict
```

### レガシーエンジンのプリセット名を使用

```python
cfg = DOFConfig.from_engine('S3S1Engine')    # S³×S¹ squash
cfg = DOFConfig.from_engine('T3S1Engine')    # T³×S¹
cfg = DOFConfig.from_engine('Nil3S1Engine')  # Nil³×S¹ squash
cfg.torsion_mode = Mode.MX
engine = UnifiedEngine(cfg)
engine.run()
```

### 証明スクリプトの実行

```bash
# 定理 1：AX/VT ドロップアウト（C²_EC = C²_LC）
python scripts/proofs/ax_vt_dropout_proof.py

# 定理 3：δ セクターヌル（MIXING で C_δ ≡ 0）
python scripts/proofs/c_delta_null_proof.py

# Palatini 保護（VT 背景で V* = 0）
python scripts/proofs/palatini_protection_proof.py

# T³ 平坦性（C²_LC ≡ 0、EC 偽真空なし）
python scripts/proofs/t3_flatness_null_test.py

# 定理 6：Nil³ 偽真空 γ=1/2 スケーリング
python scripts/proofs/gamma_scaling_proof.py
```

### 解析スクリプトの実行

```bash
# モード辞書
python scripts/paper03ec/nil3_mode_dictionary.py
python scripts/paper03ec/t3_mode_dictionary.py

# スペクトル解析
python scripts/paper03ec/nil3_spin2_quintet_splitting.py
python scripts/paper03ec/s3_vt_spin_masses.py
python scripts/paper03ec/squash_shear_cross_term.py

# MX 真空探索（全トポロジー）
python scripts/paper03ec/nil3_mx_vacuum_search.py
python scripts/paper03ec/s3_mx_vacuum_search.py
python scripts/paper03ec/t3_mx_vacuum_search.py
```

---

## paper03-EC の主要結果と対応スクリプト

| 結果 | 内容（証明種別） | 対応スクリプト |
|-----|----------------|--------------|
| **定理 1：AX/VT ドロップアウト** | C²_EC = C²_LC — 単成分トーションでは Weyl 二乗が不活性 — **記号的証明** | `proofs/ax_vt_dropout_proof.py` |
| **定理 3：δ セクターヌル** | 全トポロジーで MIXING セクターの C_δ ≡ 0 — **記号的証明** | `proofs/c_delta_null_proof.py` |
| **Palatini 保護** | VT 背景でオンシェル V\* = 0 — **解析的** | `proofs/palatini_protection_proof.py`、`paper03ec/s3_vt_spin_masses.py` |
| **T³ 平坦性** | C²_LC ≡ 0 — 平坦トポロジーでは EC 偽真空なし — **数値 + 記号** | `proofs/t3_flatness_null_test.py`、`paper03ec/t3_mode_dictionary.py` |
| **Nil³ spin-2 五重項分裂** | q₃ = ゼロモード；4+1 分裂 — **数値** | `paper03ec/nil3_spin2_quintet_splitting.py`、`paper03ec/nil3_mode_dictionary.py` |
| **Nil³ 偽真空（EFT）** | α 誘起真空（η=V=0）；γ=1/2 スケーリング — **解析的 + 数値** | `paper03ec/nil3_false_vacuum_eft.py`、`proofs/gamma_scaling_proof.py` |
| **ε-s 交差項分解** | T³=kinematic / S³=0 / Nil³=curvature — **記号的** | `paper03ec/squash_shear_cross_term.py` |
| **EC 3点頂点** | θ-cubic 保存；α-cubic は EC 固有 — **記号的** | `paper03ec/ec_cubic_vertex.py` |
| **MX 真空存在** | 各トポロジーでの MX モード真空構造 — **数値スキャン** | `paper03ec/nil3_mx_vacuum_search.py`、`paper03ec/s3_mx_vacuum_search.py`、`paper03ec/t3_mx_vacuum_search.py` |

---

## トポロジー一覧

| トポロジー | Lie 群 | 背景曲率 | Koszul 公式 |
|-----------|-------|---------|-----------|
| S³×S¹ | SU(2) | +24/r² | 簡略（bi-invariant）|
| T³×S¹ | U(1)³ | 0（平坦）| 自明 |
| Nil³×S¹ | Heisenberg 群 | −1/(2R²) | 一般（bi-invariant 非成立）|

---

## チェックポイント機能

長時間計算の中断・再開が可能：

```bash
# チェックポイント付きで実行
python scripts/paper03ec/nil3_mx_vacuum_search.py \
    --output-dir output/ --checkpoint-dir ./checkpoints

# 中断後の再開（チェックポイントディレクトリを再指定するだけ）
python scripts/paper03ec/nil3_mx_vacuum_search.py \
    --output-dir output/ --checkpoint-dir ./checkpoints
```

---

## ライセンス

リポジトリルートの LICENSE ファイルを参照してください。

---

**Author**: Muacca
**Version**: DPPUv2 Engine v5
**Date**: 2026-03
