"""Score plot + loading plot : même PCA, protéines alignées, couleurs selon le design."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import deimos                                                     # noqa: E402

PARAMS = dict(volcano_use_padj=True, volcano_p_thresh=0.05, volcano_ratio_min=2.0,
              volcano_lfc_min=1.0)


def _data(conds=("A", "B", "C"), n=300, seed=0):
    rng = np.random.default_rng(seed)
    labels = [f"{c}{r}" for c in conds for r in range(3)]
    cond = [c for c in conds for _ in range(3)]
    X = rng.normal(20, 0.3, (n, len(labels)))
    X[:30, [i for i, c in enumerate(cond) if c == conds[-1]]] += 3      # up in the last one
    mat = pd.DataFrame(X, columns=labels)
    meta = pd.DataFrame({"name": [f"P{i}" for i in range(n)], "Genes": [f"G{i}" for i in range(n)]})
    design = pd.DataFrame({"label": labels, "condition": cond})
    return mat, meta, design


def test_loadings_match_the_score_plot_pca(tmp_path):
    mat, meta, design = _data()
    fit = deimos._pca_fit(mat)
    f, L, S = deimos.plot_pca_loadings(mat, meta, design, pd.DataFrame({"name": meta["name"]}),
                                    ["B_vs_A", "C_vs_A"], [(f"P{i}", "Cluster_1") for i in range(30)],
                                    PARAMS, str(tmp_path))
    assert Path(f).exists() and len(L) == len(mat)
    L = L.set_index("name").loc[meta["name"]]
    np.testing.assert_allclose(L["PC1_loading"].values, fit["loadings"][:, 0])
    # samples of C and the proteins up in C lie on the same side of PC1
    sc = fit["coords"][design["condition"].values == "C", 0].mean()
    assert np.sign(L["PC1_loading"].values[:30].mean()) == np.sign(sc)
    assert set(L["group"].values[:30]) == {"C"}               # peak condition
    assert (L["group"].values[30:] == "Not significant (ANOVA)").all()
    assert (L["ANOVA_cluster"].values[:30] == "Cluster_1").all()


def test_single_contrast_uses_the_volcano_rule(tmp_path):
    mat, meta, design = _data(conds=("A", "B"))
    d = np.zeros(len(meta)); d[:30] = 3.0
    q = np.ones(len(meta)); q[:30] = 1e-4
    res = pd.DataFrame({"name": meta["name"], "B_vs_A_diff": d, "B_vs_A_p.val": q,
                        "B_vs_A_p.adj": q})
    f, L, S = deimos.plot_pca_loadings(mat, meta, design, res, ["B_vs_A"], [], PARAMS, str(tmp_path))
    g = L.set_index("name").loc[meta["name"], "group"].values
    assert set(g[:30]) == {"Higher in B"}
    assert (g[30:] == "No change (volcano rule)").all()


def test_drivers_cover_the_four_directions(tmp_path):
    mat, meta, design = _data()
    f, L, S = deimos.plot_pca_loadings(mat, meta, design, pd.DataFrame({"name": meta["name"]}),
                                       ["B_vs_A", "C_vs_A"], [], {**PARAMS, "pca_n_drivers": 20},
                                       str(tmp_path))
    tags = ";".join(L["driver"])
    for d in ("PC1+", "PC1-", "PC2+", "PC2-"):
        assert tags.count(d) == 5
    top_pc1 = L.sort_values("PC1_loading", ascending=False).head(5)
    assert top_pc1["driver"].str.contains("PC1\\+").all()
    assert list(S.columns[:5]) == ["label", "condition", "PC1", "PC2", "PC3"]
    assert len(S) == mat.shape[1]
    fit = deimos._pca_fit(mat)
    np.testing.assert_allclose(S["PC1"].values, fit["coords"][:, 0])


def test_no_driver_when_disabled(tmp_path):
    mat, meta, design = _data()
    f, L, S = deimos.plot_pca_loadings(mat, meta, design, pd.DataFrame({"name": meta["name"]}),
                                       ["B_vs_A", "C_vs_A"], [], {**PARAMS, "pca_n_drivers": 0},
                                       str(tmp_path))
    assert (L["driver"] == "").all() and Path(f).exists()


def test_dashboard_reads_the_deimos_pca(tmp_path):
    """Le dashboard relit PCA_Scores / PCA_Loadings : mêmes coordonnées, protéines
    motrices regroupées par direction, gris sous-échantillonné mais motrices gardées."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "bd", Path(__file__).resolve().parents[1] / "build_dashboardv7.py")
    bd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bd)
    mat, meta, design = _data(n=600)
    f, L, S = deimos.plot_pca_loadings(mat, meta, design, pd.DataFrame({"name": meta["name"]}),
                                       ["B_vs_A", "C_vs_A"], [(f"P{i}", "Cluster_1") for i in range(30)],
                                       {**PARAMS, "pca_n_drivers": 8}, str(tmp_path))
    xl = tmp_path / "w.xlsx"
    with pd.ExcelWriter(xl) as w:
        L.to_excel(w, sheet_name="PCA_Loadings", index=False)
        S.to_excel(w, sheet_name="PCA_Scores", index=False)
    d = bd.build_pca_loadings(str(xl), max_grey=100)
    assert len(d["scores"]) == mat.shape[1]
    assert d["scores"][0]["x"] == round(float(S["PC1"].iloc[0]), 3)
    assert sorted(d["drivers"]) == ["PC1+", "PC1-", "PC2+", "PC2-"]
    assert all(len(v) == 2 for v in d["drivers"].values())
    shown = {p["n"] for p in d["proteins"]}
    assert {n for v in d["drivers"].values() for n in v} <= shown
    assert d["n_grey_shown"] <= 100 + 8 and d["n_grey"] == 570
    assert d["sides"]["PC1+"] or d["sides"]["PC1-"]
    empty = tmp_path / "none.xlsx"
    pd.DataFrame({"a": [1]}).to_excel(empty, index=False)
    assert bd.build_pca_loadings(str(empty)) is None


def test_driver_labels_do_not_overlap():
    """Étiquettes des protéines motrices : hors du nuage, sans chevauchement entre elles,
    toutes dans les axes, une par protéine, même quand les points sont tassés au bord."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(1)
    L = rng.normal(scale=0.01, size=(3000, 2))
    L[:6] = [[-0.04, 0.0], [-0.04, 0.0005], [-0.0401, -0.0003], [-0.0399, 0.0002],
             [-0.04, 0.0001], [-0.0398, -0.0001]]                      # PC1- tassées
    coords = rng.normal(size=(12, 2))
    conds = np.repeat(["A", "B", "C", "D"], 3)
    drv = deimos.pca_drivers(L, coords, conds, 24)
    drv["PC2+"]["idx"].append(drv["PC1-"]["idx"][0])                  # doublon entre directions
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(L[:, 0], L[:, 1], s=3)
    fig.tight_layout()
    labels = np.array([f"GENE_{i:05d}" for i in range(len(L))])
    anns = deimos._place_driver_labels(ax, L, drv, labels, {})
    n_unique = len({i for v in drv.values() for i in v["idx"]})
    assert len(anns) == n_unique
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    from matplotlib.text import Text                     # boîte du texte seul (sans le trait)
    boxes = [Text.get_window_extent(a, renderer=r) for a in anns]
    axb = ax.get_window_extent(renderer=r)
    for k, b in enumerate(boxes):
        assert axb.x0 - 1 <= b.x0 and b.x1 <= axb.x1 + 1 and axb.y0 - 1 <= b.y0 and b.y1 <= axb.y1 + 1
        for b2 in boxes[k + 1:]:
            assert not b.overlaps(b2)
    plt.close(fig)


def test_missing_gene_and_description_do_not_crash(tmp_path):
    """Gènes / descriptions vides (NaN) : avec pandas >= 3, astype(str) les laisse en
    NaN (float) ; le loading plot doit les traiter comme des chaînes vides (repli sur
    l'accession pour l'étiquette), y compris pour les protéines motrices."""
    mat, meta, design = _data()
    meta["Genes"] = meta["Genes"].astype(object)
    meta.loc[::2, "Genes"] = np.nan
    meta["First.Protein.Description"] = [np.nan if i % 3 else f"Protein {i}" for i in range(len(meta))]
    f, L, S = deimos.plot_pca_loadings(mat, meta, design, pd.DataFrame({"name": meta["name"]}),
                                       ["B_vs_A", "C_vs_A"], [(f"P{i}", "Cluster_1") for i in range(30)],
                                       {**PARAMS, "pca_n_drivers": 40}, str(tmp_path))
    assert Path(f).exists()
    assert not L["Genes"].isna().any() and not L["Description"].isna().any()
    no_gene = L.set_index("name").loc[["P0", "P2"], "Genes"]
    assert list(no_gene) == ["P0", "P2"]                       # repli sur l'accession
