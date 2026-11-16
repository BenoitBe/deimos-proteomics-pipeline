"""Les volcanos vont dans un sous-dossier de out_dir (défaut "Volcanos")."""
import os
import re

import numpy as np
import pandas as pd
import pytest

import deimos
from deimos import VOLCANO_DEFAULT_SUBDIR, volcano_dir
from limma_ebayes import _bh_correction

PARAMS = dict(volcano_use_padj=True, volcano_p_thresh=0.05, volcano_ratio_min=2.0,
              volcano_lfc_min=1.0)


def _df(n=400, n_sig=20, seed=0):
    rng = np.random.default_rng(seed)
    p = rng.uniform(size=n)
    p[:n_sig] = rng.uniform(0, 1e-4, n_sig)
    lfc = rng.normal(0, 0.5, n)
    lfc[:n_sig] = rng.choice([-2, 2], n_sig) + rng.normal(0, 0.2, n_sig)
    out = {"name": [f"P{i}" for i in range(n)]}
    for c in ("A_vs_B", "C_vs_D"):
        out[f"{c}_diff"] = lfc
        out[f"{c}_p.val"] = p
        out[f"{c}_p.adj"] = _bh_correction(p)
    return pd.DataFrame(out)


def test_the_default_is_a_volcanos_folder(tmp_path):
    out = str(tmp_path)
    d = volcano_dir({}, out)
    assert os.path.basename(d) == "Volcanos"
    assert os.path.dirname(os.path.abspath(d)) == os.path.abspath(out)
    assert os.path.isdir(d), "le dossier doit être créé, pas seulement nommé"
    assert VOLCANO_DEFAULT_SUBDIR == "Volcanos"


def test_an_empty_value_keeps_the_old_behaviour(tmp_path):
    """Comportement d'avant le 12/11/2026 : volcanos à la racine de out_dir."""
    out = str(tmp_path)
    for value in ("", False, None, "   ", "..."):
        d = volcano_dir({"volcano_subdir": value}, out)
        assert os.path.abspath(d) == os.path.abspath(out), value
    assert not os.path.isdir(os.path.join(out, "Volcanos"))


def test_a_custom_name_is_honoured(tmp_path):
    d = volcano_dir({"volcano_subdir": "figures_volcan"}, str(tmp_path))
    assert os.path.basename(d) == "figures_volcan"
    assert os.path.isdir(d)


def test_a_path_cannot_escape_out_dir(tmp_path):
    """Même garde-fou que dashboard_filename : une valeur de config ne doit pas
    pouvoir faire écrire ailleurs que dans out_dir."""
    out = tmp_path / "sortie"
    out.mkdir()
    for value in ("../../ailleurs", "/etc/ailleurs", "C:\\\\Windows\\\\ailleurs",
                  "a/b/ailleurs", "..\\\\..\\\\ailleurs"):
        d = volcano_dir({"volcano_subdir": value}, str(out))
        assert os.path.basename(d) == "ailleurs", value
        assert os.path.dirname(os.path.abspath(d)) == os.path.abspath(str(out)), value
        assert ".." not in d


def test_create_false_does_not_touch_the_disk(tmp_path):
    """Le nettoyage interroge le dossier sans vouloir le créer."""
    d = volcano_dir({}, str(tmp_path), create=False)
    assert os.path.basename(d) == "Volcanos"
    assert not os.path.exists(d)


def test_plot_volcanoes_writes_into_that_folder():
    """Les deux savefig (individuels + facet) doivent viser volc_dir, pas out_dir :
    c'est exactement ce qui était faux avant le correctif."""
    import inspect

    import deimos
    src = inspect.getsource(deimos.plot_volcanoes)
    assert "volc_dir = volcano_dir(params, out_dir)" in src
    joins = re.findall(r"os\.path\.join\((\w+), *[\"f]", src)
    assert joins, "aucun os.path.join trouvé dans plot_volcanoes"
    assert set(joins) == {"volc_dir"}, (
        f"des PNG sont encore écrits ailleurs que dans volc_dir : {set(joins)}")


def test_the_cleanup_sweeps_the_resized_copies_in_the_subfolder():
    """_resize_png écrit la copie _xl À CÔTÉ de l'original, donc dans le
    sous-dossier : la boucle de nettoyage, qui ne descend pas, doit le balayer
    explicitement, sinon le dossier Volcanos se remplit de copies."""
    src = open("deimos.py", encoding="utf-8").read()
    i = src.index("[CLEAN] Removing temporary figures")
    block = src[i:i + 2200]
    assert "volcano_dir(params, out_dir, create=False)" in block
    assert "os.listdir(volc_dir)" in block
    # et il ne doit supprimer QUE les copies redimensionnées
    assert '"_xl" in fn' in block


def test_the_facet_is_no_longer_deleted():
    """Le facet était supprimé par le nettoyage parce qu'il vivait à la racine
    et ne commençait pas par volc_. Dans le sous-dossier, il survit."""
    import inspect

    import deimos
    src = inspect.getsource(deimos.plot_volcanoes)
    assert 'os.path.join(volc_dir, "all_volcano_plots_grid.png")' in src


def test_the_config_default_matches_the_code_default():
    from config import _DEFAULTS
    assert _DEFAULTS["volcano_subdir"] == VOLCANO_DEFAULT_SUBDIR


# --- Vérifications sur les fichiers réellement produits -----------------------
# Les tests ci-dessus lisent le source. Ceux-ci exécutent : c'est un test de
# motif qui avait laissé passer un dashboard écrit sous l'ancien nom.

def test_a_real_run_puts_every_png_in_the_subfolder(tmp_path):
    out = str(tmp_path)
    files, facet = deimos.plot_volcanoes(_df(), ["A_vs_B", "C_vs_D"], PARAMS, out)

    vdir = os.path.join(out, "Volcanos")
    assert os.path.isdir(vdir)
    for f in list(files) + [facet]:
        assert os.path.isfile(f), f
        assert os.path.dirname(os.path.abspath(f)) == os.path.abspath(vdir), f

    assert sorted(os.listdir(out)) == ["Volcanos"], \
        "out_dir doit être vide de PNG : tout est descendu dans Volcanos/"
    assert sorted(os.listdir(vdir)) == ["all_volcano_plots_grid.png",
                                        "volc_A_vs_B.png", "volc_C_vs_D.png"]


def test_an_empty_value_still_writes_at_the_root(tmp_path):
    out = str(tmp_path)
    files, facet = deimos.plot_volcanoes(_df(), ["A_vs_B"],
                                         {**PARAMS, "volcano_subdir": ""}, out)
    assert os.path.dirname(os.path.abspath(files[0])) == os.path.abspath(out)
    assert os.path.dirname(os.path.abspath(facet)) == os.path.abspath(out)
    assert not os.path.isdir(os.path.join(out, "Volcanos"))


def test_the_cleanup_removes_the_resized_copies_but_keeps_the_volcanos(tmp_path):
    """Reproduit ce que fait réellement l'étape 16 sur le sous-dossier."""
    out = str(tmp_path)
    deimos.plot_volcanoes(_df(), ["A_vs_B"], PARAMS, out)
    vdir = os.path.join(out, "Volcanos")

    # copies telles que _resize_png les nomme, à côté de leur original
    for n in ("volc_A_vs_B_xl560x440.png", "all_volcano_plots_grid_xl1100x1400.png"):
        open(os.path.join(vdir, n), "wb").close()

    removed = 0
    for fn in os.listdir(vdir):
        if not fn.lower().endswith((".png", ".tiff", ".tif")):
            continue
        if not ("_xl" in fn and "x" in fn.split("_xl")[-1]):
            continue
        os.remove(os.path.join(vdir, fn))
        removed += 1

    assert removed == 2
    assert sorted(os.listdir(vdir)) == ["all_volcano_plots_grid.png", "volc_A_vs_B.png"]
