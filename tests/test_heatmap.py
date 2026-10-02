import itertools

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")

from seqalign.heatmap import THEMES, build_pairwise_matrices, get_theme, plot_heatmap  # noqa: E402


def _matrix(n, seed=0):
    rng = np.random.default_rng(seed)
    ident = rng.uniform(20, 95, size=(n, n))
    ident = (ident + ident.T) / 2
    sim = np.minimum(ident + 10, 100)
    labels = [f"seq{i}" for i in range(n)]
    return labels, build_pairwise_matrices(labels, ident, sim)


def _layout_groups(fig):
    """Bounding boxes (display coords) of each separately-placed element of the heatmap."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    (heat_ax,) = fig.axes  # no colorbars: values are written in the cells

    def boxes(artists):
        return [a.get_window_extent(renderer) for a in artists if a.get_visible() and a.get_text()]

    return {
        "title": boxes(fig.texts),
        "heatmap": [heat_ax.get_window_extent(renderer)],
        "xlabels": boxes(heat_ax.get_xticklabels()),
        "ylabels": boxes(heat_ax.get_yticklabels()),
    }


def _assert_no_overlap(fig):
    groups = _layout_groups(fig)
    for (name_a, boxes_a), (name_b, boxes_b) in itertools.combinations(groups.items(), 2):
        for a in boxes_a:
            for b in boxes_b:
                overlap_w = min(a.x1, b.x1) - max(a.x0, b.x0)
                overlap_h = min(a.y1, b.y1) - max(a.y0, b.y0)
                assert overlap_w <= 0.5 or overlap_h <= 0.5, f"{name_a} overlaps {name_b}"
    # Everything must lie inside the saved figure.
    fig_box = fig.bbox
    for name, boxes in groups.items():
        for b in boxes:
            assert b.x0 >= -0.5 and b.y0 >= -0.5 and b.x1 <= fig_box.x1 + 0.5 and b.y1 <= fig_box.y1 + 0.5, name


def test_build_pairwise_matrices_triangles():
    labels = ["a", "b", "c"]
    ident = [[0, 10, 20], [10, 0, 30], [20, 30, 0]]
    sim = [[0, 15, 25], [15, 0, 35], [25, 35, 0]]
    m = build_pairwise_matrices(labels, ident, sim)
    assert m[0, 1] == 15 and m[1, 0] == 10  # upper = similarity, lower = identity
    assert m[1, 2] == 35 and m[2, 1] == 30
    assert (np.diag(m) == 100).all()


@pytest.mark.parametrize("theme", list(THEMES))
def test_every_theme_renders(theme, tmp_path):
    labels, combined = _matrix(4)
    out = tmp_path / f"{theme}.png"
    plot_heatmap(labels, combined, str(out), dpi=60, theme=theme)
    assert out.stat().st_size > 0


def test_triangles_use_different_colormaps(tmp_path):
    labels, combined = _matrix(5)
    fig = plot_heatmap(labels, combined, str(tmp_path / "h.png"), dpi=60, return_figure=True)
    images = fig.axes[0].get_images()
    assert len(images) == 2
    assert images[0].cmap.name != images[1].cmap.name
    # Each image shows exactly one triangle; the other cells are masked.
    upper, lower = (np.ma.getmaskarray(img.get_array()) for img in images)
    n = len(labels)
    assert (~upper == np.triu(np.ones((n, n), bool), 1)).all()
    assert (~lower == np.tril(np.ones((n, n), bool), -1)).all()


def test_theme_overrides():
    theme = get_theme("classic", upper_cmap="viridis")
    assert theme.upper_cmap == "viridis" and theme.lower_cmap == THEMES["classic"].lower_cmap
    with pytest.raises(ValueError):
        get_theme("no-such-theme")


def test_short_labels_are_horizontal_long_labels_are_angled(tmp_path):
    labels, combined = _matrix(4)
    fig = plot_heatmap(labels, combined, str(tmp_path / "a.png"), dpi=60, return_figure=True)
    assert all(t.get_rotation() == 0 for t in fig.axes[0].get_xticklabels())

    long_labels = [f"Organism_with_a_rather_long_name_{i}" for i in range(4)]
    fig = plot_heatmap(long_labels, combined, str(tmp_path / "b.png"), dpi=60, return_figure=True)
    assert all(t.get_rotation() == 45 for t in fig.axes[0].get_xticklabels())

    fig = plot_heatmap(long_labels, combined, str(tmp_path / "c.png"), dpi=60,
                       label_rotation=90, return_figure=True)
    assert all(t.get_rotation() == 90 for t in fig.axes[0].get_xticklabels())


@pytest.mark.parametrize(
    "n,label_len,rotation,title",
    [
        (3, 5, None, "Short"),
        (4, 60, None, "Pairwise Sequence Comparison"),
        (15, 25, None, "Fifteen organisms"),
        (15, 40, 90, "Fifteen organisms, vertical labels"),
        (30, 12, None, "Thirty sequences"),
        (3, 8, None, "A very long title that is much wider than a three-by-three heatmap would be on its own"),
    ],
)
def test_layout_has_no_overlaps(n, label_len, rotation, title, tmp_path):
    rng = np.random.default_rng(n)
    labels = ["".join(rng.choice(list("abcdefghij_"), size=label_len)) + str(i) for i in range(n)]
    _, combined = _matrix(n)
    fig = plot_heatmap(labels, combined, str(tmp_path / "h.png"), dpi=72, title=title,
                       label_rotation=rotation, return_figure=True)
    _assert_no_overlap(fig)
