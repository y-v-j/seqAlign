"""
heatmap.py
==========
When more than two sequences are supplied, the CLI (mirroring the original
`-a` exhaustive all-vs-all loop in main()) computes every pairwise
alignment once and renders it as a single annotated matrix:

    * upper triangle -> percent SIMILARITY  (its own color map)
    * lower triangle -> percent IDENTITY    (its own color map)
    * diagonal       -> a sequence vs. itself (left blank)

Every cell is annotated with its value, so no colorbar is drawn: the
colors only make patterns easy to spot at a glance.

This is new functionality relative to align.c (which only ever printed
plain-text pairwise results); it's provided here because comparing more
than a couple of sequences as scrolling text output doesn't scale, and a
heatmap is the standard way this kind of all-vs-all matrix is presented in
bioinformatics tooling (cf. MatGAT, Campanella et al. 2003).

Color themes
------------
Each theme pairs two color maps, one per triangle, so the two measures can
never be confused. See ``THEMES`` / ``seqalign --list-themes``; individual
color maps can also be overridden with ``--upper-cmap`` / ``--lower-cmap``.

Layout
------
The figure is laid out explicitly in inches rather than relying on
``tight_layout``: everything is drawn once, the real rendered extents of
the tick labels and title are measured, and then the
axes are positioned so that no two elements can overlap - however long
the sequence names are or however many sequences are compared.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence

import numpy as np


@dataclass(frozen=True)
class HeatmapTheme:
    name: str
    upper_cmap: str  # similarity
    lower_cmap: str  # identity
    description: str
    diagonal_color: str = "#f2f2f2"
    grid_color: str = "white"


THEMES: Dict[str, HeatmapTheme] = {
    t.name: t
    for t in [
        HeatmapTheme("sunset-ocean", "YlOrRd", "PuBu", "warm yellow-red vs. cool purple-blue (default)"),
        HeatmapTheme("classic", "Reds", "Blues", "red vs. blue"),
        HeatmapTheme("forest-berry", "YlGn", "RdPu", "yellow-green vs. pink-purple"),
        HeatmapTheme("perceptual", "viridis", "magma", "perceptually uniform viridis vs. magma",
                     diagonal_color="white"),
        HeatmapTheme("colorblind", "Oranges", "Purples", "orange vs. purple, safe for common color-vision deficiencies"),
        HeatmapTheme("grayscale", "Greys", "Greys_r", "print-friendly; dark-high vs. light-high grays",
                     diagonal_color="white", grid_color="#bbbbbb"),
    ]
}
DEFAULT_THEME = "sunset-ocean"


def describe_themes() -> str:
    """Human-readable list of themes (for --list-themes)."""
    lines = ["Heatmap color themes (upper triangle = similarity / lower triangle = identity):"]
    for t in THEMES.values():
        lines.append(f"  {t.name:<14} {t.upper_cmap:>8} / {t.lower_cmap:<8} {t.description}")
    lines.append("Override either half with --upper-cmap / --lower-cmap using any matplotlib colormap name.")
    return "\n".join(lines)


def get_theme(name: str = DEFAULT_THEME, upper_cmap: Optional[str] = None,
              lower_cmap: Optional[str] = None) -> HeatmapTheme:
    """Look up a theme by name, optionally overriding its color maps."""
    key = name.lower()
    if key not in THEMES:
        raise ValueError(f"Unknown heatmap theme '{name}'. Available: {', '.join(THEMES)}")
    theme = THEMES[key]
    if upper_cmap or lower_cmap:
        theme = HeatmapTheme(
            name=f"{theme.name} (custom)",
            upper_cmap=upper_cmap or theme.upper_cmap,
            lower_cmap=lower_cmap or theme.lower_cmap,
            description=theme.description,
            diagonal_color=theme.diagonal_color,
            grid_color=theme.grid_color,
        )
    return theme


def build_pairwise_matrices(
    labels: Sequence[str], identity: Sequence[Sequence[float]], similarity: Sequence[Sequence[float]]
) -> "np.ndarray":
    """
    Combine separately-computed identity/similarity values into one matrix
    for plotting: lower triangle = identity, upper triangle = similarity,
    diagonal = 100.
    """
    n = len(labels)
    combined = np.full((n, n), np.nan, dtype=np.float64)
    for i in range(n):
        combined[i, i] = 100.0
        for j in range(i + 1, n):
            combined[i, j] = similarity[i][j]  # upper: similarity
            combined[j, i] = identity[i][j]  # lower: identity
    return combined


# Font fallback chain: Roboto if the user's system/environment has it
# installed, otherwise matplotlib silently falls through to the next
# available match (Helvetica/Arial-family, then matplotlib's built-in
# DejaVu Sans) so the plot never errors out over a missing font file.
_FONT_FALLBACK_CHAIN = ["Roboto", "Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]

DEFAULT_TITLE = "Pairwise Sequence Comparison"

# Fixed spacing, in inches.
_MARGIN = 0.15
_GAP = 0.15  # between any two neighbouring elements


def _auto_rotation(text_widths_in: Sequence[float], cell_in: float) -> float:
    """Angle x labels once any of them is wider than its column."""
    return 45.0 if max(text_widths_in, default=0.0) > cell_in * 0.95 else 0.0


def plot_heatmap(
    labels: Sequence[str],
    combined_matrix: "np.ndarray",
    output_path: str,
    dpi: int = 150,
    title: str = DEFAULT_TITLE,
    theme: "HeatmapTheme | str" = DEFAULT_THEME,
    label_rotation: Optional[float] = None,
    annotate: bool = True,
    return_figure: bool = False,
):
    """
    Render the combined identity/similarity matrix as a heatmap.

    :param theme: a ``HeatmapTheme`` or theme name (see ``THEMES``).
    :param label_rotation: angle (degrees) for the x-axis labels; ``None``
        picks 0 for short names and 45 when any name is wider than its column.
    :param annotate: write the percentage inside each cell.
    :param return_figure: return the (still open) matplotlib Figure instead
        of closing it - used by the tests to inspect the layout.
    :returns: ``output_path`` (or the Figure if ``return_figure``).

    Raises a clear, actionable ImportError if matplotlib is missing.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")  # safe for headless/CI environments
        import matplotlib.pyplot as plt
        from matplotlib.colors import to_rgba
    except ImportError as exc:
        raise ImportError(
            "Heatmap generation requires matplotlib. Install it with:\n"
            "    pip install matplotlib"
        ) from exc

    if isinstance(theme, str):
        theme = get_theme(theme)

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = _FONT_FALLBACK_CHAIN

    n = len(labels)
    data = np.asarray(combined_matrix, dtype=np.float64)
    upper_mask = np.triu(np.ones((n, n), dtype=bool), k=1)
    lower_mask = np.tril(np.ones((n, n), dtype=bool), k=-1)

    # Cell size shrinks as n grows, within readable bounds.
    cell_in = float(np.clip(7.0 / max(n, 1), 0.42, 0.9))
    heat_in = cell_in * n
    tick_fs = float(np.clip(cell_in * 17, 7, 11))
    annot_fs = float(np.clip(cell_in * 15, 5.5, 10))

    fig = plt.figure(figsize=(heat_in + 4, heat_in + 3))
    ax = fig.add_axes([0.1, 0.1, 0.5, 0.5])

    images = {}
    for key, mask, cmap_name in (
        ("sim", upper_mask, theme.upper_cmap),
        ("ident", lower_mask, theme.lower_cmap),
    ):
        # Cells outside this triangle are masked and drawn fully transparent.
        cmap = plt.get_cmap(cmap_name).with_extremes(bad=(0, 0, 0, 0))
        images[key] = ax.imshow(
            np.ma.masked_array(data, mask=~mask), cmap=cmap, vmin=0, vmax=100,
            interpolation="nearest", aspect="equal",
        )

    ax.set_facecolor(theme.diagonal_color)  # only the diagonal shows through
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=tick_fs)
    ax.tick_params(length=0, pad=4)
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color=theme.grid_color, linewidth=1.2)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if annotate:
        for (i, j), value in np.ndenumerate(data):
            if i == j or np.isnan(value):
                continue
            img = images["sim"] if j > i else images["ident"]
            r, g, b, _ = to_rgba(img.cmap(img.norm(value)))
            luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
            ax.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=annot_fs,
                    color="white" if luminance < 0.5 else "#1a1a1a")

    subtitle = "upper triangle: % similarity   |   lower triangle: % identity"
    title_text = fig.text(0.5, 0.95, f"{title}\n{subtitle}" if title else subtitle,
                          ha="center", va="top", fontsize=tick_fs + 1.5, linespacing=1.5)

    # --- measure x labels to choose a rotation --------------------------------
    renderer = fig.canvas.get_renderer()
    ax.set_xticklabels(labels, fontsize=tick_fs)
    widths = [t.get_window_extent(renderer).width / fig.dpi for t in ax.get_xticklabels()]
    rotation = _auto_rotation(widths, cell_in) if label_rotation is None else float(label_rotation)
    if rotation:
        ax.set_xticklabels(labels, fontsize=tick_fs, rotation=rotation,
                           ha="right", va="top", rotation_mode="anchor")

    # --- measure every element, then lay out in inches -------------------------
    def extent(artists):
        boxes = [a.get_window_extent(renderer) for a in artists if a.get_visible() and a.get_text()]
        if not boxes:
            return 0.0, 0.0
        x0 = min(b.x0 for b in boxes); x1 = max(b.x1 for b in boxes)
        y0 = min(b.y0 for b in boxes); y1 = max(b.y1 for b in boxes)
        return (x1 - x0) / fig.dpi, (y1 - y0) / fig.dpi

    ylabel_w, _ = extent(ax.get_yticklabels())
    _, xlabel_h = extent(ax.get_xticklabels())
    pad_in = 4 / 72  # tick pad in points -> inches

    # A label anchored at column i's centre extends (its width) to the left;
    # make room if the leftmost one would stick out past the axes.
    xlabel_widths = [t.get_window_extent(renderer).width / fig.dpi for t in ax.get_xticklabels()]
    left_overhang = max(
        [w - (i + 0.5) * cell_in for i, w in enumerate(xlabel_widths)] + [0.0]
    ) if rotation else 0.0

    title_w, title_h = extent([title_text])

    ax_left = _MARGIN + max(ylabel_w + pad_in, left_overhang)
    ax_bottom = _MARGIN + xlabel_h + pad_in
    # Rotated labels can poke slightly right of their column centre.
    content_w = ax_left + heat_in + _GAP + _MARGIN
    fig_w = max(content_w, title_w + 2 * _MARGIN)
    title_top = ax_bottom + heat_in + _GAP * 2 + title_h
    fig_h = title_top + _MARGIN

    # Centre the content if the title made the figure wider.
    shift = (fig_w - content_w) / 2
    fig.set_size_inches(fig_w, fig_h)
    ax.set_position([(ax_left + shift) / fig_w, ax_bottom / fig_h, heat_in / fig_w, heat_in / fig_h])
    title_text.set_position((0.5, title_top / fig_h))

    fig.savefig(output_path, dpi=dpi, facecolor="white")
    if return_figure:
        return fig
    plt.close(fig)
    return output_path
