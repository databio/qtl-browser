"""Quick-look figures over the derived tables (not part of the build).

    uv run python -m pipeline.figures manhattan      # eQTL and sQTL lead-variant Manhattans -> data/figures/
"""
import sys

import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .common import CHROMS, Config, log, register_search_index

# dataviz reference palette: slot-1 blue for the series, two grays for alternating chromosomes
BLUE, GRAY_A, GRAY_B, INK, GRID = "#2a78d6", "#b8b7b1", "#d6d5cf", "#52514e", "#e1e0d9"


def chrom_offsets(cfg: Config) -> dict[str, tuple[int, int]]:
    """chr -> (offset, length) from the GRCh38 assembly report already in data/raw."""
    lengths = {}
    for line in cfg.assembly_report.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[1] == "assembled-molecule" and f[9] in CHROMS:
            lengths[f[9]] = int(f[8])
    out, pos = {}, 0
    for c in CHROMS:
        out[c] = (pos, lengths[c])
        pos += lengths[c]
    return out


def manhattan(cfg: Config) -> None:
    con = duckdb.connect()
    offs = chrom_offsets(cfg)
    outdir = cfg.derived.parent / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    specs = [
        ("eqtl", "cis-eQTL: one point per gene at its lead variant",
         f"SELECT chr, lead_position AS pos, pval_beta AS p, is_egene AS sig FROM '{cfg.tables / 'genes.parquet'}' WHERE tested"),
        ("sqtl", "cis-sQTL: one point per splice phenotype at its lead variant",
         f"SELECT chr, lead_position AS pos, pval_beta AS p, is_sqtl AS sig FROM '{cfg.tables / 'splice_phenotypes.parquet'}'"),
    ]
    for name, title, sql in specs:
        df = con.execute(sql).df()
        df = df[df.chr.isin(offs)]
        df["x"] = df.chr.map(lambda c: offs[c][0]) + df.pos
        df["y"] = -np.log10(df.p.clip(lower=1e-300))
        fig, ax = plt.subplots(figsize=(14, 4.2), dpi=150)
        for i, c in enumerate(CHROMS):
            o, L = offs[c]
            ax.axvspan(o, o + L, color=GRAY_B if i % 2 else "white", lw=0, zorder=0)
        ns = df[~df.sig.fillna(False)]
        s = df[df.sig.fillna(False)]
        ax.scatter(ns.x, ns.y, s=3, c=GRAY_A, lw=0, alpha=0.7, zorder=2, label=f"not significant ({len(ns):,})")
        ax.scatter(s.x, s.y, s=3, c=BLUE, lw=0, alpha=0.8, zorder=3, label=f"permutation p < 0.05 ({len(s):,})")
        ax.set_xlim(0, max(o + L for o, L in offs.values()))
        ax.set_xticks([o + L / 2 for o, L in offs.values()])
        ax.set_xticklabels([c.replace("chr", "") for c in CHROMS], fontsize=8, color=INK)
        ax.set_ylabel("−log10 permutation p (beta-approximated)", fontsize=9, color=INK)
        ax.set_title(title, fontsize=11, color=INK, loc="left")
        ax.tick_params(axis="y", labelsize=8, colors=INK)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.yaxis.grid(True, color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=8, loc="upper right", markerscale=3)
        fig.tight_layout()
        out = outdir / f"manhattan_{name}.png"
        fig.savefig(out)
        plt.close(fig)
        log(f"{out}  ({len(df):,} points; y max {df.y.max():.1f})")


VIOLET = "#7c5cbf"


def density(cfg: Config, bin_bp: int = 1_000_000) -> None:
    """Miami-style mock for the landing track: chromosome bar on top, downward density of
    eGenes and significant sQTL introns per window beneath it. From search_index only."""
    con = duckdb.connect()
    offs = chrom_offsets(cfg)
    outdir = cfg.derived.parent / "figures"
    register_search_index(cfg, con)
    df = con.execute("SELECT chr, tss, is_egene, n_sqtl_sig FROM search_index WHERE tested").df()
    df = df[df.chr.isin(offs)]
    df["x"] = df.chr.map(lambda c: offs[c][0]) + df.tss
    total = max(o + L for o, L in offs.values())
    edges = np.arange(0, total + bin_bp, bin_bp)
    e_counts, _ = np.histogram(df.x[df.is_egene.fillna(False)], bins=edges)
    s_counts, _ = np.histogram(df.x, bins=edges, weights=df.n_sqtl_sig)
    centers = (edges[:-1] + edges[1:]) / 2

    fig, ax = plt.subplots(figsize=(14, 3.2), dpi=150)
    bar_y, bar_h = 0.0, 1.2
    for i, c in enumerate(CHROMS):
        o, L = offs[c]
        ax.add_patch(plt.Rectangle((o, bar_y), L, bar_h, color=GRAY_A if i % 2 == 0 else GRAY_B, lw=0))
        ax.text(o + L / 2, bar_y + bar_h + 0.6, c.replace("chr", ""), ha="center", va="bottom", fontsize=7, color=INK)
    ax.bar(centers, -e_counts, width=bin_bp, color=BLUE, alpha=0.75, lw=0, bottom=bar_y - 0.4, label="eGenes per Mb")
    ax.bar(centers, -s_counts, width=bin_bp, color=VIOLET, alpha=0.45, lw=0, bottom=bar_y - 0.4, label="significant sQTL introns per Mb")
    ax.set_xlim(0, total)
    ymin = -(max(e_counts.max(), s_counts.max()) + 2)
    ax.set_ylim(ymin, bar_y + bar_h + 3)
    ax.set_xticks([])
    yt = [0, 25, 50, 100, 150, 200, 300, 400, 600]
    ax.set_yticks([bar_y - 0.4 - t for t in yt if -t > ymin])
    ax.set_yticklabels([str(t) for t in yt if -t > ymin], fontsize=7, color=INK)
    for sp in ("top", "right", "bottom"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color(GRID)
    ax.yaxis.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title(f"Landing-track mock: chromosome bar with downward QTL density ({bin_bp // 1_000_000} Mb bins)", fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    out = outdir / f"track_density_mock_{bin_bp // 1_000_000}mb.png"
    fig.savefig(out)
    plt.close(fig)
    log(f"{out}  (max eGenes/Mb {e_counts.max()}, max sQTL introns/Mb {s_counts.max():.0f})")


RED = "#e34948"


def gwas_bins(cfg: Config, bin_bp: int = 5_000_000):
    """Strongest DCM GWAS signal per window, straight from the raw Jurgens meta-analysis."""
    con = duckdb.connect()
    src = cfg.raw / "dcm_gwas_jurgens2024" / "DCM_GWAS" / "Jurgens_DCM_GWAS_META.tsv.gz"
    return con.execute(f"""
        SELECT 'chr' || split_part(CHRBP_B38, ':', 1) AS chr,
               (split_part(CHRBP_B38, ':', 2)::BIGINT // {bin_bp}) * {bin_bp} AS bin_start,
               min(P) AS min_p, count(*) FILTER (WHERE P < 5e-8) AS n_gws, count(*) AS n
        FROM read_csv('{src}', delim='\\t', header=true, columns={{'CHRBP_B37':'VARCHAR','CHRBP_B38':'VARCHAR','ID_B38':'VARCHAR','CHR':'VARCHAR','POS':'BIGINT','EA':'VARCHAR','NEA':'VARCHAR','BETA':'DOUBLE','SE':'DOUBLE','P':'DOUBLE','EAFREQ':'DOUBLE','HetDf':'INTEGER','HetPVal':'DOUBLE','N':'BIGINT','N_cases':'BIGINT','N_controls':'BIGINT','rsID':'VARCHAR','ID_B37':'VARCHAR','INDEL':'VARCHAR'}})
        WHERE CHRBP_B38 IS NOT NULL AND CHRBP_B38 <> '' AND P > 0
        GROUP BY 1, 2 ORDER BY 1, 2""").df()


def gwas_density(cfg: Config, bin_bp: int = 5_000_000) -> None:
    """Miami mock: chromosome bar, DCM GWAS strongest -log10 p per window pointing down."""
    offs = chrom_offsets(cfg)
    outdir = cfg.derived.parent / "figures"
    df = gwas_bins(cfg, bin_bp)
    df = df[df.chr.isin(offs)]
    df["x"] = df.chr.map(lambda c: offs[c][0]) + df.bin_start + bin_bp / 2
    df["y"] = -np.log10(df.min_p.clip(lower=1e-300))
    total = max(o + L for o, L in offs.values())
    fig, ax = plt.subplots(figsize=(14, 3.2), dpi=150)
    bar_y, bar_h = 0.0, 0.5
    for i, c in enumerate(CHROMS):
        o, L = offs[c]
        ax.add_patch(plt.Rectangle((o, bar_y), L, bar_h, color=GRAY_A if i % 2 == 0 else GRAY_B, lw=0))
        ax.text(o + L / 2, bar_y + bar_h + 0.3, c.replace("chr", ""), ha="center", va="bottom", fontsize=7, color=INK)
    gws = df.n_gws > 0
    ax.bar(df.x[~gws], -df.y[~gws], width=bin_bp * 0.9, color=GRAY_A, lw=0, bottom=bar_y - 0.2, label="strongest variant in window")
    ax.bar(df.x[gws], -df.y[gws], width=bin_bp * 0.9, color=RED, lw=0, bottom=bar_y - 0.2, label="window holds a genome-wide significant hit")
    thr = -np.log10(5e-8)
    ax.axhline(bar_y - 0.2 - thr, color=RED, lw=0.6, ls=(0, (3, 3)), alpha=0.6)
    ax.set_xlim(0, total)
    ax.set_ylim(-(df.y.max() + 2), bar_y + bar_h + 1.5)
    ax.set_xticks([])
    yt = [0, 5, 10, 20, 30, 40]
    ax.set_yticks([bar_y - 0.2 - t for t in yt if t < df.y.max() + 2])
    ax.set_yticklabels([str(t) for t in yt if t < df.y.max() + 2], fontsize=7, color=INK)
    ax.set_ylabel("−log10 p", fontsize=8, color=INK)
    for sp in ("top", "right", "bottom"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color(GRID)
    ax.yaxis.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title(f"Landing-track mock: DCM GWAS (Jurgens 2024) strongest signal per {bin_bp // 1_000_000} Mb window", fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    out = outdir / f"track_gwas_mock_{bin_bp // 1_000_000}mb.png"
    fig.savefig(out)
    plt.close(fig)
    log(f"{out}  ({len(df)} windows, {int(gws.sum())} with GWS hits, y max {df.y.max():.1f})")


THEMES = {
    # name: (light tokens, dark tokens); tokens = surface, surface2, border, ink, muted, primary, secondary, accent, error
    "1 Myocardium": (
        dict(bg="#fbf8f5", bg2="#f3ede8", border="#e4dbd4", ink="#2b1d1a", muted="#7a6b66", primary="#7a1f2b", secondary="#3f5f8a", accent="#c9788a", error="#c0392b"),
        dict(bg="#1f1615", bg2="#191110", border="#332624", ink="#f1e9e4", muted="#a3948e", primary="#d97b88", secondary="#8fb0dc", accent="#e3a3b2", error="#e3706a"),
    ),
    "2 Pâtisserie": (
        dict(bg="#fcf9f1", bg2="#f5efe2", border="#e6dcc8", ink="#2a211b", muted="#7b6f62", primary="#9a5b1f", secondary="#1f6f6b", accent="#d9a441", error="#c0392b"),
        dict(bg="#1e1a15", bg2="#181410", border="#332c22", ink="#f3ede2", muted="#a89c8c", primary="#e0a060", secondary="#6fc2bd", accent="#e8c36a", error="#e3706a"),
    ),
    "3 Clinical": (
        dict(bg="#fafbfb", bg2="#f1f4f4", border="#dfe5e6", ink="#1f2428", muted="#6b7378", primary="#0f6b6e", secondary="#4a4fb0", accent="#e0694f", error="#c0392b"),
        dict(bg="#151a1c", bg2="#101416", border="#25302f", ink="#e9eeee", muted="#97a2a6", primary="#6ec6c8", secondary="#a3a7f0", accent="#f19a82", error="#e3706a"),
    ),
    "4 Ink & ochre": (
        dict(bg="#f6f4ef", bg2="#eeebe3", border="#dcd7cc", ink="#22211f", muted="#6f6c65", primary="#b5761a", secondary="#1e4b73", accent="#7d3c98", error="#c0392b"),
        dict(bg="#1a1917", bg2="#141311", border="#2d2b27", ink="#efece5", muted="#a19d94", primary="#e6a34c", secondary="#7fb0e0", accent="#c48ad8", error="#e3706a"),
    ),
}


def theme_mocks(cfg: Config) -> None:
    """One mock card per scheme, light and dark: wordmark, title, description, marker colors,
    badges, a primary button, the GWAS error bar. For eyeballing, not the real UI."""
    from matplotlib.patches import Polygon, Rectangle
    outdir = cfg.derived.parent / "figures"
    n = len(THEMES)
    fig, axes = plt.subplots(n, 2, figsize=(13, 2.6 * n), dpi=150)
    for row, (name, (light, dark)) in enumerate(THEMES.items()):
        for col, t in enumerate((light, dark)):
            ax = axes[row][col]
            ax.set_xlim(0, 100); ax.set_ylim(0, 40); ax.axis("off")
            ax.add_patch(Rectangle((0, 0), 100, 40, color=t["bg"], lw=0))
            # navbar
            ax.add_patch(Rectangle((0, 34), 100, 6, color=t["bg"], lw=0))
            ax.plot([0, 100], [34, 34], color=t["border"], lw=0.8)
            ax.text(3, 36.4, "topchef", color=t["muted"], fontsize=10, va="center", fontweight="light")
            ax.text(13.5, 36.4, ".", color=t["border"], fontsize=10, va="center")
            ax.text(14.6, 36.4, "qtl", color=t["primary"], fontsize=10, va="center")
            for i, lab in enumerate(("Search", "Genes", "About")):
                ax.text(72 + i * 9, 36.4, lab, color=t["primary"] if i == 0 else t["muted"], fontsize=7.5, va="center")
            # hero
            ax.text(3, 29.5, "TOPCHeF", color=t["ink"], fontsize=15, va="center", fontweight="light")
            ax.text(3, 25.3, "Expression and splicing QTL mapped in left-ventricle tissue from failing and non-failing human hearts.", color=t["muted"], fontsize=6.3, va="center")
            ax.text(3, 23.0, "Search a gene, variant, or region for summary statistics, fine-mapping, and locus views.", color=t["muted"], fontsize=6.3, va="center")
            ax.text(66.5, 23.0, "Preprint", color=t["primary"], fontsize=6.3, va="center")
            ax.text(72.5, 23.0, "Zenodo", color=t["primary"], fontsize=6.3, va="center")
            # search bar with primary circular button
            ax.add_patch(Rectangle((3, 16.5), 60, 4.2, facecolor=t["bg"], edgecolor=t["primary"], lw=1.0, alpha=0.9))
            ax.text(5, 18.6, "Gene symbol, Ensembl ID, rsID, chr:pos, or chr:start-end", color=t["muted"], fontsize=6, va="center", alpha=0.7)
            ax.add_patch(plt.Circle((61.2, 18.6), 1.5, color=t["primary"], lw=0))
            # markers + chromosome bar
            ax.add_patch(Rectangle((3, 9.6), 60, 0.9, color=t["border"], lw=0))
            ax.add_patch(Rectangle((23, 9.6), 20, 0.9, color=t["muted"], alpha=0.45, lw=0))
            for x, key, lab in ((12, "primary", "eQTL"), (33, "secondary", "sQTL"), (52, "accent", "both")):
                ax.add_patch(Polygon([(x - 1.2, 12.6), (x + 1.2, 12.6), (x, 10.8)], color=t[key], lw=0))
                ax.text(x, 13.4, lab, color=t["ink"], fontsize=5.5, ha="center", va="bottom", alpha=0.8)
            # downward GWAS bars
            rng = np.random.default_rng(3)
            for i in range(60):
                h = rng.gamma(1.2, 0.6)
                sig = h > 2.4
                ax.add_patch(Rectangle((3 + i, 9.3 - min(h, 4)), 0.8, min(h, 4), color=t["error"] if sig else t["muted"], alpha=1 if sig else 0.35, lw=0))
            # badges + button
            for i, (key, lab) in enumerate((("primary", "eGene"), ("secondary", "3 sQTL introns"), ("accent", "DCM coloc · eQTL"))):
                x = 68 + i * 10.5
                ax.add_patch(Rectangle((x, 27.6), 9.6, 2.8, color=t[key], lw=0))
                ax.text(x + 4.8, 29.0, lab, color=t["bg"], fontsize=4.6, ha="center", va="center", fontweight="bold")
            ax.add_patch(Rectangle((68, 16.5), 12, 4.2, color=t["primary"], lw=0))
            ax.text(74, 18.6, "New query", color=t["bg"], fontsize=6, ha="center", va="center", fontweight="bold")
            ax.add_patch(Rectangle((81.5, 16.5), 8, 4.2, facecolor=t["bg2"], edgecolor=t["border"], lw=0.8))
            ax.text(85.5, 18.6, "CSV", color=t["ink"], fontsize=6, ha="center", va="center")
            # table sample
            ax.add_patch(Rectangle((68, 3), 29, 10, facecolor=t["bg"], edgecolor=t["border"], lw=0.8))
            ax.add_patch(Rectangle((68, 8), 29, 2.5, color=t["bg2"], lw=0))
            for i, (a, b) in enumerate((("Gene", "Perm p"), ("FLNC", "3.1e-3"), ("SKI", "1.0e-4"), ("CAMK2D", "1.0e-4"))):
                y = 11.8 - i * 2.5
                ax.text(69, y, a, color=t["primary"] if i else t["muted"], fontsize=5.5, va="center")
                ax.text(96, y, b, color=t["ink"], fontsize=5.5, va="center", ha="right")
            ax.text(3, 2.2, f"{name}  ·  {'light' if col == 0 else 'dark'}", color=t["muted"], fontsize=6.5, va="center")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, hspace=0.06, wspace=0.03)
    out = outdir / "theme_mocks.png"
    fig.savefig(out)
    plt.close(fig)
    log(str(out))


if __name__ == "__main__":
    cfg = Config()
    cmd = sys.argv[1] if len(sys.argv) > 1 else "manhattan"
    if cmd == "themes":
        theme_mocks(cfg)
        sys.exit(0)
    if cmd == "manhattan":
        manhattan(cfg)
    elif cmd == "density":
        for mb in (int(a) for a in (sys.argv[2:] or ["1"])):
            density(cfg, bin_bp=mb * 1_000_000)
    elif cmd == "gwas":
        for mb in (int(a) for a in (sys.argv[2:] or ["5"])):
            gwas_density(cfg, bin_bp=mb * 1_000_000)
    else:
        sys.exit(f"unknown command {cmd}")
