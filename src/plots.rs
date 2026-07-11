//! Optional PDF plots for the complexity metrics.
//!
//! Rendered with [`kuva`] when `--complexity-metrics` is set, alongside the TSVs
//! in [`crate::complexity`] and [`crate::counts`]. Everything here is plot
//! presentation only — the numbers come from the TSV row types.
//!
//! Two plots:
//! - **Duplication-sampled ladder** ([`write_duplicate_ladder_pdf`]): the marginal
//!   (per-window) duplication rate vs. depth, one line per library on a shared
//!   axis. The cumulative rate is left to the TSV.
//! - **Group-size histogram (η_k)** ([`write_count_histogram_pdfs`]): the fraction
//!   of molecules and of reads/pairs seen exactly `k` times (log-y / linear-x),
//!   with the heavy tail trimmed and summarized in the x-label. Always one file;
//!   with several libraries it is a faceted grid, one panel per library.
//!
//! Colours follow the Fulcrum Genomics brand palette.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, anyhow};
use kuva::backend::pdf::PdfBackend;
use kuva::plot::{LegendPosition, LinePlot, MarkerShape, ScatterPlot};
use kuva::render::figure::Figure;
use kuva::render::layout::{Layout, TickFormat};
use kuva::render::plots::Plot;

use crate::complexity::DuplicateLadderRow;
use crate::counts::CountHistogramRow;

/// Fulcrum Genomics categorical palette (hex for kuva), in fixed assignment
/// order. One colour per library; cycled only in the rare case of more libraries
/// than entries (a dedup run with >8 libraries).
const FG_PALETTE: [&str; 8] = [
    "#26a8e0", // FG Blue
    "#38b44a", // FG Green
    "#160052", // FG Navy
    "#1693b9", // FG Pacific
    "#315848", // FG Pine
    "#2fae99", // FG Teal
    "#4dcc68", // FG Emerald
    "#269e2a", // FG Forest
];

/// The two η_k series colours: Molecules (FG Blue), Reads/Pairs (FG Green).
const FG_BLUE: &str = FG_PALETTE[0];
const FG_GREEN: &str = FG_PALETTE[1];

// Standard plot dimensions (points; 8 × 6 inches at 72 DPI).
const PLOT_WIDTH: f64 = 800.0;
const PLOT_HEIGHT: f64 = 600.0;

/// x-axis is scaled to millions of templates so tick labels read "50", "100", …
/// (the axis label carries the unit) instead of `5e7` scientific notation.
const TEMPLATES_PER_UNIT: f64 = 1e6;

/// The η_k x-axis is trimmed at `min(this, <last k where [1..k] is ≥50% populated>)`
/// so heavy-tailed libraries (e.g. RNA) stay readable; the tail is folded into the
/// x-axis label rather than clipped away.
const HIST_CUTOFF_CAP: f64 = 500.0;

/// `<prefix>.duplication-sampled.pdf` (sibling to the `.tsv`).
pub fn ladder_plot_path(prefix: &Path) -> PathBuf {
    let mut name = prefix.as_os_str().to_owned();
    name.push(".duplication-sampled.pdf");
    PathBuf::from(name)
}

/// The marginal (per-window) series for one library's rows (ascending `total`),
/// as `(templates_in_millions, window_frac_duplicates)`. The final point is
/// dropped when it covers a short (partial) window so the line doesn't crash at
/// the right edge.
fn marginal_series(rows: &[&DuplicateLadderRow]) -> Vec<(f64, f64)> {
    let n = rows.len();
    if n == 0 {
        return Vec::new();
    }
    // The first snapshot lands at exactly one interval; use it as the reference.
    let interval = rows[0].total;
    let end = if n >= 2 && rows[n - 1].total - rows[n - 2].total < interval { n - 1 } else { n };
    rows[..end]
        .iter()
        .map(|r| (r.total as f64 / TEMPLATES_PER_UNIT, r.window_frac_duplicates))
        .collect()
}

/// Write the duplication-sampled ladder as a single PDF: the marginal (per-window)
/// duplication rate as one line per library, on a shared fraction-duplicated
/// y-axis vs. templates processed (in millions). A legend is shown only when
/// there is more than one library (a lone library is named by the title). `rows`
/// are as [`crate::complexity::LadderRecorder`] emits them (one category per
/// library, sorted by library then `total`). Does nothing (returns `Ok`) when
/// `rows` is empty.
///
/// # Errors
/// Returns an error if PDF rendering fails or the file cannot be written.
pub fn write_duplicate_ladder_pdf(rows: &[DuplicateLadderRow], path: &Path) -> Result<()> {
    let Some(first) = rows.first() else {
        return Ok(());
    };

    // Distinct libraries in first-seen (sorted) order.
    let mut libraries: Vec<&str> = Vec::new();
    for r in rows {
        if !libraries.contains(&r.library.as_str()) {
            libraries.push(&r.library);
        }
    }
    let multi = libraries.len() > 1;

    let mut plots: Vec<Plot> = Vec::new();
    let (mut x_min, mut x_max) = (f64::INFINITY, 0.0_f64);
    for (i, &lib) in libraries.iter().enumerate() {
        let lib_rows: Vec<&DuplicateLadderRow> =
            rows.iter().filter(|r| r.library.as_str() == lib).collect();
        let series = marginal_series(&lib_rows);
        if let Some(&(x, _)) = series.first() {
            x_min = x_min.min(x);
        }
        if let Some(&(x, _)) = series.last() {
            x_max = x_max.max(x);
        }
        let mut line = LinePlot::new()
            .with_data(series)
            .with_color(FG_PALETTE[i % FG_PALETTE.len()])
            .with_stroke_width(2.5);
        if multi {
            line = line.with_legend(lib);
        }
        plots.push(Plot::Line(line));
    }

    let title = plot_title(&first.sample, "Duplication Rate vs. Depth");
    let mut layout = Layout::auto_from_plots(&plots)
        .with_width(PLOT_WIDTH)
        .with_height(PLOT_HEIGHT)
        .with_title(&title)
        .with_x_label("Templates Processed (Millions)")
        .with_y_label("Marginal Fraction Duplicated")
        .with_y_tick_format(TickFormat::Percent)
        .with_y_axis_min(0.0);
    // Pin the x range to the data (first point → deepest total) to drop the empty
    // bands auto-ranging leaves. Guard degenerate single-point input.
    if x_max > x_min {
        layout = layout.with_x_axis_min(x_min).with_x_axis_max(x_max);
    }
    if multi {
        layout = layout.with_legend_position(LegendPosition::InsideTopRight);
    }

    let pdf_bytes = kuva::render_to_pdf(plots, layout).map_err(|e| anyhow!("{e}"))?;
    std::fs::write(path, pdf_bytes)
        .with_context(|| format!("writing duplication-sampled plot to {}", path.display()))?;
    Ok(())
}

// ── Group-size histogram (η_k) ──────────────────────────────────────────────

/// `<prefix>.duplication-spectrum.pdf` (sibling to the `.tsv`). One file for the whole
/// run — multiple libraries are faceted into it, never split across files.
fn count_histogram_plot_path(prefix: &Path) -> PathBuf {
    let mut name = prefix.as_os_str().to_owned();
    name.push(".duplication-spectrum.pdf");
    PathBuf::from(name)
}

/// x-axis cutoff `k`: extend while ≥50% of bins in `[1..k]` are populated, cap at
/// [`HIST_CUTOFF_CAP`] but never below the smallest observed `k`, round up to a
/// multiple of 10. `rows` are ascending in `n_observations` (as
/// [`crate::counts::histogram_rows`] emits per library).
///
/// The cap trims a heavy tail while keeping the head. If *every* signature was
/// observed more than the cap (a pathologically low-complexity library — say a
/// tiny amplicon panel where nothing is seen ≤500×), capping at the cap would
/// fold all points into the tail and leave nothing to plot; the smallest-`k`
/// floor keeps at least the first point in range so the plot is never empty.
fn histogram_cutoff(rows: &[&CountHistogramRow]) -> f64 {
    let Some(first) = rows.first() else {
        return 0.0;
    };
    let first_k = first.n_observations as f64;
    let mut raw = first_k;
    for (i, r) in rows.iter().enumerate().skip(1) {
        let k = r.n_observations as f64;
        if (i as f64 + 1.0) / k >= 0.5 {
            raw = k;
        } else {
            break;
        }
    }
    let capped = raw.min(HIST_CUTOFF_CAP).max(first_k);
    (capped / 10.0).ceil() * 10.0
}

/// Plain-decimal percent for a log-axis tick value (fraction → "%"); kuva's
/// built-in [`TickFormat::Percent`] is fixed at one decimal, which collapses the
/// small decades on a log axis to "0.0%".
fn pct_tick(v: f64) -> String {
    let p = v * 100.0;
    if p <= 0.0 {
        return "0%".to_string();
    }
    let e = p.log10().round() as i32;
    if e >= 0 { format!("{p:.0}%") } else { format!("{:.*}%", (-e) as usize, p) }
}

/// A fraction as a compact percentage string (plain decimal); values below
/// 0.001% are floored to "<0.001%" so a negligible tail doesn't render as
/// "0.0000041%".
fn pct_note(frac: f64) -> String {
    let p = frac * 100.0;
    if p > 0.0 && p < 0.001 {
        return "<0.001%".to_string();
    }
    if p >= 10.0 {
        return format!("{p:.0}%");
    }
    if p >= 1.0 {
        return format!("{p:.1}%");
    }
    let decimals = ((-p.log10().floor()) as usize) + 1;
    let mut s = format!("{p:.decimals$}");
    while s.ends_with('0') {
        s.pop();
    }
    if s.ends_with('.') {
        s.pop();
    }
    format!("{s}%")
}

/// Build the η_k panel for one library's `rows`: point series for the fraction of
/// **molecules** and of **reads/pairs** (`k·η_k`) observed exactly `k` times, as
/// cross markers on a log-y / linear-x percent plot titled `title`. The x-axis is
/// trimmed at [`histogram_cutoff`]; molecules/reads beyond it are summarized in the
/// x-axis label.
///
/// Returns the `(plots, layout)` to render — via [`kuva::render_to_pdf`] for a lone
/// library, or as one cell of a [`Figure`] grid for several — or `None` when there
/// is nothing plottable (empty rows, or every point folded past the cutoff). `rows`
/// must be a single `(library, category)` series ascending in `n_observations`.
fn count_histogram_panel(rows: &[&CountHistogramRow], title: &str) -> Option<(Vec<Plot>, Layout)> {
    rows.first()?;
    let cutoff = histogram_cutoff(rows);
    let total_mol: f64 = rows.iter().map(|r| r.n_molecules as f64).sum();
    let total_rp: f64 = rows.iter().map(|r| r.n_observations as f64 * r.n_molecules as f64).sum();
    if total_mol == 0.0 || total_rp == 0.0 {
        return None;
    }

    // Points up to the cutoff; molecules/reads beyond it fold into the label.
    let mut mol: Vec<(f64, f64)> = Vec::new();
    let mut rp: Vec<(f64, f64)> = Vec::new();
    let (mut mol_tail, mut rp_tail) = (0.0_f64, 0.0_f64);
    for r in rows {
        let (k, n) = (r.n_observations as f64, r.n_molecules as f64);
        if k <= cutoff {
            mol.push((k, n / total_mol));
            rp.push((k, k * n / total_rp));
        } else {
            mol_tail += n;
            rp_tail += k * n;
        }
    }

    // Defensive: `histogram_cutoff` floors at the smallest `k`, so at least one
    // point is always in range — but never hand kuva empty series (fmin would be
    // +∞, producing an infinite log-y bound and a render panic) if that ever
    // changes.
    if mol.is_empty() {
        return None;
    }

    let fmin = mol.iter().chain(rp.iter()).map(|&(_, f)| f).fold(f64::INFINITY, f64::min);
    let fmax = mol.iter().chain(rp.iter()).map(|&(_, f)| f).fold(0.0, f64::max);
    let y_min = 10f64.powf(fmin.log10() - 0.4);
    let y_max = 10f64.powf(fmax.log10() + 0.25);

    let cutoff_i = cutoff as u64;
    let x_label = if mol_tail > 0.0 {
        format!(
            "Times Observed (k)   ({} molecules, {} reads @ k>{cutoff_i})",
            pct_note(mol_tail / total_mol),
            pct_note(rp_tail / total_rp),
        )
    } else {
        "Times Observed (k)".to_string()
    };

    let points = |data: Vec<(f64, f64)>, color: &str, label: &str| {
        Plot::Scatter(
            ScatterPlot::new()
                .with_data(data)
                .with_color(color)
                .with_size(5.0)
                .with_marker(MarkerShape::Cross)
                .with_legend(label),
        )
    };
    let plots: Vec<Plot> =
        vec![points(mol, FG_BLUE, "Molecules"), points(rp, FG_GREEN, "Reads/Pairs")];

    // `width`/`height` are honored by `render_to_pdf` (lone library) and overridden
    // by `Figure` when this is one facet cell — set once here, correct for both.
    let mut layout = Layout::auto_from_plots(&plots)
        .with_width(PLOT_WIDTH)
        .with_height(PLOT_HEIGHT)
        .with_title(title)
        .with_x_label(&x_label)
        .with_y_label("Percent")
        .with_log_y()
        .with_x_axis_min(0.0)
        .with_x_axis_max(cutoff)
        .with_y_axis_min(y_min)
        .with_minor_ticks(9)
        .with_show_minor_grid(true)
        .with_y_tick_format(TickFormat::Custom(Arc::new(pct_tick)))
        .with_legend_position(LegendPosition::InsideTopRight);
    if y_max > y_min {
        layout = layout.with_y_axis_max(y_max);
    }
    Some((plots, layout))
}

/// A placeholder [`Layout`] for an empty facet cell (a padded trailing grid slot):
/// no data, it just satisfies the index-aligned `Figure::with_layouts` contract.
fn blank_cell_layout() -> Layout {
    Layout::new((0.0, 1.0), (0.0, 1.0))
}

/// Plot/figure title, dropping the `sample — ` prefix when there is no sample
/// name (e.g. a BAM with no `@RG SM:`), so it never renders a stray leading dash.
fn plot_title(sample: &str, what: &str) -> String {
    if sample.is_empty() { what.to_string() } else { format!("{sample} — {what}") }
}

/// Write the group-size histogram (η_k) as a single `<prefix>.duplication-spectrum.pdf`
/// (rows as [`crate::counts::histogram_rows`] emits them: one category per library,
/// grouped by library, ascending in `n_observations`). A lone library renders one
/// plot; several libraries render a **faceted grid**, one panel per library, on the
/// single page. Does nothing when `rows` is empty.
///
/// # Errors
/// Returns an error if PDF rendering fails or the file cannot be written.
pub fn write_count_histogram_pdfs(rows: &[CountHistogramRow], prefix: &Path) -> Result<()> {
    // Distinct libraries in first-seen order (matches the TSV row order).
    let mut libraries: Vec<&str> = Vec::new();
    for r in rows {
        if !libraries.contains(&r.library.as_str()) {
            libraries.push(&r.library);
        }
    }
    if libraries.is_empty() {
        return Ok(());
    }
    let path = count_histogram_plot_path(prefix);
    let sample = rows.first().map_or("", |r| r.sample.as_str());
    let title = plot_title(sample, "Duplication Spectrum (η_k)");

    // Lone library: a single plot, exactly as before.
    if libraries.len() == 1 {
        let series: Vec<&CountHistogramRow> = rows.iter().collect();
        if let Some((plots, layout)) = count_histogram_panel(&series, &title) {
            let pdf = kuva::render_to_pdf(plots, layout).map_err(|e| anyhow!("{e}"))?;
            std::fs::write(&path, pdf).with_context(|| {
                format!("writing duplication-spectrum plot to {}", path.display())
            })?;
        }
        return Ok(());
    }

    // Several libraries: one faceted page, a panel per library titled by library.
    let cols = (libraries.len() as f64).sqrt().ceil() as usize;
    let n_rows = libraries.len().div_ceil(cols);
    let cell_count = n_rows * cols;
    let mut cell_plots: Vec<Vec<Plot>> = Vec::with_capacity(cell_count);
    let mut cell_layouts: Vec<Layout> = Vec::with_capacity(cell_count);
    for &lib in &libraries {
        let series: Vec<&CountHistogramRow> =
            rows.iter().filter(|r| r.library.as_str() == lib).collect();
        match count_histogram_panel(&series, lib) {
            Some((plots, layout)) => {
                cell_plots.push(plots);
                cell_layouts.push(layout);
            }
            None => {
                cell_plots.push(Vec::new());
                cell_layouts.push(blank_cell_layout());
            }
        }
    }
    // Pad the final row's trailing cells so the grid is exactly `n_rows * cols`.
    while cell_plots.len() < cell_count {
        cell_plots.push(Vec::new());
        cell_layouts.push(blank_cell_layout());
    }

    let scene = Figure::new(n_rows, cols)
        .with_plots(cell_plots)
        .with_layouts(cell_layouts)
        .with_title(&title)
        .with_cell_size(PLOT_WIDTH, PLOT_HEIGHT)
        .render();
    let pdf = PdfBackend::new().render_scene(&scene).map_err(|e| anyhow!("{e}"))?;
    std::fs::write(&path, pdf)
        .with_context(|| format!("writing duplication-spectrum plot to {}", path.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(
        library: &str,
        total: u64,
        duplicates: u64,
        window_duplicates: u64,
    ) -> DuplicateLadderRow {
        // Test rows assume the default 1M-template window.
        let window_total = 1_000_000;
        DuplicateLadderRow {
            sample: "s".to_string(),
            library: library.to_string(),
            category: "pairs",
            total,
            unique: total - duplicates,
            duplicates,
            frac_duplicates: duplicates as f64 / total as f64,
            window_total,
            window_unique: window_total - window_duplicates,
            window_duplicates,
            window_frac_duplicates: window_duplicates as f64 / window_total as f64,
        }
    }

    #[test]
    fn ladder_plot_path_appends_pdf_suffix() {
        assert_eq!(
            ladder_plot_path(Path::new("out/p")),
            Path::new("out/p.duplication-sampled.pdf")
        );
    }

    #[test]
    fn renders_a_valid_pdf_single_library() {
        let rows = vec![
            row("lib", 1_000_000, 50_000, 50_000),
            row("lib", 2_000_000, 90_000, 40_000),
            row("lib", 2_500_000, 95_000, 5_000), // partial final window
        ];
        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_duplicate_ladder_pdf(&rows, tmp.path()).unwrap();
        let bytes = std::fs::read(tmp.path()).unwrap();
        assert!(bytes.starts_with(b"%PDF"), "expected PDF magic bytes");
    }

    #[test]
    fn renders_one_pdf_for_multiple_libraries() {
        let rows = vec![
            row("libA", 1_000_000, 10_000, 10_000),
            row("libA", 2_000_000, 20_000, 10_000),
            row("libB", 1_000_000, 5_000, 5_000),
            row("libB", 2_000_000, 9_000, 4_000),
        ];
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.duplication-sampled.pdf");
        write_duplicate_ladder_pdf(&rows, &path).unwrap();
        // One file for the whole run (not one per library).
        assert!(path.exists());
        let entries: Vec<_> = std::fs::read_dir(dir.path()).unwrap().collect();
        assert_eq!(entries.len(), 1, "exactly one plot file");
        assert!(std::fs::read(&path).unwrap().starts_with(b"%PDF"));
    }

    #[test]
    fn empty_rows_writes_no_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.pdf");
        write_duplicate_ladder_pdf(&[], &path).unwrap();
        assert!(!path.exists(), "no file should be written for empty rows");
    }

    #[test]
    fn marginal_series_drops_partial_final_window() {
        let rows = [
            row("lib", 1_000_000, 10_000, 10_000),
            row("lib", 2_000_000, 20_000, 10_000),
            row("lib", 2_300_000, 22_000, 2_000), // 0.3M partial window
        ];
        let refs: Vec<&DuplicateLadderRow> = rows.iter().collect();
        let series = marginal_series(&refs);
        assert_eq!(series.len(), 2, "partial final point dropped");
        assert_eq!(series.last().unwrap().0, 2.0); // ends at the last full window
    }

    // ── Group-size histogram (η_k) ──────────────────────────────────────────

    fn hrow(library: &str, k: u32, n: u64) -> CountHistogramRow {
        CountHistogramRow {
            sample: "s".to_string(),
            library: library.to_string(),
            category: "pairs",
            n_observations: k,
            n_molecules: n,
        }
    }

    #[test]
    fn count_histogram_path_appends_suffix() {
        assert_eq!(
            count_histogram_plot_path(Path::new("out/p")),
            Path::new("out/p.duplication-spectrum.pdf")
        );
    }

    #[test]
    fn histogram_cutoff_caps_heavy_tails() {
        // Contiguous k=1..600 stays 100% populated -> raw 600 -> capped to 500.
        let rows: Vec<CountHistogramRow> = (1..=600).map(|k| hrow("lib", k, 10)).collect();
        let refs: Vec<&CountHistogramRow> = rows.iter().collect();
        assert_eq!(histogram_cutoff(&refs), 500.0);
    }

    #[test]
    fn histogram_cutoff_rounds_up_to_ten() {
        // Contiguous k=1..24 -> raw 24 -> rounds up to 30.
        let rows: Vec<CountHistogramRow> = (1..=24).map(|k| hrow("lib", k, 5)).collect();
        let refs: Vec<&CountHistogramRow> = rows.iter().collect();
        assert_eq!(histogram_cutoff(&refs), 30.0);
    }

    #[test]
    fn histogram_cutoff_handles_rows_starting_above_one() {
        // No singleton bucket (every signature seen >=2x), so rows start at k=2.
        // The (i+1)/k heuristic still holds — at k=20, 19 of [1..20] are
        // populated (95% >= 50%) — so raw stays 20 and rounds to 20.
        let rows: Vec<CountHistogramRow> = (2..=20).map(|k| hrow("lib", k, 5)).collect();
        let refs: Vec<&CountHistogramRow> = rows.iter().collect();
        assert_eq!(histogram_cutoff(&refs), 20.0);
    }

    #[test]
    fn histogram_cutoff_never_drops_below_smallest_k() {
        // Every signature seen far more than the cap: capping at 500 would fold
        // all points into the tail, so the cutoff must floor at the smallest k.
        let rows = [hrow("lib", 600, 3), hrow("lib", 900, 1)];
        let refs: Vec<&CountHistogramRow> = rows.iter().collect();
        assert_eq!(histogram_cutoff(&refs), 600.0);
    }

    #[test]
    fn count_histogram_renders_valid_pdf() {
        // A sparse tail point (k=600) trims the cutoff and folds into the label.
        let rows = [
            hrow("lib", 1, 1_000_000),
            hrow("lib", 2, 50_000),
            hrow("lib", 3, 2_000),
            hrow("lib", 600, 1),
        ];
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("out");
        write_count_histogram_pdfs(&rows, &prefix).unwrap();
        let pdf = std::fs::read(dir.path().join("out.duplication-spectrum.pdf")).unwrap();
        assert!(pdf.starts_with(b"%PDF"));
    }

    #[test]
    fn count_histogram_renders_when_all_signatures_exceed_cap() {
        // Regression: one signature seen 600x (> the 500 cutoff) with no
        // singletons used to fold the only point into the tail, leaving empty
        // series and an infinite log-y bound that panicked kuva. Must render.
        let rows = [hrow("lib", 600, 1)];
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("out");
        write_count_histogram_pdfs(&rows, &prefix).unwrap();
        let pdf = std::fs::read(dir.path().join("out.duplication-spectrum.pdf")).unwrap();
        assert!(pdf.starts_with(b"%PDF"));
    }

    #[test]
    fn count_histogram_multi_library_is_one_faceted_file() {
        let rows = vec![
            hrow("libA", 1, 100),
            hrow("libA", 2, 10),
            hrow("libB", 1, 50),
            hrow("libB", 2, 5),
        ];
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("out");
        write_count_histogram_pdfs(&rows, &prefix).unwrap();
        // One faceted file for the whole run — never one per library.
        let pdf = std::fs::read(dir.path().join("out.duplication-spectrum.pdf")).unwrap();
        assert!(pdf.starts_with(b"%PDF"));
        assert!(!dir.path().join("out.duplication-spectrum.libA.pdf").exists());
        assert!(!dir.path().join("out.duplication-spectrum.libB.pdf").exists());
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 1, "exactly one plot file");
    }

    #[test]
    fn count_histogram_single_library_unqualified_name() {
        let rows = vec![hrow("only", 1, 100), hrow("only", 2, 10)];
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("out");
        write_count_histogram_pdfs(&rows, &prefix).unwrap();
        assert!(dir.path().join("out.duplication-spectrum.pdf").exists());
    }

    #[test]
    fn pct_note_floors_tiny_and_formats_plain() {
        assert_eq!(pct_note(1e-8), "<0.001%");
        assert_eq!(pct_note(0.15), "15%");
        assert_eq!(pct_note(0.000012), "0.0012%");
    }

    #[test]
    fn plot_title_omits_dash_when_no_sample() {
        assert_eq!(plot_title("", "Duplication Spectrum (η_k)"), "Duplication Spectrum (η_k)");
        assert_eq!(
            plot_title("NA12878", "Duplication Spectrum (η_k)"),
            "NA12878 — Duplication Spectrum (η_k)"
        );
    }
}
