/** Methods / "How the islands work" page content. */

export function methodsHTML(): string {
  return `
  <article class="methods">
    <p class="eyebrow">Methods</p>
    <h1>How AlphaMissense <em>islands</em> work</h1>
    <p class="lead">
      An <b>AlphaMissense island</b> is a contiguous stretch of protein residues that
      AlphaMissense predicts to be consistently intolerant to missense change. Because
      AlphaMissense is derived from a protein language model and structural context — not
      from clinical databases — islands are an <b>orthogonal, non-circular</b> signal for
      the ACMG <b>PM1</b> criterion (“located in a mutational hot spot / critical
      functional domain”), unlike the traditional approach of counting ClinVar variants in
      a fixed window.
    </p>

    <h2>How an island is defined</h2>
    <ol class="steps">
      <li>
        <span class="step-k">Per-residue score</span>
        For each residue, take the <b>mean AlphaMissense pathogenicity</b> across all
        possible substitutions at that position.
      </li>
      <li>
        <span class="step-k">Disorder mask <span class="opt">optional</span></span>
        Residues with AlphaFold <b>pLDDT &lt; 50</b> are set to 0, so islands terminate in
        disordered regions where the structural context — and thus the prediction — is
        unreliable.
      </li>
      <li>
        <span class="step-k">Smoothing</span>
        A <b>median filter (window = 5 residues)</b> removes isolated spikes/dips, yielding
        longer, more robust segments.
      </li>
      <li>
        <span class="step-k">Segmentation</span>
        Keep contiguous runs where the smoothed score exceeds <b>0.564</b> (the
        AlphaMissense “ambiguous” boundary). A run breaks at the threshold or at any gap in
        residue numbering, and runs shorter than <b>35 residues</b> are discarded.
      </li>
    </ol>

    <figure class="svgfig">
      <canvas id="mfig" class="mfig" width="760" height="372"
        role="img" aria-label="Island creation: the per-residue mean of the AlphaMissense substitution heatmap gives a spiky signal; median-smoothing and a 0.564 threshold with a minimum length define islands"></canvas>
      <figcaption>The per-residue <b>mean</b> of the AlphaMissense substitution heatmap (bottom)
        is a spiky signal (grey); median-smoothing (blue) and the 0.564 threshold with a
        minimum-length filter turn sustained high-signal stretches into islands, while short
        spikes are discarded.</figcaption>
    </figure>

    <p class="src">Reference implementation: <code>computeAlphaMissenseIslands.py</code>
      (acmgscore-v2). Each island is emitted with its residue span, genomic coordinates and
      mean AM score — exactly the <code>islands</code> track you see in the viewer.</p>

    <h2>From islands to a PM1 prediction</h2>
    <p>
      At classification time a variant's residue is intersected with the island set. A
      variant that falls <b>inside an island</b> satisfies the AM-defined hot-spot condition;
      combining that with ClinVar pathogenic enrichment <i>within the island boundaries</i>
      gives a PM1 call that is both structurally grounded and clinically anchored. In this
      viewer the same islands are additionally <b>colored by their ClinVar overlap</b>
      (red = pathogenic-heavy, blue = benign-heavy), so you can see at a glance which
      predicted hot-spots are corroborated by real variants.
    </p>

    <h2>Grading PM1 by cumulative evidence</h2>
    <p>
      PM1 need not be all-or-nothing. Two independent lines of evidence can be combined and
      the criterion's <b>strength dialed</b> on the ACMG points scale (Supporting = 1,
      Moderate = 2, Strong = 4): the <b>AM island</b> (a structural/predictive hot-spot) and
      the <b>ClinVar pathogenic density inside that island</b> (empirical corroboration).
      The more the two agree — and the cleaner the island is of benign variation — the
      stronger the call.
    </p>
    <div class="grades">
      <table class="bench-table grades-table">
        <thead><tr><th>Evidence at the residue</th><th>PM1 strength</th></tr></thead>
        <tbody>
          <tr><td>Not in an AM island</td><td><span class="pill none">Not met</span></td></tr>
          <tr><td>In an AM island, little/no ClinVar corroboration</td><td><span class="pill sup">Supporting</span> <span class="pts">+1</span></td></tr>
          <tr><td>In an AM island <b>and</b> the island is pathogenic-enriched (several P/LP, no benign)</td><td><span class="pill mod">Moderate</span> <span class="pts">+2</span></td></tr>
          <tr><td>In a dense, well-established hot-spot island (many P/LP, zero benign)</td><td><span class="pill str">Strong</span> <span class="pts">+4</span></td></tr>
          <tr><td>Island contradicted by benign variation (B/LB-heavy)</td><td><span class="pill none">Downweight / not met</span></td></tr>
        </tbody>
      </table>
    </div>
    <p class="src">
      This is the framework Archipelago's data supports — the island track supplies the
      hot-spot boundary, and each island's red/blue ClinVar colouring supplies the
      corroboration. The exact thresholds and the maximum allowed strength are ultimately
      gene- and VCEP-specific (per ClinGen SVI recommendations), so a lab tunes the cut-offs;
      the viewer gives the two signals needed to make that call.
    </p>

    <h2>Benchmark — PM1</h2>
    <p>
      On an independent ACMG benchmark (support = 1,270 variants), the AlphaMissense-island
      PM1 is dramatically <b>more precise</b> than ClinVar-window approaches — ~3.6× the
      precision of InterVar and ~2× Franklin — and achieves the <b>best overall F1</b>.
    </p>
    <div class="bench">
      <table class="bench-table">
        <thead>
          <tr><th>Tool</th><th>Precision</th><th>Recall</th><th>F1</th></tr>
        </thead>
        <tbody>
          <tr class="win"><td>SeqOne <span class="tag">AM islands</span></td><td>0.70</td><td>0.50</td><td>0.59</td></tr>
          <tr><td>Franklin</td><td>0.33</td><td>0.86</td><td>0.47</td></tr>
          <tr><td>InterVar</td><td>0.19</td><td>0.84</td><td>0.31</td></tr>
        </tbody>
      </table>
      <figure class="bench-fig">
        <img src="/assets/benchmark_pm1.png" alt="PM1 precision, recall and F1 for SeqOne, InterVar and Franklin" />
        <figcaption>PM1 — mutational hotspot / functional domain. Precision / Recall / F1.</figcaption>
      </figure>
    </div>
    <p class="src">
      The high recall of window-based tools comes at a steep precision cost (they flag broad
      regions as hot-spots); the island approach trades a little recall for far fewer false
      positives — the right trade-off for an ACMG <i>moderate</i> criterion.
    </p>
    <p class="bench-cta"><a href="?page=benchmark" data-benchmark>Explore every PM1 call (TP / FP / FN) →</a></p>

    <p class="src">Cancer-hotspot track: cancerhotspots.org (Chang et al., 2018).</p>

    <p class="back"><a href="?page=methods" data-home>← Back to search</a></p>
  </article>`;
}

// ---- realistic island-construction figure (canvas) ----
type RGB = [number, number, number];
const css = (v: string) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
function hexToRgb(h: string): RGB {
  h = h.replace("#", "");
  if (h.length === 3) h = h.split("").map((c) => c + c).join("");
  return [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16)];
}
const mix = (a: RGB, b: RGB, t: number): RGB =>
  [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
const rgb = (c: RGB) => `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`;

/** Deterministic PRNG so the schematic is stable across renders. */
function mulberry32(seed: number): () => number {
  return () => {
    seed |= 0; seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function median5(a: number[]): number[] {
  const out = a.slice();
  for (let i = 0; i < a.length; i++) {
    const lo = Math.max(0, i - 2), hi = Math.min(a.length - 1, i + 2);
    out[i] = a.slice(lo, hi + 1).sort((x, y) => x - y)[Math.floor((hi - lo) / 2)];
  }
  return out;
}

const N = 104, ROWS = 20, THRESH = 0.564, MIN_RUN = 16;

function buildFigureData() {
  const rnd = mulberry32(20240714);
  // [start, end, height, rampWidth] — flat-topped plateaus (two wide, one narrow spike)
  const bumps: [number, number, number, number][] = [
    [12, 49, 0.64, 6],   // wide → island
    [60, 69, 0.42, 5],   // narrow spike → rejected (never reaches its flat top)
    [76, 101, 0.64, 6],  // wide → island
  ];
  const base = 0.27;
  const raw: number[] = [];
  for (let i = 0; i < N; i++) {
    let v = base;
    for (const [s, e, add, ramp] of bumps) {
      if (i >= s && i <= e) v += add * Math.min(1, Math.min(i - s, e - i) / ramp);
    }
    v += (rnd() - 0.5) * 0.13;                        // spiky per-residue noise
    if (rnd() < 0.09) v -= 0.12;                       // occasional tolerant residue
    raw.push(Math.max(0.02, Math.min(0.98, v)));
  }
  const smooth = median5(raw);
  // islands = smoothed runs above threshold; long ones pass, short ones rejected
  const runs: { s: number; e: number; ok: boolean }[] = [];
  let i = 0;
  while (i < N) {
    if (smooth[i] > THRESH) {
      let j = i; while (j < N && smooth[j] > THRESH) j++;
      runs.push({ s: i, e: j - 1, ok: j - i >= MIN_RUN }); i = j;
    } else i++;
  }
  // heatmap cells: column mean tracks `raw`; per-row bias gives realistic texture
  const rowBias = Array.from({ length: ROWS }, (_, r) => ((r / (ROWS - 1)) - 0.5) * 0.34);
  const grid: number[][] = [];
  const refRow: number[] = [];
  for (let c = 0; c < N; c++) {
    const col: number[] = [];
    for (let r = 0; r < ROWS; r++) {
      col.push(Math.max(0, Math.min(1, raw[c] + rowBias[r] + (rnd() - 0.5) * 0.4)));
    }
    grid.push(col);
    refRow.push(Math.floor(rnd() * ROWS));
  }
  return { raw, smooth, runs, grid, refRow };
}

let FIG: ReturnType<typeof buildFigureData> | null = null;
let figObserver: MutationObserver | null = null;

export function drawMethodFigure(): void {
  const cv = document.getElementById("mfig") as HTMLCanvasElement | null;
  if (!cv) { figObserver?.disconnect(); figObserver = null; return; }
  if (!FIG) FIG = buildFigureData();
  const { raw, smooth, runs, grid, refRow } = FIG;

  const W = 760, H = 372, dpr = Math.min(window.devicePixelRatio || 1, 2);
  cv.width = W * dpr; cv.height = H * dpr;
  const ctx = cv.getContext("2d")!;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, W, H);

  const padL = 52, padR = 16, plotW = W - padL - padR, colW = plotW / N;
  const sy0 = 16, sy1 = 176;                          // signal panel (value 1 → 0)
  const hy0 = 232, rowH = 6, hy1 = hy0 + ROWS * rowH; // heatmap panel
  const Y = (v: number) => sy1 - v * (sy1 - sy0);
  const X = (i: number) => padL + i * colW;

  const stops = ["--am0", "--am1", "--am2", "--am3", "--am4"].map((v) => hexToRgb(css(v)));
  const ramp = (v: number): RGB => {
    const x = Math.max(0, Math.min(1, v)) * 4, k = Math.min(3, Math.floor(x));
    return mix(stops[k], stops[k + 1], x - k);
  };
  const ink = css("--ink"), ink2 = css("--ink2"), muted = css("--muted");
  const border = css("--border"), islp = css("--isl-path"), accent = css("--accent");
  const refCol = css("--reference");

  // heatmap
  for (let c = 0; c < N; c++) {
    for (let r = 0; r < ROWS; r++) {
      ctx.fillStyle = r === refRow[c] ? refCol : rgb(ramp(grid[c][r]));
      ctx.fillRect(X(c), hy0 + r * rowH, Math.ceil(colW) - 0.3, rowH - 0.3);
    }
  }

  // island bands (both panels) + rejected run
  for (const run of runs) {
    const x = X(run.s), w = (run.e - run.s + 1) * colW;
    if (run.ok) {
      ctx.fillStyle = islp; ctx.globalAlpha = 0.15;
      ctx.fillRect(x, sy0, w, sy1 - sy0);
      ctx.globalAlpha = 0.1; ctx.fillRect(x, hy0, w, hy1 - hy0); ctx.globalAlpha = 1;
    } else {
      ctx.strokeStyle = muted; ctx.globalAlpha = 0.9; ctx.setLineDash([3, 3]); ctx.lineWidth = 1;
      ctx.strokeRect(x + 0.5, sy0 + 0.5, w - 1, sy1 - sy0 - 1); ctx.setLineDash([]); ctx.globalAlpha = 1;
    }
  }

  // baseline + threshold
  ctx.strokeStyle = border; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(padL, sy1 + 0.5); ctx.lineTo(W - padR, sy1 + 0.5); ctx.stroke();
  ctx.strokeStyle = muted; ctx.setLineDash([5, 4]);
  ctx.beginPath(); ctx.moveTo(padL, Y(THRESH)); ctx.lineTo(W - padR, Y(THRESH)); ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = muted; ctx.font = "11px ui-monospace,monospace"; ctx.textBaseline = "alphabetic";
  ctx.textAlign = "left"; ctx.fillText("0.564", padL + 2, Y(THRESH) - 4);
  ctx.textAlign = "right";
  ctx.fillText("1.0", padL - 6, sy0 + 8); ctx.fillText("0", padL - 6, sy1);

  // raw (spiky) + smoothed lines
  const line = (arr: number[], color: string, w: number, alpha = 1) => {
    ctx.strokeStyle = color; ctx.globalAlpha = alpha; ctx.lineWidth = w;
    ctx.lineJoin = "round"; ctx.beginPath();
    arr.forEach((v, i) => { const x = X(i) + colW / 2, y = Y(v); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); });
    ctx.stroke(); ctx.globalAlpha = 1;
  };
  line(raw, muted, 1, 0.6);
  line(smooth, accent, 2.4);

  // labels
  ctx.textAlign = "center"; ctx.font = "600 11.5px ui-sans-serif,system-ui";
  for (const run of runs) {
    const cx = X((run.s + run.e) / 2) + colW / 2;
    ctx.fillStyle = run.ok ? islp : muted;
    ctx.fillText(run.ok ? "island ✓ ≥35 aa" : "✗ too short", cx, sy0 + 12);
  }
  ctx.fillStyle = ink; ctx.textAlign = "left"; ctx.font = "600 11px ui-sans-serif,system-ui";
  ctx.fillText("Mean AlphaMissense per residue", padL, sy0 - 4 + 0);
  // connector: heatmap → column mean
  ctx.fillStyle = ink2; ctx.font = "11px ui-sans-serif,system-ui";
  ctx.fillText("AlphaMissense heatmap — 20 substitutions per residue (column mean ↑ drives the signal)", padL, hy1 + 16);
  ctx.fillStyle = muted; ctx.textAlign = "right";
  ctx.fillText("residue →", W - padR, hy1 + 16);

  // redraw on theme change while the figure is on screen
  if (!figObserver) {
    figObserver = new MutationObserver(() => drawMethodFigure());
    figObserver.observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  }
}
