/** Canvas track-stack viewer: AM heatmap + mean + islands + ClinVar + domains. */

export interface Island { s: number; e: number; m: number; plp?: number; blb?: number; }
export interface Domain { name: string; s: number; e: number; }
export interface CVar { p: number; ref: string; alt: string; sig: string; n: number; st?: number; }
export interface Bundle {
  gene: string; uniprot: string; ensembl: string; refseq: string;
  length: number; aa_order: string;
  mean: (number | null)[]; ref: string; grid: string;
  islands: Island[]; domains: Domain[];
  clinvar: CVar[]; clinvar_benign: CVar[];
  /** ClinVar variant counts by consequence: {missense|lof|other: [P/LP, B/LB]}. */
  cv_counts?: Record<"missense" | "lof" | "other", [number, number]>;
  /** cancerhotspots.org residue hotspots (present only for genes with entries). */
  hotspots?: { p: number; r: string; q?: number; n?: number; onc: number }[];
  /** gnomAD regional missense constraint: observed/expected per sub-region. */
  rmc?: { s: number; e: number; oe: number; obs?: number; exp?: number; p?: number; whole?: number }[];
  /** exon-start residues (1-based) for multi-exon transcripts. */
  exons?: number[];
}

const ICON: Record<string, string> = {
  alert: `<svg class="kpi-i" viewBox="0 0 24 24"><path d="M12 3.5 21 19.5H3z"/><path d="M12 9.5v4.5"/><circle cx="12" cy="16.8" r=".3"/></svg>`,
  shield: `<svg class="kpi-i" viewBox="0 0 24 24"><path d="M12 3l7 3v5c0 4.4-3 7.5-7 9-4-1.5-7-4.6-7-9V6z"/><path d="M8.7 12l2.2 2.2L15.4 10"/></svg>`,
  pulse: `<svg class="kpi-i" viewBox="0 0 24 24"><path d="M3 12h4l3 7 4-15 3 8h4"/></svg>`,
  island: `<svg class="kpi-i" viewBox="0 0 24 24"><path d="M7.5 13.5C9 8 10.5 6 12 6s3 2 4.5 7.5"/><path d="M3 16.5c1.8-1.8 3.6-1.8 5.4 0s3.6 1.8 5.4 0 3.6-1.8 5.4 0"/><path d="M3 20c1.8-1.8 3.6-1.8 5.4 0s3.6 1.8 5.4 0 3.6-1.8 5.4 0"/></svg>`,
  layers: `<svg class="kpi-i" viewBox="0 0 24 24"><path d="M12 3l8.5 4.5L12 12 3.5 7.5z"/><path d="M4 12.5l8 4.2 8-4.2"/></svg>`,
  expand: `<svg class="kpi-x" viewBox="0 0 24 24"><path d="M9 4H4v5"/><path d="M15 20h5v-5"/><path d="M4 4l6 6"/><path d="M20 20l-6-6"/></svg>`,
  target: `<svg class="kpi-i" viewBox="0 0 24 24"><circle cx="12" cy="12" r="8"/><circle cx="12" cy="12" r="3.4"/></svg>`,
  constraint: `<svg class="kpi-i" viewBox="0 0 24 24"><path d="M9 5v4H5"/><path d="M15 5v4h4"/><path d="M9 19v-4H5"/><path d="M15 19v-4h4"/></svg>`,
};

const B_REF = 251, B_MISSING = 252, Q_MAX = 250;
const CUTOFF_B = 0.34, CUTOFF_P = 0.564;
const AXIS = 22, DOM = 34, ISL = 26, RMC = 24, HOTS = 26, CVH = 66, SCORE = 76, ROWH = 15;

const css = (v: string) =>
  getComputedStyle(document.documentElement).getPropertyValue(v).trim();

function hexToRgb(h: string): [number, number, number] {
  h = h.replace("#", "");
  if (h.length === 3) h = h.split("").map((c) => c + c).join("");
  return [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16)];
}
function mix(a: [number, number, number], b: [number, number, number], t: number): string {
  return `rgb(${Math.round(a[0] + (b[0] - a[0]) * t)},${Math.round(a[1] + (b[1] - a[1]) * t)},${Math.round(a[2] + (b[2] - a[2]) * t)})`;
}
function rampLUT(): string[] {
  const stops = [css("--am0"), css("--am1"), css("--am2"), css("--am3"), css("--am4")].map(hexToRgb);
  const lut: string[] = [];
  for (let i = 0; i <= Q_MAX; i++) {
    const x = (i / Q_MAX) * (stops.length - 1);
    const k = Math.min(stops.length - 2, Math.floor(x));
    lut[i] = mix(stops[k], stops[k + 1], x - k);
  }
  return lut;
}
function divergeColor(frac: number): string {
  // frac in [-1,1]: +1 pathogenic (red), 0 mixed (neutral), -1 benign (blue)
  const path = hexToRgb(css("--isl-path"));
  const mid = hexToRgb(css("--isl-mixed"));
  const ben = hexToRgb(css("--isl-benign"));
  return frac >= 0 ? mix(mid, path, frac) : mix(mid, ben, -frac);
}

function decodeGrid(grid: string, n: number, naa: number): Uint8Array {
  const bin = atob(grid);
  const out = new Uint8Array(n * naa);
  for (let i = 0; i < out.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}

function cvModalHTML(b: Bundle): string {
  const c = b.cv_counts;
  if (!c) return "";
  const rows: [string, [number, number]][] = [
    ["Missense", c.missense], ["Truncating / LoF", c.lof], ["Other", c.other],
  ];
  const tot = rows.reduce((s, [, v]) => [s[0] + v[0], s[1] + v[1]], [0, 0]);
  if (tot[0] + tot[1] === 0) return "";
  return `<div class="modal-backdrop" id="cv-modal" hidden>
    <div class="modal" role="dialog" aria-modal="true" aria-label="ClinVar variants by consequence">
      <div class="modal-head">
        <h3>ClinVar variants by consequence</h3>
        <button class="modal-x" data-cvclose aria-label="Close">×</button>
      </div>
      <p class="modal-note">${b.gene} · ${b.refseq || "transcript"} — Pathogenic/Likely-path.
        and Benign/Likely-benign ClinVar records, grouped by molecular consequence.
        The heatmap and islands use <b>missense only</b>.</p>
      <table class="cvtable">
        <thead><tr><th>Consequence</th><th>P / LP</th><th>B / LB</th></tr></thead>
        <tbody>
          ${rows.map(([k, v]) => `<tr><td>${k}</td><td class="p">${v[0]}</td><td class="b">${v[1]}</td></tr>`).join("")}
          <tr class="tot"><td>Total</td><td class="p">${tot[0]}</td><td class="b">${tot[1]}</td></tr>
        </tbody>
      </table>
    </div>
  </div>`;
}

const CVCOL: Record<string, string> = {};
const CVNAME: Record<string, string> = {
  P: "Pathogenic", PLP: "Pathogenic / Likely-path.", LP: "Likely pathogenic",
  B: "Benign", BLB: "Benign / Likely-benign", LB: "Likely benign",
};

export function renderViewer(root: HTMLElement, b: Bundle): void {
  const AA = b.aa_order;
  const N = b.length;
  const HEAT = ROWH * AA.length;
  const grid = decodeGrid(b.grid, N, AA.length);
  const cvMap = new Map<number, CVar>(); b.clinvar.forEach((v) => cvMap.set(v.p, v));
  const cvbMap = new Map<number, CVar>(); b.clinvar_benign.forEach((v) => cvbMap.set(v.p, v));

  // Headline = ClinVar missense *variant* count (matches the modal breakdown's
  // Missense row); the track/lollipops are the deduped *residues* shown as sublabel.
  const ms = b.cv_counts?.missense;
  const plpVariants = ms ? ms[0] : b.clinvar.reduce((s, v) => s + v.n, 0);
  const blbVariants = ms ? ms[1] : b.clinvar_benign.reduce((s, v) => s + v.n, 0);
  const scored = b.mean.filter((m) => m != null) as number[];
  const amPathoPct = scored.length
    ? Math.round((100 * scored.filter((m) => m >= CUTOFF_P).length) / scored.length)
    : 0;
  const islandResidues = b.islands.reduce((s, i) => s + (i.e - i.s + 1), 0);
  const hasBd = !!(b.cv_counts && Object.values(b.cv_counts).some((v) => v[0] + v[1] > 0));
  const hasHot = !!(b.hotspots && b.hotspots.length);
  const hasRmc = !!(b.rmc && b.rmc.length);
  const minOE = hasRmc ? Math.min(...b.rmc!.map((r) => r.oe)) : null;

  root.innerHTML = `
    <div class="vhead">
      <h2><em>${b.gene}</em> AlphaMissense saturation map
        <span class="mono" style="font-size:.6em;color:var(--muted)">${b.uniprot}</span></h2>
      <div class="vmeta">
        <span><b>UniProt</b> <span class="mono">${b.uniprot}</span></span>
        <span><b>Transcript</b> <span class="mono">${b.ensembl}${b.refseq ? " · " + b.refseq : ""}</span></span>
        <span><b>Length</b> <span class="mono">${N}</span> aa</span>
      </div>
      <div class="kpis">
        <${hasBd ? "button" : "div"} class="kpi${hasBd ? " kpi-btn" : ""}"${hasBd ? ' data-cvmodal title="View ClinVar variants by consequence (missense · LoF · other)"' : ""}>
          <div class="kpi-head">${ICON.alert}<span>ClinVar P/LP</span>${hasBd ? `<span class="kpi-more">${ICON.expand}</span>` : ""}</div>
          <div class="kpi-n">${plpVariants}</div>
          <div class="kpi-sub">missense · ${b.clinvar.length} residue${b.clinvar.length === 1 ? "" : "s"}</div>
        </${hasBd ? "button" : "div"}>
        <${hasBd ? "button" : "div"} class="kpi${hasBd ? " kpi-btn" : ""}"${hasBd ? ' data-cvmodal title="View ClinVar variants by consequence (missense · LoF · other)"' : ""}>
          <div class="kpi-head">${ICON.shield}<span>ClinVar B/LB</span>${hasBd ? `<span class="kpi-more">${ICON.expand}</span>` : ""}</div>
          <div class="kpi-n">${blbVariants}</div>
          <div class="kpi-sub">missense · ${b.clinvar_benign.length} residue${b.clinvar_benign.length === 1 ? "" : "s"}</div>
        </${hasBd ? "button" : "div"}>
        <div class="kpi">
          <div class="kpi-head">${ICON.pulse}<span>AM-pathogenic</span></div>
          <div class="kpi-n">${amPathoPct}<span class="pct">%</span></div>
          <div class="kpi-sub">of residues (mean ≥ ${CUTOFF_P})</div>
        </div>
        <div class="kpi">
          <div class="kpi-head">${ICON.island}<span>AM islands</span></div>
          <div class="kpi-n">${b.islands.length}</div>
          <div class="kpi-sub">${islandResidues} residues</div>
        </div>
        <div class="kpi">
          <div class="kpi-head">${ICON.layers}<span>Domains</span></div>
          <div class="kpi-n">${b.domains.length}</div>
          <div class="kpi-sub">UniProt features</div>
        </div>
        ${hasHot ? `<div class="kpi">
          <div class="kpi-head">${ICON.target}<span>Cancer hotspots</span></div>
          <div class="kpi-n">${b.hotspots!.length}</div>
          <div class="kpi-sub">cancerhotspots.org</div>
        </div>` : ""}
        ${hasRmc ? `<div class="kpi">
          <div class="kpi-head">${ICON.constraint}<span>Missense o/e</span></div>
          <div class="kpi-n">${minOE!.toFixed(2)}</div>
          <div class="kpi-sub">${b.rmc![0].whole ? "gnomAD gene-level" : "gnomAD region (min)"}</div>
        </div>` : ""}
      </div>
    </div>
    <div class="card">
      <div class="toolbar">
        <div class="legend">
          <span class="key"><i class="sw grad"></i>AM score 0→1</span>
          <span class="key"><i class="sw reference"></i>Reference</span>
          <span class="sep"></span>
          <span class="key"><i class="sw dot cvp"></i>ClinVar P/LP ↑</span>
          <span class="key"><i class="sw dot cvblb"></i>B/LB ↓</span>
          <span class="sep"></span>
          <span class="key"><i class="sw islp"></i>island: P-rich</span>
          <span class="key"><i class="sw islb"></i>B-rich</span>
          <span class="key"><i class="sw isln"></i>no ClinVar</span>
          <span class="key"><i class="sw domain"></i>Domain</span>
          ${hasHot ? `<span class="key"><i class="sw hotspot"></i>Cancer hotspot</span>` : ""}
          ${hasRmc ? `<span class="key"><i class="sw rmc"></i>Missense-constrained (low o/e)</span>` : ""}
        </div>
        <div class="zoom">
          <button id="v-png" title="Download a PNG of the whole plot">Save PNG</button>
          <button id="v-fit">Fit</button>
          <label for="v-z">Zoom</label>
          <input id="v-z" type="range" min="0.35" max="60" step="0.05" />
        </div>
      </div>
      <div class="plot">
        <div class="gutter" id="v-gutter"></div>
        <div class="scroll" id="v-scroll">
          <div class="stack" id="v-stack">
            <canvas id="v-axt"></canvas>
            <canvas id="v-dom"></canvas>
            <canvas id="v-isl"></canvas>
            ${hasRmc ? `<canvas id="v-rmc"></canvas>` : ""}
            ${hasHot ? `<canvas id="v-hot"></canvas>` : ""}
            <canvas id="v-cv"></canvas>
            <canvas id="v-score"></canvas>
            <canvas id="v-heat"></canvas>
            <canvas id="v-axb"></canvas>
            <div class="xhair" id="v-xhair"></div>
            <div class="zoomsel" id="v-zoomsel"></div>
          </div>
        </div>
      </div>
      <div class="hint">Heatmap cell = AlphaMissense pathogenicity for that substitution.
        Islands are colored by ClinVar overlap (red = pathogenic-heavy, blue = benign-heavy,
        amber = none). ClinVar lollipops: P/LP up, B/LB down — head size ∝ #records, pin height ∝ review stars (0–4★). Hover for details; drag across the plot to zoom to a region, Fit to reset.</div>
    </div>
    <p class="vfoot">AlphaMissense (DeepMind, CC BY 4.0) · ClinVar P/LP + B/LB missense on
      <code>${b.refseq || "RefSeq"}</code> · UniProt domains · AM islands (score&gt;0.564, ≥35 aa).</p>
    <div class="tip" id="v-tip"></div>
    ${cvModalHTML(b)}`;

  const modal = root.querySelector<HTMLElement>("#cv-modal");
  if (modal) {
    root.querySelectorAll("[data-cvmodal]").forEach((el) =>
      el.addEventListener("click", () => { modal.hidden = false; }));
    modal.addEventListener("click", (e) => {
      if ((e.target as HTMLElement).closest("[data-cvclose]") || e.target === modal) modal.hidden = true;
    });
    root.addEventListener("keydown", (e) => { if ((e as KeyboardEvent).key === "Escape") modal.hidden = true; });
  }

  const $ = (id: string) => root.querySelector<HTMLElement>("#" + id)!;
  const scroll = $("v-scroll") as HTMLDivElement;
  const stack = $("v-stack");
  const tip = $("v-tip");
  const xhair = $("v-xhair");
  const gutter = $("v-gutter");
  const heat = $("v-heat") as HTMLCanvasElement;
  const zoomsel = $("v-zoomsel");
  const MAXW = 60;
  let dragX: number | null = null;
  const plotX = (e: MouseEvent) => e.clientX - scroll.getBoundingClientRect().left + scroll.scrollLeft;

  const LANES = [
    { id: "v-axt", h: AXIS, label: b.exons && b.exons.length ? "Exons" : "" },
    { id: "v-dom", h: DOM, label: "Domains" },
    { id: "v-isl", h: ISL, label: "AM islands" },
    ...(hasRmc ? [{ id: "v-rmc", h: RMC, label: "Missense o/e" }] : []),
    ...(hasHot ? [{ id: "v-hot", h: HOTS, label: "Cancer hotspots" }] : []),
    { id: "v-cv", h: CVH, label: "ClinVar" },
    { id: "v-score", h: SCORE, label: "Mean AM" },
    { id: "v-heat", h: HEAT, label: "" },
    { id: "v-axb", h: AXIS, label: "Residue" },
  ];
  const laneCV = LANES.find((l) => l.id === "v-cv")!;
  // ClinVar lollipops feel cramped when zoomed in; scale the whole lane (height +
  // lollipop geometry) with the zoom level, up to CV_ZMAX×.
  const CV_ZMAX = 2.4;
  const cvZoom = () => Math.min(CV_ZMAX, Math.max(1, colW / 4));
  const dpr = Math.max(1, window.devicePixelRatio || 1);
  // Browsers cap canvas backing-store dimensions (~32k px in Firefox, higher in Chrome).
  // A long protein zoomed in makes contentW × dpr exceed that and the canvas fails to
  // allocate (blank/crash). Clamp the backing width; the CSS width stays contentW so
  // scroll/hover geometry is unchanged — only the pixel buffer downscales past the limit.
  const MAX_CANVAS_DIM = 30000;
  let colW = 1;
  // Virtualised rendering: canvases are only viewport-wide and redrawn on scroll, so the
  // backing store never approaches the browser's max-canvas-dimension no matter how far a
  // long protein is zoomed. `viewW` = visible plot width, `scrollX` = horizontal offset,
  // both in CSS px; draw functions work in full-protein coordinates and are translated.
  let viewW = 1;
  let scrollX = 0;
  const visRange = (): [number, number] => [
    Math.max(0, Math.floor(scrollX / colW) - 1),
    Math.min(N - 1, Math.ceil((scrollX + viewW) / colW) + 1),
  ];
  let LUT = rampLUT();
  const refreshColors = () => {
    LUT = rampLUT();
    Object.assign(CVCOL, {
      P: css("--cv-p"), PLP: css("--cv-plp"), LP: css("--cv-lp"),
      B: css("--cv-b"), BLB: css("--cv-blb"), LB: css("--cv-lb"),
    });
  };
  refreshColors();

  const contentW = () => Math.max(1, Math.round(N * colW));
  const X = (p: number) => (p - 1) * colW;
  function ctxFor(id: string, h: number): CanvasRenderingContext2D {
    const cv = $(id) as HTMLCanvasElement;
    // Canvas spans only the visible width; sx = dpr on screen (viewW small), and clamps in
    // the full-width export path where viewW is set to contentW.
    const sx = Math.min(dpr, MAX_CANVAS_DIM / Math.max(1, viewW));
    cv.style.width = viewW + "px"; cv.style.height = h + "px";
    cv.width = Math.max(1, Math.round(viewW * sx)); cv.height = Math.round(h * dpr);
    const ctx = cv.getContext("2d")!;
    // draw in full-protein coordinates (X(p) = (p-1)*colW); shift the visible window into view
    ctx.setTransform(sx, 0, 0, dpr, -Math.round(scrollX * sx), 0);
    return ctx;
  }

  function buildGutter(): void {
    gutter.innerHTML = "";
    for (const L of LANES) {
      const d = document.createElement("div");
      d.style.height = L.h + "px";
      if (L.id === "v-heat") {
        d.className = "grows";
        for (let i = 0; i < AA.length; i++) {
          const r = document.createElement("div");
          r.className = "grow"; r.style.top = i * ROWH + ROWH / 2 + "px";
          r.innerHTML = `<span style="letter-spacing:.5px">${AA[i]}</span>`;
          d.appendChild(r);
        }
      } else {
        d.className = "glane";
        d.innerHTML = L.id === "v-score"
          ? `<span>Mean AM<br><span style="color:var(--muted);font-size:10px">1.0 – 0.0</span></span>`
          : L.id === "v-cv"
            ? `<span>ClinVar<br><span style="color:var(--muted);font-size:10px">P/LP ↑ · B/LB ↓<br>height ∝ review ★</span></span>`
            : `<span>${L.label}</span>`;
      }
      gutter.appendChild(d);
    }
  }

  function drawHeat(): void {
    const ctx = ctxFor("v-heat", HEAT);
    const gap = colW >= 4 ? 1 : 0;
    const rw = Math.max(colW - gap, colW < 1 ? colW : 0.6);
    const refCol = css("--reference");
    const [lo, hi] = visRange();
    for (let pi = lo; pi <= hi; pi++) {
      const base = pi * AA.length; const x = pi * colW;
      for (let r = 0; r < AA.length; r++) {
        const val = grid[base + r];
        if (val === B_MISSING) continue;
        ctx.fillStyle = val === B_REF ? refCol : LUT[val];
        ctx.fillRect(x, r * ROWH, rw, ROWH - gap);
      }
    }
  }

  const scoreY = (v: number) => SCORE - 4 - v * (SCORE - 8);
  function drawScore(): void {
    const ctx = ctxFor("v-score", SCORE);
    ctx.strokeStyle = css("--grid"); ctx.lineWidth = 1;
    for (const v of [CUTOFF_B, CUTOFF_P]) {
      ctx.beginPath(); ctx.moveTo(0, scoreY(v) + 0.5); ctx.lineTo(contentW(), scoreY(v) + 0.5); ctx.stroke();
    }
    const [lo, hi] = visRange();
    const area = new Path2D(); let st = false, lastX = 0;
    for (let pi = lo; pi <= hi; pi++) {
      const m = b.mean[pi]; const x = pi * colW + colW / 2;
      if (m != null) { const y = scoreY(m); if (!st) { area.moveTo(x, scoreY(0)); area.lineTo(x, y); st = true; } else area.lineTo(x, y); lastX = x; }
    }
    if (st) { area.lineTo(lastX, scoreY(0)); area.closePath(); ctx.fillStyle = css("--score-fill"); ctx.fill(area); }
    ctx.beginPath(); st = false;
    for (let pi = lo; pi <= hi; pi++) {
      const m = b.mean[pi]; const x = pi * colW + colW / 2;
      if (m != null) { const y = scoreY(m); if (!st) { ctx.moveTo(x, y); st = true; } else ctx.lineTo(x, y); }
    }
    ctx.strokeStyle = css("--score-stroke"); ctx.lineWidth = 1.2; ctx.stroke();
  }

  function lolli(ctx: CanvasRenderingContext2D, x: number, mid: number, r: number, up: boolean, color: string, ext = 0, scale = 1): void {
    const y = up ? mid - 6 * scale - r - ext : mid + 6 * scale + r + ext;
    ctx.strokeStyle = css("--muted"); ctx.globalAlpha = 0.5; ctx.lineWidth = Math.max(1, scale);
    ctx.beginPath(); ctx.moveTo(x, mid); ctx.lineTo(x, y); ctx.stroke(); ctx.globalAlpha = 1;
    ctx.fillStyle = color; ctx.beginPath(); ctx.arc(x, y, r, 0, 6.2832); ctx.fill();
    ctx.strokeStyle = css("--panel"); ctx.lineWidth = 1; ctx.stroke();
  }
  const STAR_EXT = 4;  // stem length added per ClinVar gold star (0..4)
  function drawClinvar(): void {
    const H = laneCV.h, zf = cvZoom();
    const ctx = ctxFor("v-cv", H); const mid = H / 2;
    ctx.strokeStyle = css("--grid"); ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(0, mid + 0.5); ctx.lineTo(contentW(), mid + 0.5); ctx.stroke();
    const rad = (n: number) => (2.2 + Math.min(n - 1, 4) * 0.5) * zf;
    const pr: Record<string, number> = { P: 0, PLP: 1, LP: 2 };
    [...b.clinvar].sort((a, c) => pr[a.sig] - pr[c.sig]).reverse().forEach((v) =>
      lolli(ctx, X(v.p) + Math.max(colW / 2, 0.5), mid, rad(v.n), true, CVCOL[v.sig], (v.st || 0) * STAR_EXT * zf, zf));
    const br: Record<string, number> = { B: 0, BLB: 1, LB: 2 };
    [...b.clinvar_benign].sort((a, c) => br[a.sig] - br[c.sig]).reverse().forEach((v) =>
      lolli(ctx, X(v.p) + Math.max(colW / 2, 0.5), mid, rad(v.n), false, CVCOL[v.sig], (v.st || 0) * STAR_EXT * zf, zf));
  }
  function islandColor(i: Island): string {
    const plp = i.plp || 0, blb = i.blb || 0;
    if (plp + blb === 0) return css("--isl-none");
    return divergeColor((plp - blb) / (plp + blb));
  }
  function drawHotspots(): void {
    const ctx = ctxFor("v-hot", HOTS); const y = HOTS / 2;
    for (const h of b.hotspots!) {
      const x = X(h.p) + Math.max(colW / 2, 0.5);
      const r = 3 + Math.min((h.n || 1) / 20, 3);
      ctx.beginPath();
      ctx.moveTo(x, y - r); ctx.lineTo(x + r, y); ctx.lineTo(x, y + r); ctx.lineTo(x - r, y); ctx.closePath();
      if (h.onc) { ctx.fillStyle = css("--hotspot"); ctx.fill(); }
      else { ctx.strokeStyle = css("--hotspot"); ctx.lineWidth = 1.5; ctx.stroke(); }
    }
  }
  function drawRmc(): void {
    const ctx = ctxFor("v-rmc", RMC);
    const col = hexToRgb(css("--rmc"));
    for (const r of b.rmc!) {
      const x = X(r.s), w = (r.e - r.s + 1) * colW;
      const constr = Math.max(0, Math.min(1, (0.85 - r.oe) / 0.85));
      // gene-level (no sub-regional constraint) is drawn subtler than a real sub-region
      const a = (r.whole ? 0.1 + 0.45 * constr : 0.15 + 0.7 * constr);
      ctx.fillStyle = `rgba(${col[0]},${col[1]},${col[2]},${a.toFixed(3)})`;
      ctx.fillRect(x, 4, Math.max(w, 1), RMC - 8);
      ctx.strokeStyle = "rgba(0,0,0,.12)"; ctx.lineWidth = 1;
      if (r.whole) ctx.setLineDash([3, 3]);
      ctx.strokeRect(x + 0.5, 4.5, Math.max(w, 1), RMC - 9);
      ctx.setLineDash([]);
      if (w > 40) {
        ctx.fillStyle = a > 0.5 ? "#fff" : css("--ink2");
        ctx.font = "10px ui-monospace,monospace"; ctx.textBaseline = "middle"; ctx.textAlign = "center";
        ctx.fillText(r.oe.toFixed(2) + (r.whole ? " · gene" : ""), x + w / 2, RMC / 2);
      }
    }
  }
  function drawIsl(): void {
    const ctx = ctxFor("v-isl", ISL);
    for (const i of b.islands) {
      const x = X(i.s), w = (i.e - i.s + 1) * colW;
      ctx.fillStyle = islandColor(i); ctx.globalAlpha = 0.9;
      ctx.fillRect(x, 4, Math.max(w, 1), ISL - 8); ctx.globalAlpha = 1;
      ctx.strokeStyle = "rgba(0,0,0,.15)"; ctx.lineWidth = 1;
      ctx.strokeRect(x + 0.5, 4.5, Math.max(w, 1), ISL - 9);
    }
  }
  function roundRect(ctx: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, r: number): void {
    r = Math.min(r, w / 2, h / 2); ctx.beginPath();
    ctx.moveTo(x + r, y); ctx.arcTo(x + w, y, x + w, y + h, r); ctx.arcTo(x + w, y + h, x, y + h, r);
    ctx.arcTo(x, y + h, x, y, r); ctx.arcTo(x, y, x + w, y, r); ctx.closePath();
  }
  function drawDom(): void {
    const ctx = ctxFor("v-dom", DOM);
    for (const d of b.domains) {
      const x = X(d.s), w = (d.e - d.s + 1) * colW;
      ctx.fillStyle = css("--domain"); roundRect(ctx, x, 6, Math.max(w, 2), DOM - 12, 4); ctx.fill();
      if (w > 60) {
        ctx.fillStyle = "#fff"; ctx.font = "600 11px ui-sans-serif,system-ui";
        ctx.textBaseline = "middle"; ctx.textAlign = "center"; ctx.fillText(d.name, x + w / 2, DOM / 2);
      }
    }
  }
  function drawAxis(id: string, top: boolean): void {
    const ctx = ctxFor(id, AXIS);
    ctx.fillStyle = css("--muted"); ctx.strokeStyle = css("--border");
    ctx.textAlign = "center";
    const baseline = top ? AXIS - 1 : 1;
    ctx.beginPath(); ctx.moveTo(0, baseline + 0.5); ctx.lineTo(contentW(), baseline + 0.5); ctx.stroke();

    if (top) {
      if (b.exons && b.exons.length) {
        const ex = b.exons;
        ctx.strokeStyle = css("--muted"); ctx.globalAlpha = 0.5; ctx.lineWidth = 1;
        for (const st of ex) {
          const x = X(st);
          ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, AXIS); ctx.stroke();
        }
        ctx.globalAlpha = 1;
        ctx.fillStyle = css("--muted"); ctx.font = "9px ui-sans-serif,system-ui";
        ctx.textBaseline = "middle";
        for (let k = 0; k < ex.length; k++) {
          const s = ex[k], e = k + 1 < ex.length ? ex[k + 1] : N + 1;
          if ((e - s) * colW > 16) ctx.fillText(String(k + 1), (X(s) + X(e)) / 2, AXIS / 2);
        }
      }
      return;
    }

    ctx.font = "10px ui-monospace,monospace"; ctx.textBaseline = "bottom";
    let step = Math.max(1, Math.round(90 / colW));
    for (const s of [1, 2, 5, 10, 25, 50, 100, 200, 250, 500, 1000]) { if (s >= step) { step = s; break; } }
    for (let p = step; p <= N; p += step) {
      const x = X(p) + colW / 2;
      ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, 5); ctx.stroke();
      ctx.fillText(String(p), x, AXIS - 2);
    }
  }

  function render(): void {                   // repaint canvases for the current window
    viewW = Math.max(1, scroll.clientWidth);
    // Clamp the scroll offset to the current content: after zooming out while scrolled
    // right, scroll.scrollLeft can momentarily exceed the (now smaller) contentW, which
    // would leave the visible window past the protein and blank every lane. Correct the
    // DOM offset and the value we draw with so the view can never disappear.
    const maxScroll = Math.max(0, contentW() - viewW);
    if (scroll.scrollLeft > maxScroll) scroll.scrollLeft = maxScroll;
    scrollX = Math.min(scroll.scrollLeft, maxScroll);
    drawHeat(); drawScore(); drawClinvar(); drawIsl(); drawDom();
    if (hasHot) drawHotspots();
    if (hasRmc) drawRmc();
    drawAxis("v-axt", true); drawAxis("v-axb", false);
  }

  function drawAll(): void {                   // full layout pass (zoom / fit / resize / theme)
    stack.style.width = contentW() + "px";
    ($("v-z") as HTMLInputElement).value = String(colW);
    laneCV.h = Math.round(CVH * cvZoom());     // keep gutter + canvas in sync with zoom
    buildGutter();
    xhair.style.height = LANES.reduce((s, l) => s + l.h, 0) + "px";
    render();
  }

  const CLS = (s: number) => (s >= CUTOFF_P ? "likely pathogenic" : s <= CUTOFF_B ? "likely benign" : "ambiguous");
  const clsColor = (s: number) => (s >= CUTOFF_P ? css("--pathogenic") : s <= CUTOFF_B ? css("--benign") : css("--ambiguous"));
  function cvRow(v: CVar, p: number): string {
    const st = v.st || 0;
    const stars = "★".repeat(st) + "☆".repeat(4 - st);
    return `<div class="row" style="margin-top:2px"><span>ClinVar</span><span class="pill" style="background:${CVCOL[v.sig]}">${CVNAME[v.sig]}</span></div>`
      + `<div class="row"><span>&nbsp;</span><b>${v.ref}${p}${v.alt} · ${v.n} record${v.n > 1 ? "s" : ""}</b></div>`
      + `<div class="row"><span>review</span><span title="${st}/4 gold stars">${stars}</span></div>`;
  }
  function onMove(e: MouseEvent): void {
    if (dragX !== null) {                        // rubber-band zoom selection
      const x = plotX(e);
      zoomsel.style.left = Math.min(x, dragX) + "px";
      zoomsel.style.width = Math.abs(x - dragX) + "px";
      tip.style.opacity = "0"; xhair.style.opacity = "0";
      return;
    }
    const rect = scroll.getBoundingClientRect();
    const xIn = e.clientX - rect.left + scroll.scrollLeft;
    const pi = Math.floor(xIn / colW);
    if (pi < 0 || pi >= N) { onLeave(); return; }
    const p = pi + 1; const refAA = b.ref[pi] || "?";
    xhair.style.opacity = "0.35"; xhair.style.left = pi * colW + Math.max(colW / 2, 0.5) + "px";
    const hr = heat.getBoundingClientRect(); let rowHtml = "";
    if (e.clientY >= hr.top && e.clientY <= hr.bottom) {
      const r = Math.floor((e.clientY - hr.top) / ROWH);
      if (r >= 0 && r < AA.length) {
        const val = grid[pi * AA.length + r]; const alt = AA[r];
        let name: string, col: string;
        if (val === B_MISSING) { name = "not scored"; col = "var(--muted)"; }
        else if (val === B_REF) { name = "reference"; col = css("--reference"); }
        else { const sc = val / Q_MAX; name = `${sc.toFixed(3)} · ${CLS(sc)}`; col = clsColor(sc); }
        rowHtml = `<div class="row" style="margin-top:5px"><span>${refAA}${p}&rarr;${alt}</span><span class="pill" style="background:${col}">${name}</span></div>`;
      }
    }
    const m = b.mean[pi];
    const isl = b.islands.find((i) => p >= i.s && p <= i.e);
    const dom = b.domains.find((d) => p >= d.s && p <= d.e);
    const hot = b.hotspots?.find((h) => h.p === p);
    const rmc = b.rmc?.find((x) => p >= x.s && p <= x.e);
    const cvp = cvMap.get(p), cvb = cvbMap.get(p);
    tip.innerHTML = `<div class="th">${refAA} ${p}</div>`
      + `<div class="row"><span>Mean AM</span><b>${m != null ? m.toFixed(3) : "—"}</b></div>`
      + (cvp ? cvRow(cvp, p) : "") + (cvb ? cvRow(cvb, p) : "")
      + (hot ? `<div class="row"><span>Cancer hotspot</span><b>${hot.r}${hot.p}${hot.onc ? " · OncoKB" : ""}</b></div>`
        + `<div class="row"><span>&nbsp;cancerhotspots</span><b>${hot.n != null ? hot.n + " mutations" : ""}${hot.q != null ? ` · q=${hot.q.toExponential(1)}` : ""}</b></div>` : "")
      + (rmc ? `<div class="row"><span>Missense o/e${rmc.whole ? " (gene)" : ""}</span><b>${rmc.oe.toFixed(2)}${rmc.obs != null ? ` (${rmc.obs}/${rmc.exp} obs/exp)` : ""}</b></div>` : "")
      + (dom ? `<div class="row"><span>Domain</span><b>${dom.name}</b></div>` : "")
      + (isl ? `<div class="row"><span>AM island</span><b>${isl.s}–${isl.e} (${isl.m.toFixed(2)})</b></div>`
        + `<div class="row"><span>&nbsp;in island</span><b>${isl.plp || 0} P/LP · ${isl.blb || 0} B/LB</b></div>` : "")
      + rowHtml;
    tip.style.opacity = "1";
    let tx = e.clientX + 14, ty = e.clientY + 14;
    const tw = tip.offsetWidth, th = tip.offsetHeight;
    if (tx + tw > innerWidth - 8) tx = e.clientX - tw - 14;
    if (ty + th > innerHeight - 8) ty = e.clientY - th - 14;
    tip.style.left = tx + "px"; tip.style.top = ty + "px";
  }
  function onLeave(): void { tip.style.opacity = "0"; xhair.style.opacity = "0"; }

  function exportPNG(): void {
    const S = 2;                       // supersample for crisp text
    const GUT = 110, HEADER = 64, FOOT = 20, PAD = 16;
    const w = contentW();
    const lanesH = LANES.reduce((s, l) => s + l.h, 0);
    const W = GUT + w + PAD;
    // render every lane at full protein width so the composite below captures it all
    // (backing store is clamped in ctxFor); the viewport is repainted before returning.
    const savedVW = viewW, savedSX = scrollX;
    viewW = w; scrollX = 0; render();

    // Legend mirrors the on-screen one, including the conditional tracks; it wraps
    // to as many rows as the export width needs so nothing is clipped.
    const legendEntries: { kind: string; col?: string; label: string }[] = [
      { kind: "grad", label: "AM score 0→1" },
      { kind: "square", col: css("--reference"), label: "Reference" },
      { kind: "dot", col: css("--cv-p"), label: "ClinVar P/LP ↑" },
      { kind: "dot", col: css("--cv-blb"), label: "B/LB ↓" },
      { kind: "dot", col: css("--isl-path"), label: "island: P-rich" },
      { kind: "dot", col: css("--isl-benign"), label: "B-rich" },
      { kind: "dot", col: css("--isl-none"), label: "no ClinVar" },
      { kind: "square", col: css("--domain"), label: "Domain" },
      ...(hasHot ? [{ kind: "diamond", col: css("--hotspot"), label: "Cancer hotspot" }] : []),
      ...(hasRmc ? [{ kind: "rmc", col: css("--rmc"), label: "Missense-constrained (low o/e)" }] : []),
    ];
    const LEG_FONT = "11px ui-sans-serif,system-ui";
    const SW = 13, SWGAP = 6, ENTRY_GAP = 18, ROW_H = 18;
    function layoutLegend(c: CanvasRenderingContext2D, y0: number, maxX: number, draw: boolean): number {
      c.font = LEG_FONT; c.textBaseline = "middle"; c.textAlign = "left";
      let x = PAD, row = 0;
      for (const e of legendEntries) {
        const swW = e.kind === "grad" ? 44 : SW;
        const tw = c.measureText(e.label).width;
        const need = swW + SWGAP + tw;
        if (x + need > maxX && x > PAD) { row++; x = PAD; }
        const cy = y0 + row * ROW_H + ROW_H / 2;
        if (draw) {
          if (e.kind === "grad") {
            const g = c.createLinearGradient(x, 0, x + 44, 0);
            g.addColorStop(0, css("--am0")); g.addColorStop(0.5, css("--am2")); g.addColorStop(1, css("--am4"));
            c.fillStyle = g; c.fillRect(x, cy - 5.5, 44, 11);
          } else if (e.kind === "dot") {
            c.fillStyle = e.col!; c.beginPath(); c.arc(x + SW / 2, cy, 5, 0, 6.2832); c.fill();
          } else if (e.kind === "diamond") {
            c.fillStyle = e.col!; c.save(); c.translate(x + SW / 2, cy); c.rotate(Math.PI / 4);
            c.fillRect(-4.5, -4.5, 9, 9); c.restore();
          } else if (e.kind === "rmc") {
            const g = c.createLinearGradient(x, 0, x + SW, 0);
            g.addColorStop(0, css("--panel")); g.addColorStop(1, e.col!);
            c.fillStyle = g; c.fillRect(x, cy - 5.5, SW, 11);
          } else {
            c.fillStyle = e.col!; c.fillRect(x, cy - 5.5, SW, 11);
          }
          c.fillStyle = css("--muted"); c.fillText(e.label, x + swW + SWGAP, cy);
        }
        x += need + ENTRY_GAP;
      }
      return row + 1;
    }
    const measure = document.createElement("canvas").getContext("2d")!;
    const legendRows = layoutLegend(measure, 0, W - PAD, false);
    const LEGEND = legendRows * ROW_H + 8;
    const H = HEADER + LEGEND + lanesH + FOOT + PAD;
    const out = document.createElement("canvas");
    const S2 = Math.min(S, MAX_CANVAS_DIM / W);       // clamp export width to the browser limit
    out.width = Math.round(W * S2); out.height = Math.round(H * S2);
    const ctx = out.getContext("2d")!;
    ctx.scale(S2, S2);
    const ink = css("--ink"), ink2 = css("--ink2"), muted = css("--muted");

    ctx.fillStyle = css("--panel"); ctx.fillRect(0, 0, W, H);

    // header
    ctx.textBaseline = "alphabetic"; ctx.textAlign = "left";
    ctx.fillStyle = ink;
    ctx.font = '600 20px "Iowan Old Style",Palatino,Georgia,serif';
    ctx.fillText(`${b.gene} — AlphaMissense saturation map`, PAD, 30);
    ctx.fillStyle = ink2; ctx.font = "12px ui-monospace,monospace";
    ctx.fillText(
      `${b.uniprot} · ${b.ensembl}${b.refseq ? " · " + b.refseq : ""} · ${N} aa · `
      + `${b.islands.length} islands · ${b.clinvar.length} P/LP · ${b.clinvar_benign.length} B/LB`,
      PAD, 50);

    // legend (mirrors the on-screen one; wraps to fit the export width)
    layoutLegend(ctx, HEADER + 4, W - PAD, true);

    // tracks + gutter
    let y = HEADER + LEGEND;
    for (const L of LANES) {
      const cv = $(L.id) as HTMLCanvasElement;
      ctx.drawImage(cv, GUT, y, w, L.h);
      ctx.textAlign = "right"; ctx.textBaseline = "middle";
      if (L.id === "v-heat") {
        ctx.fillStyle = ink2; ctx.font = "10px ui-monospace,monospace";
        for (let i = 0; i < AA.length; i++) ctx.fillText(AA[i], GUT - 10, y + i * ROWH + ROWH / 2);
      } else if (L.label) {
        ctx.fillStyle = ink2; ctx.font = "11px ui-sans-serif,system-ui";
        ctx.fillText(L.id === "v-score" ? "Mean AM" : L.label, GUT - 10, y + L.h / 2);
      }
      y += L.h;
    }

    // footer credit
    ctx.textAlign = "left"; ctx.textBaseline = "alphabetic";
    ctx.fillStyle = muted; ctx.font = "10px ui-sans-serif,system-ui";
    ctx.fillText("AlphaMissense (DeepMind, CC BY 4.0) · ClinVar · UniProt domains · am-islands-viewer",
      PAD, H - 8);

    viewW = savedVW; scrollX = savedSX; render();   // restore the on-screen viewport
    out.toBlob((blob) => {
      if (!blob) return;
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url; a.download = `${b.gene}_alphamissense.png`;
      a.click(); URL.revokeObjectURL(url);
    }, "image/png");
  }

  function startDrag(e: MouseEvent): void {
    if (e.button !== 0) return;
    dragX = plotX(e);
    zoomsel.style.left = dragX + "px"; zoomsel.style.width = "0px";
    zoomsel.style.height = LANES.reduce((s, l) => s + l.h, 0) + "px";
    zoomsel.style.opacity = "1";
    e.preventDefault();
  }
  function endDrag(e: MouseEvent): void {
    if (dragX === null) return;
    const lo = Math.min(plotX(e), dragX), hi = Math.max(plotX(e), dragX);
    dragX = null; zoomsel.style.opacity = "0"; zoomsel.style.width = "0px";
    if (hi - lo < 5) return;                     // treat as a click, not a zoom
    const resLo = Math.max(1, Math.floor(lo / colW) + 1);
    const resHi = Math.min(N, Math.ceil(hi / colW));
    const span = Math.max(1, resHi - resLo + 1);
    colW = Math.min(MAXW, Math.max(0.35, (scroll.clientWidth - 2) / span));
    drawAll();
    scroll.scrollLeft = (resLo - 1) * colW;
  }
  let scrollRAF = 0;
  scroll.addEventListener("scroll", () => {
    if (scrollRAF) return;
    scrollRAF = requestAnimationFrame(() => { scrollRAF = 0; render(); });
  });
  scroll.addEventListener("mousedown", startDrag);
  scroll.addEventListener("mousemove", onMove);
  scroll.addEventListener("mouseup", endDrag);
  scroll.addEventListener("mouseleave", (e) => { if (dragX !== null) endDrag(e); else onLeave(); });
  function fit(): void { colW = Math.max(0.35, (scroll.clientWidth - 2) / N); drawAll(); scroll.scrollLeft = 0; }
  $("v-fit").addEventListener("click", fit);
  $("v-png").addEventListener("click", exportPNG);
  ($("v-z") as HTMLInputElement).addEventListener("input", (e) => {
    colW = parseFloat((e.target as HTMLInputElement).value); drawAll();
  });
  let rt: number | undefined;
  const onResize = () => { window.clearTimeout(rt); rt = window.setTimeout(drawAll, 150); };
  window.addEventListener("resize", onResize);
  const mo = new MutationObserver(() => { refreshColors(); drawAll(); });
  mo.observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });

  buildGutter();
  fit();
}
