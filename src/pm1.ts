/** PM1 prediction for a single variant (GENE:p.… or GENE:c.…), fully client-side. */
import type { Bundle, Island } from "./viewer";

const AA3: Record<string, string> = {
  ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C", GLN: "Q", GLU: "E", GLY: "G",
  HIS: "H", ILE: "I", LEU: "L", LYS: "K", MET: "M", PHE: "F", PRO: "P", SER: "S",
  THR: "T", TRP: "W", TYR: "Y", VAL: "V", TER: "*",
};

export interface ParsedVariant {
  gene: string;
  kind: "p" | "c";
  aaPos: number | null;   // null when non-coding (intron/UTR)
  ref?: string;           // 1-letter ref aa (p.)
  alt?: string;           // 1-letter alt aa (p.)
  raw: string;
  note?: string;          // e.g. non-coding / nonsense
}

function aa1(tok: string | undefined): string | undefined {
  if (!tok) return undefined;
  if (tok.length === 1) return tok.toUpperCase();
  return AA3[tok.toUpperCase()];
}

/** Parse "GENE:p.Arg350Gly" / "GENE p.R350G" / "GENE:c.1048C>T" (case-insensitive
 *  p./c.). Returns null if it does not look like a parseable variant. */
export function parseVariant(input: string): ParsedVariant | null {
  const s = input.trim();
  const m = s.match(/^([A-Za-z0-9_-]+)\s*[:\s]\s*([pcPC])\.\s*(.+)$/);
  if (!m) return null;
  const gene = m[1].toUpperCase();
  const kind = m[2].toLowerCase() as "p" | "c";
  const body = m[3].replace(/\s+/g, "");

  if (kind === "p") {
    const pm = body.replace(/[()]/g, "").match(/^([A-Za-z]{3}|[A-Za-z])(\d+)([A-Za-z]{3}|[A-Za-z*=]{0,3})/);
    if (!pm) return null;
    const ref = aa1(pm[1]);
    const alt = pm[3] === "=" ? ref : aa1(pm[3]);
    const aaPos = parseInt(pm[2], 10);
    let note: string | undefined;
    if (alt === "*") note = "nonsense (stop-gain) — PM1 targets missense hotspots";
    else if (pm[3] === "=") note = "synonymous — PM1 targets missense hotspots";
    return { gene, kind, aaPos, ref, alt, raw: s, note };
  }
  // c.
  const cm = body.match(/^(-?\*?)(\d+)([+-]\d+)?/);
  if (!cm) return null;
  if (cm[1] === "-" || cm[1] === "*")
    return { gene, kind, aaPos: null, raw: s, note: "UTR (non-coding) — PM1 not applicable" };
  if (cm[3])
    return { gene, kind, aaPos: null, raw: s, note: "intronic (non-coding) — PM1 not applicable" };
  const cpos = parseInt(cm[2], 10);
  return { gene, kind, aaPos: Math.ceil(cpos / 3), raw: s, note: `codon from c.${cpos} (≈ ⌈c/3⌉)` };
}

/** Heuristic: does the query look like an attempt at a variant (so a parse failure
 *  should show a variant error, not a gene-not-found)? */
export function looksLikeVariant(q: string): boolean {
  return /[:\s][pcPC]\./.test(q) || /\b[pcPC]\.\s*[A-Za-z0-9]/.test(q) || /^[A-Za-z0-9_-]+\s*:/.test(q);
}

const css = (v: string) =>
  getComputedStyle(document.documentElement).getPropertyValue(v).trim();

type Grade = { label: string; points: number; cls: "none" | "sup" | "mod" | "str" };

export function evaluatePM1(b: Bundle, aaPos: number): {
  grade: Grade; island: Island | null; meanAt: number | null; domain: string | null;
} {
  const island = b.islands.find((i) => aaPos >= i.s && aaPos <= i.e) || null;
  const meanAt = aaPos >= 1 && aaPos <= b.length ? b.mean[aaPos - 1] : null;
  const domain = (b.domains.find((d) => aaPos >= d.s && aaPos <= d.e) || null)?.name ?? null;

  let grade: Grade;
  if (!island) {
    grade = { label: "Not met", points: 0, cls: "none" };
  } else {
    const plp = island.plp || 0, blb = island.blb || 0;
    if (blb > 0 && blb >= plp) grade = { label: "Not met — contradicted by benign variation", points: 0, cls: "none" };
    else if (plp >= 10 && blb === 0) grade = { label: "Strong", points: 4, cls: "str" };
    else if (plp >= 3 && blb === 0) grade = { label: "Moderate", points: 2, cls: "mod" };
    else grade = { label: "Supporting", points: 1, cls: "sup" };
  }
  return { grade, island, meanAt, domain };
}

function nearestIsland(b: Bundle, pos: number): Island | null {
  let best: Island | null = null, bd = Infinity;
  for (const i of b.islands) {
    const d = pos < i.s ? i.s - pos : pos > i.e ? pos - i.e : 0;
    if (d < bd) { bd = d; best = i; }
  }
  return best;
}

function decodeCell(b: Bundle, pos: number, alt: string): number | null {
  if (!b.grid) return null;
  const naa = b.aa_order.length;
  const idx = b.aa_order.indexOf(alt);
  if (idx < 0) return null;
  const bin = atob(b.grid);
  const v = bin.charCodeAt((pos - 1) * naa + idx);
  if (v === 251 || v === 252 || Number.isNaN(v)) return null;
  return v / 250;
}

export function renderPM1(root: HTMLElement, pv: ParsedVariant, b: Bundle): void {
  const aaPos = pv.aaPos;
  const outOfRange = aaPos != null && (aaPos < 1 || aaPos > b.length);
  const refAt = aaPos && !outOfRange ? b.ref[aaPos - 1] : "";
  const refMismatch = pv.ref && refAt && pv.ref !== refAt;

  const ev = aaPos && !outOfRange ? evaluatePM1(b, aaPos) : null;
  const g = ev?.grade;
  const isl = ev?.island || null;
  const cellScore = aaPos && !outOfRange && pv.alt && pv.alt !== "*" ? decodeCell(b, aaPos, pv.alt) : null;

  const posStr = aaPos != null ? `${refAt || pv.ref || "?"}${aaPos}${pv.alt && pv.alt !== "*" ? pv.alt : ""}` : "—";

  // justification
  let just = "";
  if (aaPos == null) {
    just = `This variant is ${pv.note}. PM1 assesses missense mutational hotspots, so it does not apply here.`;
  } else if (outOfRange) {
    just = `Residue ${aaPos} is outside ${pv.gene} (length ${b.length} aa) on the AlphaMissense isoform — check the transcript/notation.`;
  } else if (isl && g) {
    const dom = ev!.domain ? ` within the ${ev!.domain} domain` : "";
    const cv = `${isl.plp || 0} Pathogenic/Likely-path. and ${isl.blb || 0} Benign/Likely-benign ClinVar missense residue(s)`;
    if (g.cls === "none") {
      just = `Residue ${aaPos} lies in an AlphaMissense island (residues ${isl.s}–${isl.e}, mean AM ${isl.m.toFixed(2)})${dom}, but that island carries ${cv} — the benign signal outweighs the pathogenic one, so PM1 is not applied.`;
    } else {
      const strengthWhy = g.cls === "str" ? "a densely pathogenic, benign-free hot-spot"
        : g.cls === "mod" ? "pathogenic-enriched with no benign variation"
        : "supported by the predicted hot-spot with limited ClinVar corroboration";
      just = `Residue ${aaPos} lies within an AlphaMissense island (residues ${isl.s}–${isl.e}, mean AM ${isl.m.toFixed(2)})${dom}. The island contains ${cv} — ${strengthWhy}. PM1 is met at ${g.label} strength (+${g.points}).`;
    }
  } else {
    const ni = nearestIsland(b, aaPos);
    const near = ni ? ` The nearest island is residues ${ni.s}–${ni.e} (${ni.s > aaPos ? ni.s - aaPos : aaPos - ni.e} aa away).` : "";
    just = `Residue ${aaPos} is not within any AlphaMissense island, so PM1 is not met.${near}`;
  }

  const gradeCls = g?.cls ?? "none";
  const gradeLabel = aaPos == null || outOfRange ? "N/A" : (g?.label ?? "Not met");
  const pts = g && g.points > 0 ? `+${g.points}` : "";

  root.innerHTML = `
    <div class="pm1">
      <p class="eyebrow">PM1 prediction · AlphaMissense island</p>
      <h1><em>${pv.gene}</em> <span class="mono pm1-var">${pv.raw.replace(/^[^:\s]+[:\s]\s*/, "")}</span></h1>
      <div class="pm1-grid">
        <div class="pm1-verdict grade-${gradeCls}">
          <div class="pm1-verdict-l">PM1</div>
          <div class="pm1-verdict-n">${gradeLabel}</div>
          ${pts ? `<div class="pm1-verdict-p">${pts} point${g!.points === 1 ? "" : "s"}</div>` : ""}
        </div>
        <div class="pm1-facts">
          <div class="fact"><span>Residue</span><b class="mono">${posStr}</b></div>
          <div class="fact"><span>Mean AM at residue</span><b>${ev?.meanAt != null ? ev.meanAt.toFixed(3) : "—"}</b></div>
          ${cellScore != null ? `<div class="fact"><span>This substitution (AM)</span><b>${cellScore.toFixed(3)} · ${cellScore >= 0.564 ? "likely pathogenic" : cellScore <= 0.34 ? "likely benign" : "ambiguous"}</b></div>` : ""}
          <div class="fact"><span>In AM island</span><b>${isl ? `${isl.s}–${isl.e} (mean ${isl.m.toFixed(2)})` : "no"}</b></div>
          ${isl ? `<div class="fact"><span>Island ClinVar</span><b>${isl.plp || 0} P/LP · ${isl.blb || 0} B/LB</b></div>` : ""}
          <div class="fact"><span>Domain</span><b>${ev?.domain || "—"}</b></div>
        </div>
      </div>
      <p class="pm1-just">${just}</p>
      ${refMismatch ? `<p class="pm1-warn">⚠︎ The reference amino acid you gave (${pv.ref}) does not match residue ${aaPos} on the AlphaMissense isoform (${refAt}). The notation may be on a different transcript; the position was used as given.</p>` : ""}
      ${pv.kind === "c" && aaPos != null ? `<p class="pm1-warn subtle">Note: the amino-acid position was derived from the coding position as ⌈c/3⌉, assuming the coding coordinate matches the AlphaMissense isoform.</p>` : ""}
      ${pv.note && aaPos != null && pv.kind === "p" && (pv.note.startsWith("nonsense") || pv.note.startsWith("synonymous")) ? `<p class="pm1-warn subtle">Note: ${pv.note}.</p>` : ""}
      <figure class="pm1-fig"><canvas id="pm1-canvas"></canvas>
        <figcaption>Mean AlphaMissense (area) with islands (colored by ClinVar overlap) and ClinVar
          lollipops around the variant (▼). Dashed line = 0.564 island threshold.</figcaption>
      </figure>
      <p class="pm1-actions">
        <a class="ghost" href="?gene=${pv.gene}" data-gene-link="${pv.gene}">Open full ${pv.gene} map →</a>
        <span class="pm1-strengths">Strength scale: <span class="pill sup">Supporting +1</span> <span class="pill mod">Moderate +2</span> <span class="pill str">Strong +4</span></span>
      </p>
      <p class="pm1-caveat">PM1 strength here follows a proposed framework combining island membership with in-island ClinVar enrichment; final cut-offs are gene/VCEP-specific.</p>
    </div>`;

  if (aaPos != null && !outOfRange) drawFig(root, b, aaPos);
}

function drawFig(root: HTMLElement, b: Bundle, aaPos: number): void {
  const cv = root.querySelector<HTMLCanvasElement>("#pm1-canvas");
  if (!cv) return;
  const isl = b.islands.find((i) => aaPos >= i.s && aaPos <= i.e);
  let lo = Math.max(1, aaPos - 110), hi = Math.min(b.length, aaPos + 110);
  if (isl) { lo = Math.max(1, Math.min(lo, isl.s - 20)); hi = Math.min(b.length, Math.max(hi, isl.e + 20)); }
  const span = hi - lo + 1;

  const W = Math.min(880, root.clientWidth - 8) || 800;
  const H = 170, PADB = 26, PADT = 10, dpr = Math.max(1, window.devicePixelRatio || 1);
  cv.style.width = W + "px"; cv.style.height = H + "px";
  cv.width = W * dpr; cv.height = H * dpr;
  const ctx = cv.getContext("2d")!; ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const X = (p: number) => ((p - lo) / span) * W;
  const Y = (v: number) => PADT + (1 - v) * (H - PADB - PADT);

  // islands
  for (const i of b.islands) {
    if (i.e < lo || i.s > hi) continue;
    const plp = i.plp || 0, blb = i.blb || 0;
    const frac = plp + blb === 0 ? null : (plp - blb) / (plp + blb);
    ctx.fillStyle = frac == null ? css("--isl-none") : frac >= 0 ? css("--isl-path") : css("--isl-benign");
    ctx.globalAlpha = 0.16;
    ctx.fillRect(X(Math.max(i.s, lo)), PADT, X(Math.min(i.e, hi)) - X(Math.max(i.s, lo)), H - PADB - PADT);
    ctx.globalAlpha = 1;
  }
  // threshold
  ctx.strokeStyle = css("--muted"); ctx.setLineDash([4, 4]); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(0, Y(0.564) + 0.5); ctx.lineTo(W, Y(0.564) + 0.5); ctx.stroke();
  ctx.setLineDash([]);
  // mean-AM area
  const area = new Path2D(); let started = false;
  for (let p = lo; p <= hi; p++) {
    const m = b.mean[p - 1]; if (m == null) continue;
    const x = X(p), y = Y(m);
    if (!started) { area.moveTo(x, Y(0)); area.lineTo(x, y); started = true; } else area.lineTo(x, y);
  }
  area.lineTo(X(hi), Y(0)); area.closePath();
  ctx.fillStyle = css("--score-fill"); ctx.fill(area);
  ctx.beginPath(); started = false;
  for (let p = lo; p <= hi; p++) {
    const m = b.mean[p - 1]; if (m == null) continue;
    const x = X(p), y = Y(m); if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
  }
  ctx.strokeStyle = css("--score-stroke"); ctx.lineWidth = 1.4; ctx.stroke();
  // clinvar lollipops (small)
  const base = H - PADB;
  const dot = (p: number, up: boolean, color: string) => {
    if (p < lo || p > hi) return;
    const x = X(p), y = up ? base - 10 : base + 6;
    ctx.strokeStyle = css("--muted"); ctx.globalAlpha = 0.4;
    ctx.beginPath(); ctx.moveTo(x, base); ctx.lineTo(x, y); ctx.stroke(); ctx.globalAlpha = 1;
    ctx.fillStyle = color; ctx.beginPath(); ctx.arc(x, y, 2.4, 0, 6.2832); ctx.fill();
  };
  b.clinvar.forEach((v) => dot(v.p, true, css("--cv-plp")));
  b.clinvar_benign.forEach((v) => dot(v.p, false, css("--cv-blb")));
  // variant marker
  const vx = X(aaPos);
  ctx.strokeStyle = css("--ink"); ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(vx, PADT); ctx.lineTo(vx, base); ctx.stroke();
  ctx.fillStyle = css("--ink");
  ctx.beginPath(); ctx.moveTo(vx - 5, PADT - 1); ctx.lineTo(vx + 5, PADT - 1); ctx.lineTo(vx, PADT + 6); ctx.closePath(); ctx.fill();
  // axis ticks
  ctx.fillStyle = css("--muted"); ctx.font = "10px ui-monospace,monospace"; ctx.textBaseline = "top"; ctx.textAlign = "center";
  const step = span > 400 ? 100 : span > 200 ? 50 : 20;
  for (let p = Math.ceil(lo / step) * step; p <= hi; p += step) ctx.fillText(String(p), X(p), base + 6);
  ctx.textAlign = "center"; ctx.fillStyle = css("--ink"); ctx.font = "600 10px ui-sans-serif,system-ui";
  ctx.fillText(String(aaPos), vx, base + 6);
}
