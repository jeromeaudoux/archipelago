/** PM1 benchmark explorer: the SeqOne engine's island-PM1 vs ClinGen-curated PM1
 *  (eRepo), matched by ClinVar ID — the same computation as the published benchmark. */
import { parseVariant, renderPM1 } from "./pm1";
import type { Bundle } from "./viewer";

interface Row {
  g: string; v: string; c: string; hgvs: string; dis: string; panel: string; link: string;
  clinvar: string; assertion: string; met: string;
  truth: boolean; tstr: string; pred: boolean; pstr: string; pr: string; cat: "tp" | "fp" | "fn";
  dcat?: string; pnear?: number; bnear?: number; idist?: number; dom?: boolean; amv?: number | null;
}
interface Data {
  summary: { tp: number; fp: number; fn: number; tn: number;
    precision: number; recall: number; f1: number; matched: number;
    fn_not_missense: number; fn_no_island: number; dcats?: Record<string, number> };
  variants: Row[];
}
type LoadBundle = (gene: string) => Promise<Bundle>;

const CAT = {
  tp: { label: "True positive", cls: "tp" },
  fp: { label: "False positive", cls: "fp" },
  fn: { label: "False negative", cls: "fn" },
};

/** Honest, data-derived reasons an FP/FN diverges from ClinGen. */
const DCAT: Record<string, { label: string; blurb: string; grp: "fp" | "fn" | "na" }> = {
  fp_clinvar_hotspot: { grp: "fp", label: "ClinVar-supported hotspot",
    blurb: "Region is dense with ClinVar pathogenic missense (≥3 within ±10 aa); ClinGen simply hasn’t curated PM1 here — a defensible call rather than a clear error." },
  fp_am_only: { grp: "fp", label: "AlphaMissense-only island",
    blurb: "AlphaMissense flags the region but little/no ClinVar pathogenic support — a genuine AM-versus-panel divergence." },
  fp_benign_conflict: { grp: "fp", label: "Benign variants nearby",
    blurb: "≥2 benign/LB variants sit in the region, weakening the hot-spot claim — the engine likely over-called." },
  fn_subthreshold_hotspot: { grp: "fn", label: "Sub-threshold hotspot",
    blurb: "A real ClinVar pathogenic hotspot (≥3 P/LP within ±10 aa) where the AlphaMissense signal wasn’t strong enough to form an island — the dominant source of misses." },
  fn_near_island: { grp: "fn", label: "Just outside an island",
    blurb: "Residue lies ≤10 aa beyond an island boundary — a boundary-sensitivity miss." },
  fn_weak_am: { grp: "fn", label: "Sparse evidence",
    blurb: "Little ClinVar and sub-island AM here; the panel likely used non-AM evidence (functional assays, specific-codon rules) the island method can’t see." },
  fn_not_missense: { grp: "fn", label: "Not missense (by design)",
    blurb: "Frameshift / nonsense / splice / non-single-residue — PM1 is a missense criterion, so PM1_AM is correctly not applied." },
  unmapped: { grp: "na", label: "Isoform-ambiguous",
    blurb: "Maps to a different isoform than the AlphaMissense canonical sequence, so its neighbourhood can’t be characterised reliably." },
};
const DCAT_ORDER = ["fp_clinvar_hotspot", "fp_am_only", "fp_benign_conflict",
  "fn_subthreshold_hotspot", "fn_near_island", "fn_weak_am", "fn_not_missense"];

export function renderBenchmark(root: HTMLElement, data: Data, loadBundle: LoadBundle): void {
  const s = data.summary;
  let mode: "variant" | "gene" = "variant";
  let filter: "all" | "tp" | "fp" | "fn" = "all";
  let dfilter: string | null = null;
  let query = "";
  let gsort = { key: "n", dir: -1 };

  const dc = s.dcats || {};
  function dcatColumn(title: string, grp: "fp" | "fn"): string {
    const cats = DCAT_ORDER.filter((k) => DCAT[k].grp === grp && dc[k]);
    const tot = cats.reduce((n, k) => n + dc[k], 0) || 1;
    return `<div class="dcol"><h3 class="dcol-h ${grp}">${title}</h3>${cats.map((k) => {
      const n = dc[k], p = Math.round((n / tot) * 100);
      return `<button class="dcat" data-dcat="${k}">
        <span class="dcat-top"><span class="dcat-n ${grp}">${n}</span><span class="dcat-l">${DCAT[k].label}</span><span class="dcat-pct">${p}%</span></span>
        <span class="dcat-track"><span class="dcat-fill ${grp}" style="width:${p}%"></span></span>
        <span class="dcat-blurb">${DCAT[k].blurb}</span></button>`;
    }).join("")}</div>`;
  }
  const dcatBreakdown = () => `
    <details class="dbreak" open>
      <summary class="dbreak-sum">Where do the engine and ClinGen diverge?
        <span class="muted">— an honest categorisation of every FP and FN; click one to filter the table</span></summary>
      <div class="dcols">
        ${dcatColumn("False positives — engine applied PM1, ClinGen didn’t", "fp")}
        ${dcatColumn("False negatives — ClinGen applied PM1, engine didn’t", "fn")}
      </div>
      ${dc.unmapped ? `<p class="dcat-note"><b>${dc.unmapped}</b> further variant${dc.unmapped === 1 ? " is" : "s are"} isoform-ambiguous
        (mapped to a different isoform than the AlphaMissense canonical sequence) and left uncategorised.</p>` : ""}
    </details>`;

  // per-gene aggregates
  const geneAgg = new Map<string, { tp: number; fp: number; fn: number }>();
  for (const r of data.variants) {
    const a = geneAgg.get(r.g) || { tp: 0, fp: 0, fn: 0 };
    a[r.cat]++; geneAgg.set(r.g, a);
  }
  const genes = [...geneAgg.entries()].map(([g, a]) => ({
    g, ...a, prec: a.tp + a.fp ? a.tp / (a.tp + a.fp) : null,
    rec: a.tp + a.fn ? a.tp / (a.tp + a.fn) : null, n: a.tp + a.fp + a.fn,
  }));

  root.innerHTML = `
    <div class="bench-page">
      <p class="eyebrow">PM1 benchmark · ClinGen Evidence Repository</p>
      <h1>Explore <em>ClinGen-curated PM1</em></h1>
      <p class="bench-lead">The SeqOne engine's island-based <b>PM1</b> versus ClinGen expert-panel
        curation, matched by ClinVar ID over ${s.matched.toLocaleString()} eRepo variants — the same
        computation as the published benchmark.</p>

      <div class="bench-tiles">
        <button class="btile" data-f="tp"><span class="bt-n">${s.tp}</span><span class="bt-l">True positive</span></button>
        <button class="btile" data-f="fp"><span class="bt-n">${s.fp}</span><span class="bt-l">False positive</span></button>
        <button class="btile" data-f="fn"><span class="bt-n">${s.fn}</span><span class="bt-l">False negative</span></button>
        <div class="btile stat"><span class="bt-n">${(s.precision * 100).toFixed(0)}%</span><span class="bt-l">Precision</span></div>
        <div class="btile stat"><span class="bt-n">${(s.recall * 100).toFixed(0)}%</span><span class="bt-l">Recall</span></div>
        <div class="btile stat"><span class="bt-n">${s.f1.toFixed(2)}</span><span class="bt-l">F1</span></div>
      </div>

      ${dcatBreakdown()}

      <div class="bench-controls">
        <div class="bmode">
          <button class="bmodebtn on" data-m="variant">By variant</button>
          <button class="bmodebtn" data-m="gene">By gene (${genes.length})</button>
        </div>
        <div class="bfilters" id="bfilters">
          <button class="bfilter on" data-f="all">All (${data.variants.length})</button>
          <button class="bfilter" data-f="tp">TP ${s.tp}</button>
          <button class="bfilter" data-f="fp">FP ${s.fp}</button>
          <button class="bfilter" data-f="fn">FN ${s.fn}</button>
        </div>
        <input id="bench-search" class="bsearch" type="search" placeholder="Filter by gene, variant or disease…" spellcheck="false" />
      </div>

      <div class="btable-wrap"><table class="btable"><thead id="bench-head"></thead><tbody id="bench-body"></tbody></table></div>
      <p id="bench-count" class="bench-count"></p>
      <p class="bench-caveat">Prediction = the engine's <span class="mono">acmg_tags_v2</span> PM1 state
        (AlphaMissense-island overlap). Ground truth = PM1 (any strength) in the eRepo “Applied Evidence
        Codes (Met)”, matched by ClinVar ID and compared strength-independently, exactly as in
        <span class="mono">benchmark_acmg.py</span>. Divergence categories above are derived from each
        variant’s AlphaMissense neighbourhood (island overlap, ClinVar P/LP and B/LB density within
        ±10 aa) and are descriptive, not a re-scoring.</p>
      <div id="bench-modal"></div>
    </div>`;

  const head = root.querySelector<HTMLElement>("#bench-head")!;
  const body = root.querySelector<HTMLElement>("#bench-body")!;
  const countEl = root.querySelector<HTMLElement>("#bench-count")!;
  const modalHost = root.querySelector<HTMLElement>("#bench-modal")!;
  const search = root.querySelector<HTMLInputElement>("#bench-search")!;
  let view: Row[] = [];

  const badge = (str: string, on: boolean) =>
    `<span class="pm1badge ${on ? "on" : "off"}">${str || (on ? "PM1" : "not applied")}</span>`;
  const engineBadge = (r: Row) =>
    r.pred ? `<span class="pm1badge on">${r.pstr || "PM1"}</span>`
    : r.pr === "notmiss" ? `<span class="pm1badge na" title="PM1 applies only to missense / in-frame variants">N/A · not missense</span>`
    : `<span class="pm1badge off">not met</span>`;
  const geneLink = (g: string) => `<a href="?gene=${g}" data-gene-link="${g}" class="blink bg">${g}</a>`;

  const dcatCell = (r: Row): string => {
    if (r.cat === "tp" || !r.dcat) return `<span class="muted">—</span>`;
    const meta = DCAT[r.dcat];
    if (!meta) return `<span class="muted">—</span>`;
    const near = r.pnear != null
      ? ` — ${r.pnear} P/LP · ${r.bnear} B/LB within ±10 aa${r.amv != null ? ` · AM ${r.amv}` : ""}` : "";
    return `<span class="dpill ${meta.grp}" title="${meta.blurb}${near}">${meta.label}</span>`;
  };

  function renderVariantTable(): void {
    head.innerHTML = `<tr><th>Gene</th><th>Variant</th><th>Disease</th>
      <th>ClinGen PM1</th><th>Engine PM1</th><th>Result</th><th>Why divergent</th></tr>`;
    const q = query.toLowerCase();
    view = data.variants.filter((r) =>
      (filter === "all" || r.cat === filter)
      && (!dfilter || r.dcat === dfilter)
      && (!q || r.g.toLowerCase().includes(q) || (r.hgvs + r.v).toLowerCase().includes(q) || r.dis.toLowerCase().includes(q)));
    body.innerHTML = view.slice(0, 3000).map((r, i) => `<tr data-i="${i}">
      <td>${geneLink(r.g)}</td>
      <td class="mono bv">${r.hgvs || r.v}</td>
      <td class="bdis" title="${r.dis}">${r.dis || "—"}</td>
      <td>${badge(r.tstr, r.truth)}</td>
      <td>${engineBadge(r)}</td>
      <td><span class="rescat ${CAT[r.cat].cls}">${CAT[r.cat].label}</span></td>
      <td>${dcatCell(r)}</td>
    </tr>`).join("");
    const tag = dfilter && DCAT[dfilter] ? ` · ${DCAT[dfilter].label}` : "";
    countEl.textContent = `${view.length.toLocaleString()} variant${view.length === 1 ? "" : "s"}${tag}`
      + (view.length > 3000 ? " (showing first 3,000)" : "");
  }

  function bar(pct: number | null, kind: string): string {
    if (pct == null) return `<span class="gnum">—</span>`;
    const p = Math.round(pct * 100);
    return `<span class="gauge"><span class="gauge-bar ${kind}"><span style="width:${p}%"></span></span><span class="gnum">${p}%</span></span>`;
  }
  const SORTS: Record<string, (x: typeof genes[number]) => number | string> = {
    g: (x) => x.g, tp: (x) => x.tp, fp: (x) => x.fp, fn: (x) => x.fn,
    prec: (x) => (x.prec == null ? -1 : x.prec), rec: (x) => (x.rec == null ? -1 : x.rec), n: (x) => x.n,
  };
  function renderGeneTable(): void {
    const ind = (k: string) => gsort.key === k ? (gsort.dir < 0 ? " ▼" : " ▲") : "";
    head.innerHTML = `<tr>
      <th class="sortable" data-sort="g">Gene${ind("g")}</th>
      <th class="num sortable" data-sort="tp">TP${ind("tp")}</th>
      <th class="num sortable" data-sort="fp">FP${ind("fp")}</th>
      <th class="num sortable" data-sort="fn">FN${ind("fn")}</th>
      <th class="gcol sortable" data-sort="prec">Precision${ind("prec")}</th>
      <th class="gcol sortable" data-sort="rec">Sensitivity${ind("rec")}</th></tr>`;
    const q = query.toLowerCase();
    const keyfn = SORTS[gsort.key];
    const gs = genes.filter((x) => !q || x.g.toLowerCase().includes(q)).sort((a, b) => {
      const ka = keyfn(a), kb = keyfn(b);
      const c = ka < kb ? -1 : ka > kb ? 1 : 0;
      return (c || a.g.localeCompare(b.g)) * gsort.dir;
    });
    body.innerHTML = gs.map((x) => `<tr data-gene="${x.g}">
      <td>${geneLink(x.g)}</td>
      <td class="num">${x.tp}</td><td class="num">${x.fp}</td><td class="num">${x.fn}</td>
      <td class="gcol">${bar(x.prec, "prec")}</td>
      <td class="gcol">${bar(x.rec, "rec")}</td>
    </tr>`).join("");
    countEl.textContent = `${gs.length.toLocaleString()} gene${gs.length === 1 ? "" : "s"} · click a row for its variants`;
  }

  function render(): void { mode === "variant" ? renderVariantTable() : renderGeneTable(); }

  async function openModal(r: Row): Promise<void> {
    const meta = r.dcat ? DCAT[r.dcat] : undefined;
    const dline = meta && r.cat !== "tp" ? `<div class="bmodal-dcat ${meta.grp}">
      <span class="dpill ${meta.grp}">${meta.label}</span>
      <span class="bmodal-dblurb">${meta.blurb}${r.pnear != null
        ? ` <span class="muted">(${r.pnear} P/LP · ${r.bnear} B/LB within ±10 aa${r.amv != null ? ` · AM mean ${r.amv}` : ""}${r.dom ? " · in domain" : ""})</span>` : ""}</span>
    </div>` : "";
    const strip = `<div class="bmodal-strip">
      <span class="rescat ${CAT[r.cat].cls}">${CAT[r.cat].label}</span>
      <span class="bstrip-item">ClinGen: <b>${r.truth ? (r.tstr || "PM1") : "PM1 not applied"}</b> · ${r.assertion || "—"}</span>
      <span class="bstrip-item">${r.panel || ""}</span>
      <span class="bstrip-codes">${(r.met || "").split(",").map((c) => c.trim()).filter(Boolean)
        .map((c) => `<span class="code ${c.startsWith("PM1") ? "pm1" : ""}">${c}</span>`).join("")}</span>
      ${r.link ? `<a class="blink" href="${r.link}" target="_blank" rel="noopener">eRepo ↗</a>` : ""}
    </div>${dline}`;
    modalHost.innerHTML = `<div class="modal-backdrop" data-close>
      <div class="modal bmodal-lg" role="dialog" aria-modal="true">
        <button class="modal-x bmodal-close" data-close aria-label="Close">×</button>
        ${strip}
        <div id="bmodal-pm1"></div>
      </div>
    </div>`;
    const host = modalHost.querySelector<HTMLElement>("#bmodal-pm1")!;
    const vstr = r.v ? `p.${r.v}` : (r.c || "");
    if (vstr) {
      const pv = parseVariant(`${r.g}:${vstr}`);
      try {
        const bundle = await loadBundle(r.g);
        if (pv) renderPM1(host, pv, bundle);
        else host.innerHTML = `<p class="pm1-just">Couldn't parse ${r.g} ${r.hgvs}.</p>`;
      } catch {
        host.innerHTML = `<p class="pm1-just">Couldn't load ${r.g}.</p>`;
      }
    } else {
      host.innerHTML = `<div class="pm1"><h1><em>${r.g}</em> <span class="mono pm1-var">${r.hgvs}</span></h1>
        <p class="pm1-just">This eRepo variant is not a single-residue missense change, so the residue-level
        island view isn't shown. ClinGen applied ${r.truth ? (r.tstr || "PM1") : "no PM1"};
        the engine ${r.pred ? "met" : "did not meet"} PM1.</p>
        <p class="pm1-actions"><a class="ghost" href="?gene=${r.g}" data-gene-link="${r.g}">Open ${r.g} map →</a></p></div>`;
    }
  }
  const closeModal = () => { modalHost.innerHTML = ""; };

  head.addEventListener("click", (e) => {
    if (mode !== "gene") return;
    const th = (e.target as HTMLElement).closest<HTMLElement>("th[data-sort]");
    if (!th) return;
    const k = th.dataset.sort!;
    if (gsort.key === k) gsort.dir *= -1;
    else gsort = { key: k, dir: k === "g" ? 1 : -1 };
    renderGeneTable();
  });

  body.addEventListener("click", (e) => {
    const t = e.target as HTMLElement;
    if (t.closest("a")) return;                       // gene link → let main.ts navigate
    if (mode === "gene") {
      const gtr = t.closest<HTMLElement>("tr[data-gene]");
      if (gtr) { mode = "variant"; dfilter = null; query = gtr.dataset.gene!; search.value = query;
        root.querySelectorAll(".bmodebtn").forEach((b) => b.classList.toggle("on", (b as HTMLElement).dataset.m === "variant"));
        (root.querySelector("#bfilters") as HTMLElement).style.visibility = "visible";
        render(); }
      return;
    }
    const tr = t.closest<HTMLElement>("tr[data-i]");
    if (tr) openModal(view[parseInt(tr.dataset.i!, 10)]);
  });
  modalHost.addEventListener("click", (e) => {
    const t = e.target as HTMLElement;
    if (t.closest("[data-close]") && !t.closest("#bmodal-pm1")) closeModal();
    else if (t.closest("[data-gene-link]")) closeModal();
  });
  root.addEventListener("keydown", (e) => { if ((e as KeyboardEvent).key === "Escape") closeModal(); });

  const syncChips = () => {
    root.querySelectorAll(".bmodebtn").forEach((b) => b.classList.toggle("on", (b as HTMLElement).dataset.m === mode));
    root.querySelectorAll(".bfilter").forEach((f) => f.classList.toggle("on", (f as HTMLElement).dataset.f === filter));
    root.querySelectorAll(".dcat").forEach((d) => d.classList.toggle("on", (d as HTMLElement).dataset.dcat === dfilter));
  };
  root.querySelectorAll<HTMLElement>(".bmodebtn").forEach((el) =>
    el.addEventListener("click", () => {
      mode = (el.dataset.m as typeof mode) || "variant"; dfilter = null;
      syncChips();
      (root.querySelector("#bfilters") as HTMLElement).style.visibility = mode === "variant" ? "visible" : "hidden";
      render();
    }));
  root.querySelectorAll<HTMLElement>(".bfilter, .btile[data-f]").forEach((el) =>
    el.addEventListener("click", () => {
      filter = (el.dataset.f as typeof filter) || "all"; mode = "variant"; dfilter = null;
      syncChips();
      (root.querySelector("#bfilters") as HTMLElement).style.visibility = "visible";
      render();
    }));
  root.querySelectorAll<HTMLElement>(".dcat").forEach((el) =>
    el.addEventListener("click", () => {
      const k = el.dataset.dcat!;
      dfilter = dfilter === k ? null : k;                       // toggle
      mode = "variant";
      const grp = dfilter ? DCAT[dfilter].grp : null;
      filter = grp === "fp" || grp === "fn" ? grp : "all";
      syncChips();
      (root.querySelector("#bfilters") as HTMLElement).style.visibility = "visible";
      render();
    }));
  search.addEventListener("input", () => { query = search.value.trim(); render(); });

  render();
}
