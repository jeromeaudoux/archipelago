/** PM1 benchmark explorer: the SeqOne engine's island-PM1 vs ClinGen-curated PM1
 *  (eRepo), matched by ClinVar ID — the same computation as the published benchmark. */
import { parseVariant, renderPM1 } from "./pm1";
import type { Bundle } from "./viewer";

interface Row {
  g: string; v: string; c: string; hgvs: string; dis: string; panel: string; link: string;
  clinvar: string; assertion: string; met: string;
  truth: boolean; tstr: string; pred: boolean; pstr: string; cat: "tp" | "fp" | "fn";
}
interface Data {
  summary: { tp: number; fp: number; fn: number; tn: number;
    precision: number; recall: number; f1: number; matched: number };
  variants: Row[];
}
type LoadBundle = (gene: string) => Promise<Bundle>;

const CAT = {
  tp: { label: "True positive", cls: "tp" },
  fp: { label: "False positive", cls: "fp" },
  fn: { label: "False negative", cls: "fn" },
};

export function renderBenchmark(root: HTMLElement, data: Data, loadBundle: LoadBundle): void {
  const s = data.summary;
  let mode: "variant" | "gene" = "variant";
  let filter: "all" | "tp" | "fp" | "fn" = "all";
  let query = "";
  let gsort = { key: "n", dir: -1 };

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
        <input id="bench-search" class="bsearch" type="text" placeholder="Filter by gene, variant or disease…" spellcheck="false" />
      </div>

      <div class="btable-wrap"><table class="btable"><thead id="bench-head"></thead><tbody id="bench-body"></tbody></table></div>
      <p id="bench-count" class="bench-count"></p>
      <p class="bench-caveat">Prediction = the engine's <span class="mono">acmg_tags_v2</span> PM1 state
        (AlphaMissense-island overlap). Ground truth = PM1 (any strength) in the eRepo “Applied Evidence
        Codes (Met)”, matched by ClinVar ID and compared strength-independently, exactly as in
        <span class="mono">benchmark_acmg.py</span>.</p>
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
  const geneLink = (g: string) => `<a href="?gene=${g}" data-gene-link="${g}" class="blink bg">${g}</a>`;

  function renderVariantTable(): void {
    head.innerHTML = `<tr><th>Gene</th><th>Variant</th><th>Disease</th>
      <th>ClinGen PM1</th><th>Engine PM1</th><th>Result</th></tr>`;
    const q = query.toLowerCase();
    view = data.variants.filter((r) =>
      (filter === "all" || r.cat === filter)
      && (!q || r.g.toLowerCase().includes(q) || (r.hgvs + r.v).toLowerCase().includes(q) || r.dis.toLowerCase().includes(q)));
    body.innerHTML = view.slice(0, 3000).map((r, i) => `<tr data-i="${i}">
      <td>${geneLink(r.g)}</td>
      <td class="mono bv">${r.hgvs || r.v}</td>
      <td class="bdis" title="${r.dis}">${r.dis || "—"}</td>
      <td>${badge(r.tstr, r.truth)}</td>
      <td>${badge(r.pstr, r.pred)}</td>
      <td><span class="rescat ${CAT[r.cat].cls}">${CAT[r.cat].label}</span></td>
    </tr>`).join("");
    countEl.textContent = `${view.length.toLocaleString()} variant${view.length === 1 ? "" : "s"}`
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
    const strip = `<div class="bmodal-strip">
      <span class="rescat ${CAT[r.cat].cls}">${CAT[r.cat].label}</span>
      <span class="bstrip-item">ClinGen: <b>${r.truth ? (r.tstr || "PM1") : "PM1 not applied"}</b> · ${r.assertion || "—"}</span>
      <span class="bstrip-item">${r.panel || ""}</span>
      <span class="bstrip-codes">${(r.met || "").split(",").map((c) => c.trim()).filter(Boolean)
        .map((c) => `<span class="code ${c.startsWith("PM1") ? "pm1" : ""}">${c}</span>`).join("")}</span>
      ${r.link ? `<a class="blink" href="${r.link}" target="_blank" rel="noopener">eRepo ↗</a>` : ""}
    </div>`;
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
      if (gtr) { mode = "variant"; query = gtr.dataset.gene!; search.value = query;
        root.querySelectorAll(".bmodebtn").forEach((b) => b.classList.toggle("on", (b as HTMLElement).dataset.m === "variant"));
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

  root.querySelectorAll<HTMLElement>(".bmodebtn").forEach((el) =>
    el.addEventListener("click", () => {
      mode = (el.dataset.m as typeof mode) || "variant";
      root.querySelectorAll(".bmodebtn").forEach((b) => b.classList.toggle("on", b === el));
      (root.querySelector("#bfilters") as HTMLElement).style.visibility = mode === "variant" ? "visible" : "hidden";
      render();
    }));
  root.querySelectorAll<HTMLElement>(".bfilter, .btile[data-f]").forEach((el) =>
    el.addEventListener("click", () => {
      filter = (el.dataset.f as typeof filter) || "all"; mode = "variant";
      root.querySelectorAll(".bmodebtn").forEach((b) => b.classList.toggle("on", (b as HTMLElement).dataset.m === "variant"));
      root.querySelectorAll(".bfilter").forEach((f) => f.classList.toggle("on", (f as HTMLElement).dataset.f === filter));
      (root.querySelector("#bfilters") as HTMLElement).style.visibility = "visible";
      render();
    }));
  search.addEventListener("input", () => { query = search.value.trim(); render(); });

  render();
}
