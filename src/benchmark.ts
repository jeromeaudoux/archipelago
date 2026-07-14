/** PM1 benchmark explorer: the SeqOne engine's island-PM1 vs ClinGen-curated PM1
 *  (eRepo), matched by ClinVar ID — the same computation as the published benchmark. */

interface Row {
  g: string; v: string; hgvs: string; dis: string; panel: string; link: string;
  clinvar: string; assertion: string; met: string;
  truth: boolean; tstr: string; pred: boolean; pstr: string; cat: "tp" | "fp" | "fn";
}
interface Data {
  summary: { tp: number; fp: number; fn: number; tn: number;
    precision: number; recall: number; f1: number; matched: number };
  variants: Row[];
}

const CAT = {
  tp: { label: "True positive", cls: "tp" },
  fp: { label: "False positive", cls: "fp" },
  fn: { label: "False negative", cls: "fn" },
};

export function renderBenchmark(root: HTMLElement, data: Data): void {
  const s = data.summary;
  let filter: "all" | "tp" | "fp" | "fn" = "all";
  let query = "";

  root.innerHTML = `
    <div class="bench-page">
      <p class="eyebrow">PM1 benchmark · ClinGen Evidence Repository</p>
      <h1>Explore <em>ClinGen-curated PM1</em></h1>
      <p class="bench-lead">The SeqOne engine's island-based <b>PM1</b> call versus ClinGen expert-panel
        curation, matched by ClinVar ID over ${s.matched.toLocaleString()} eRepo variants — the same
        computation as the published benchmark (precision ${(s.precision * 100).toFixed(0)}%). Click a row
        for the full curation.</p>

      <div class="bench-tiles">
        <button class="btile" data-f="tp"><span class="bt-n">${s.tp}</span><span class="bt-l">True positive</span></button>
        <button class="btile" data-f="fp"><span class="bt-n">${s.fp}</span><span class="bt-l">False positive</span></button>
        <button class="btile" data-f="fn"><span class="bt-n">${s.fn}</span><span class="bt-l">False negative</span></button>
        <div class="btile stat"><span class="bt-n">${(s.precision * 100).toFixed(0)}%</span><span class="bt-l">Precision</span></div>
        <div class="btile stat"><span class="bt-n">${(s.recall * 100).toFixed(0)}%</span><span class="bt-l">Recall</span></div>
        <div class="btile stat"><span class="bt-n">${s.f1.toFixed(2)}</span><span class="bt-l">F1</span></div>
      </div>

      <div class="bench-controls">
        <div class="bfilters">
          <button class="bfilter on" data-f="all">All (${data.variants.length})</button>
          <button class="bfilter" data-f="tp">TP ${s.tp}</button>
          <button class="bfilter" data-f="fp">FP ${s.fp}</button>
          <button class="bfilter" data-f="fn">FN ${s.fn}</button>
        </div>
        <input id="bench-search" class="bsearch" type="text" placeholder="Filter by gene, variant or disease…" spellcheck="false" />
      </div>

      <div class="btable-wrap">
        <table class="btable">
          <thead><tr>
            <th>Gene</th><th>Variant</th><th>Disease</th>
            <th>ClinGen PM1</th><th>Engine PM1</th><th>Result</th>
          </tr></thead>
          <tbody id="bench-body"></tbody>
        </table>
      </div>
      <p id="bench-count" class="bench-count"></p>
      <p class="bench-caveat">Prediction = the engine's <span class="mono">acmg_tags_v2</span> PM1 state
        (AlphaMissense-island overlap, scored with <span class="mono">--am-islands</span>). Ground truth =
        PM1 (any strength) in the eRepo “Applied Evidence Codes (Met)”. Matched by ClinVar ID; PM1 is
        compared strength-independently, exactly as in <span class="mono">benchmark_acmg.py</span>.</p>
      <div id="bench-modal"></div>
    </div>`;

  const body = root.querySelector<HTMLElement>("#bench-body")!;
  const countEl = root.querySelector<HTMLElement>("#bench-count")!;
  const modalHost = root.querySelector<HTMLElement>("#bench-modal")!;
  let view: Row[] = [];

  const badge = (str: string, on: boolean) =>
    `<span class="pm1badge ${on ? "on" : "off"}">${str || (on ? "PM1" : "not applied")}</span>`;

  function rowHTML(r: Row, i: number): string {
    return `<tr data-i="${i}">
      <td class="bg">${r.g}</td>
      <td class="mono bv">${r.hgvs || r.v}</td>
      <td class="bdis" title="${r.dis}">${r.dis || "—"}</td>
      <td>${badge(r.tstr, r.truth)}</td>
      <td>${badge(r.pstr, r.pred)}</td>
      <td><span class="rescat ${CAT[r.cat].cls}">${CAT[r.cat].label}</span></td>
    </tr>`;
  }

  function apply(): void {
    const q = query.toLowerCase();
    view = data.variants.filter((r) =>
      (filter === "all" || r.cat === filter)
      && (!q || r.g.toLowerCase().includes(q) || (r.hgvs + r.v).toLowerCase().includes(q) || r.dis.toLowerCase().includes(q)));
    body.innerHTML = view.slice(0, 3000).map(rowHTML).join("");
    countEl.textContent = `${view.length.toLocaleString()} variant${view.length === 1 ? "" : "s"}`
      + (view.length > 3000 ? " (showing first 3,000)" : "");
  }

  function openModal(r: Row): void {
    const pm1link = r.v ? `?pm1=${encodeURIComponent(`${r.g}:p.${r.v}`)}` : `?gene=${r.g}`;
    modalHost.innerHTML = `<div class="modal-backdrop" data-close>
      <div class="modal bmodal" role="dialog" aria-modal="true">
        <div class="modal-head">
          <h3><span class="mono" style="color:var(--accent)">${r.g}</span> ${r.hgvs || r.v}</h3>
          <button class="modal-x" data-close aria-label="Close">×</button>
        </div>
        <div class="bmodal-cat"><span class="rescat ${CAT[r.cat].cls}">${CAT[r.cat].label}</span></div>
        <div class="bmodal-grid">
          <div class="fact"><span>ClinGen PM1</span><b>${r.truth ? (r.tstr || "PM1") : "not applied"}</b></div>
          <div class="fact"><span>Engine PM1</span><b>${r.pred ? (r.pstr || "PM1") + " (met)" : "not met"}</b></div>
          <div class="fact"><span>ClinGen assertion</span><b>${r.assertion || "—"}</b></div>
          <div class="fact"><span>Disease</span><b>${r.dis || "—"}</b></div>
          <div class="fact"><span>Expert panel</span><b>${r.panel || "—"}</b></div>
          <div class="fact"><span>ClinVar</span><b>${r.clinvar || "—"}</b></div>
        </div>
        <div class="bmodal-codes"><span>ClinGen applied criteria</span>
          <div class="codes">${(r.met || "").split(",").map((c) => c.trim()).filter(Boolean)
            .map((c) => `<span class="code ${c.startsWith("PM1") ? "pm1" : ""}">${c}</span>`).join("") || "—"}</div>
        </div>
        <div class="bmodal-actions">
          <a class="ghost" href="${pm1link}" data-gene-link="${r.g}">Open ${r.g} in Archipelago →</a>
          ${r.link ? `<a class="ghost" href="${r.link}" target="_blank" rel="noopener">ClinGen eRepo ↗</a>` : ""}
        </div>
      </div>
    </div>`;
  }
  function closeModal(): void { modalHost.innerHTML = ""; }

  body.addEventListener("click", (e) => {
    const tr = (e.target as HTMLElement).closest<HTMLElement>("tr[data-i]");
    if (tr) openModal(view[parseInt(tr.dataset.i!, 10)]);
  });
  modalHost.addEventListener("click", (e) => {
    const t = e.target as HTMLElement;
    if (t.closest("[data-close]")) closeModal();
    else if (t.closest("[data-gene-link]")) closeModal();  // let main.ts handle navigation
  });
  root.addEventListener("keydown", (e) => { if ((e as KeyboardEvent).key === "Escape") closeModal(); });

  root.querySelectorAll<HTMLElement>(".bfilter, .btile[data-f]").forEach((el) =>
    el.addEventListener("click", () => {
      filter = (el.dataset.f as typeof filter) || "all";
      root.querySelectorAll(".bfilter").forEach((f) =>
        f.classList.toggle("on", (f as HTMLElement).dataset.f === filter));
      apply();
    }));
  const search = root.querySelector<HTMLInputElement>("#bench-search")!;
  search.addEventListener("input", () => { query = search.value.trim(); apply(); });

  apply();
}
