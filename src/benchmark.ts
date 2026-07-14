/** PM1 benchmark explorer: ClinGen-curated PM1 (eRepo) vs Archipelago island-PM1. */

interface Row {
  g: string; v: string; p: number; dis: string; panel: string; link: string;
  clinvar: string; truth: boolean; tstr: string; pred: boolean; pstr: string;
  isl: [number, number] | null; cat: "tp" | "fp" | "fn";
}
interface Data {
  summary: {
    tp: number; fp: number; fn: number; tn: number;
    precision: number; recall: number; f1: number;
    ref_mismatch: number; skipped_non_missense: number; skipped_no_gene: number;
  };
  variants: Row[];
}

const CAT = {
  tp: { label: "True positive", cls: "tp" },
  fp: { label: "False positive", cls: "fp" },
  fn: { label: "False negative", cls: "fn" },
};

export function renderBenchmark(root: HTMLElement, data: Data): void {
  const s = data.summary;
  const curated = s.tp + s.fn;
  let filter: "all" | "tp" | "fp" | "fn" = "all";
  let query = "";

  root.innerHTML = `
    <div class="bench-page">
      <p class="eyebrow">PM1 benchmark · ClinGen Evidence Repository</p>
      <h1>Explore <em>ClinGen-curated PM1</em></h1>
      <p class="bench-lead">Every eRepo missense variant that maps onto an AlphaMissense isoform,
        classified by whether ClinGen expert panels applied <b>PM1</b> (ground truth) and whether the
        variant's residue falls in an <b>Archipelago AlphaMissense island</b> (prediction).</p>

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
            <th>ClinGen PM1</th><th>Archipelago PM1</th><th>Result</th><th></th>
          </tr></thead>
          <tbody id="bench-body"></tbody>
        </table>
      </div>
      <p id="bench-count" class="bench-count"></p>
      <p class="bench-caveat">Ground truth = PM1 in the eRepo “Applied Evidence Codes (Met)”. Prediction =
        residue in an AlphaMissense island (island-membership rule, as in the app). This app-consistent
        recomputation differs from the full SeqOne engine benchmark (which adds ClinVar enrichment and
        reaches ~0.70 precision). Excluded: ${s.skipped_non_missense.toLocaleString()} non-missense,
        ${s.ref_mismatch.toLocaleString()} isoform-discordant, ${s.skipped_no_gene.toLocaleString()} genes
        absent from AlphaMissense. Curated PM1 variants shown: ${curated.toLocaleString()} (${s.tp} caught, ${s.fn} missed).</p>
    </div>`;

  const body = root.querySelector<HTMLElement>("#bench-body")!;
  const countEl = root.querySelector<HTMLElement>("#bench-count")!;

  function rowHTML(r: Row): string {
    const pm1link = `?pm1=${encodeURIComponent(`${r.g}:p.${r.v}`)}`;
    const badge = (str: string, on: boolean) =>
      `<span class="pm1badge ${on ? "on" : "off"}">${str || (on ? "PM1" : "—")}</span>`;
    return `<tr>
      <td><a href="?gene=${r.g}" data-gene-link="${r.g}" class="blink">${r.g}</a></td>
      <td><a href="${pm1link}" data-pm1-link="${r.g}:p.${r.v}" class="blink mono">${r.v}</a></td>
      <td class="bdis" title="${r.dis}">${r.dis || "—"}</td>
      <td>${badge(r.tstr, r.truth)}</td>
      <td>${badge(r.pstr, r.pred)}</td>
      <td><span class="rescat ${CAT[r.cat].cls}">${CAT[r.cat].label}</span></td>
      <td>${r.link ? `<a href="${r.link}" target="_blank" rel="noopener" class="blink">eRepo↗</a>` : ""}</td>
    </tr>`;
  }

  function apply(): void {
    const q = query.toLowerCase();
    const rows = data.variants.filter((r) =>
      (filter === "all" || r.cat === filter)
      && (!q || r.g.toLowerCase().includes(q) || r.v.toLowerCase().includes(q) || r.dis.toLowerCase().includes(q)));
    body.innerHTML = rows.slice(0, 2000).map(rowHTML).join("");
    countEl.textContent = `${rows.length.toLocaleString()} variant${rows.length === 1 ? "" : "s"}`
      + (rows.length > 2000 ? " (showing first 2,000)" : "");
  }

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
