import "./styles.css";
import { attachAutocomplete, type GeneEntry } from "./search";
import { renderViewer, type Bundle } from "./viewer";
import { methodsHTML } from "./methods";
import { parseVariant, looksLikeVariant, renderPM1 } from "./pm1";
import { renderBenchmark } from "./benchmark";

const EXAMPLES = ["KIF1A", "TP53", "BRCA1", "SCN2A", "PTEN", "SCN1A"];
const EXAMPLE_VARIANTS = ["KIF1A:p.Arg350Gln", "KIF1A:c.760G>A"];

let INDEX: GeneEntry[] = [];
const bundleCache = new Map<string, Bundle>();

const $ = <T extends HTMLElement = HTMLElement>(id: string) =>
  document.getElementById(id) as T;

const landing = $("landing");
const viewer = $("viewer");
const about = $("about");
const viewerRoot = $("viewer-root");

// ---- theme ----
function toggleTheme(): void {
  const cur = document.documentElement.getAttribute("data-theme")
    || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
  const next = cur === "dark" ? "light" : "dark";
  document.documentElement.setAttribute("data-theme", next);
  localStorage.setItem("am-theme", next);
}
function initTheme(): void {
  const saved = localStorage.getItem("am-theme");
  if (saved) document.documentElement.setAttribute("data-theme", saved);
}

// ---- data ----
async function loadIndex(): Promise<void> {
  try {
    const res = await fetch("/data/index.json");
    if (!res.ok) throw new Error(String(res.status));
    INDEX = await res.json();
    $("index-status").textContent = `${INDEX.length.toLocaleString()} genes indexed`;
  } catch {
    $("index-status").textContent =
      "Gene index not found — run the precompute pipeline (npm run precompute).";
  }
}

async function fetchBundle(symbol: string): Promise<Bundle> {
  const cached = bundleCache.get(symbol);
  if (cached) return cached;
  const res = await fetch(`/data/genes/${encodeURIComponent(symbol)}.json.gz`);
  if (!res.ok) throw new Error(res.status === 404 ? "notfound" : String(res.status));
  const buf = await res.arrayBuffer();
  const bytes = new Uint8Array(buf);
  let text: string;
  if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
    const stream = new Response(buf).body!.pipeThrough(new DecompressionStream("gzip"));
    text = await new Response(stream).text();
  } else {
    text = new TextDecoder().decode(buf);
  }
  const bundle = JSON.parse(text) as Bundle;
  bundleCache.set(symbol, bundle);
  return bundle;
}

// ---- views ----
function show(view: "landing" | "viewer" | "about"): void {
  landing.hidden = view !== "landing";
  viewer.hidden = view !== "viewer";
  about.hidden = view !== "about";
}
function showLanding(): void { show("landing"); document.title = "Archipelago — AlphaMissense island explorer"; }

function showMethods(push = true): void {
  if (push) history.pushState({ page: "methods" }, "", "?page=methods");
  $("about-root").innerHTML = methodsHTML();
  show("about");
  document.title = "Methods — Archipelago";
  window.scrollTo(0, 0);
}

let benchData: Parameters<typeof renderBenchmark>[1] | null = null;
async function showBenchmark(push = true): Promise<void> {
  if (push) history.pushState({ page: "benchmark" }, "", "?page=benchmark");
  show("about");
  document.title = "PM1 benchmark — Archipelago";
  window.scrollTo(0, 0);
  const root = $("about-root");
  root.innerHTML = `<div class="state"><div class="spinner"></div><p>Loading PM1 benchmark…</p></div>`;
  try {
    if (!benchData) benchData = await (await fetch("/data/pm1_benchmark.json")).json();
    renderBenchmark(root, benchData!);
  } catch {
    root.innerHTML = `<div class="state"><h2>Benchmark data not found</h2>
      <p>Run <span class="mono">python3 scripts/build_pm1_benchmark.py</span>.</p></div>`;
  }
}

async function showGene(symbol: string, push = true): Promise<void> {
  symbol = symbol.toUpperCase().trim();
  if (!symbol) return;
  if (push) history.pushState({ gene: symbol }, "", `?gene=${symbol}`);
  show("viewer");
  document.title = `${symbol} · Archipelago`;
  viewerRoot.innerHTML = `<div class="state"><div class="spinner"></div><p>Loading ${symbol}…</p></div>`;
  try {
    renderViewer(viewerRoot, await fetchBundle(symbol));
  } catch (err) {
    viewerRoot.innerHTML = notFound(symbol, (err as Error).message);
  }
}

function esc(s: string): string {
  return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

function showVariantError(raw: string, push = true): void {
  if (push) history.pushState({ pm1: raw }, "", `?pm1=${encodeURIComponent(raw)}`);
  show("viewer");
  document.title = "Variant not recognized · Archipelago";
  viewerRoot.innerHTML = `<div class="state">
    <h2>Couldn't read that variant</h2>
    <p>“${esc(raw)}” isn't a variant I can parse. Use <b>GENE:p.…</b> or <b>GENE:c.…</b> —
       e.g. <span class="mono">KIF1A:p.Arg350Gln</span>, <span class="mono">KIF1A:p.R350Q</span>,
       or <span class="mono">KIF1A:c.760G&gt;A</span>.</p>
    <p><button class="chip" data-home>← Back to search</button></p>
  </div>`;
}

async function showPM1(raw: string, push = true): Promise<void> {
  const pv = parseVariant(raw);
  if (!pv) { showVariantError(raw, push); return; }
  if (push) history.pushState({ pm1: raw }, "", `?pm1=${encodeURIComponent(raw)}`);
  show("viewer");
  document.title = `${pv.gene} ${pv.kind}. · PM1 · Archipelago`;
  viewerRoot.innerHTML = `<div class="state"><div class="spinner"></div><p>Predicting PM1 for ${pv.gene}…</p></div>`;
  try {
    renderPM1(viewerRoot, pv, await fetchBundle(pv.gene));
  } catch (err) {
    viewerRoot.innerHTML = notFound(pv.gene, (err as Error).message);
  }
}

function notFound(symbol: string, msg: string): string {
  return `<div class="state">
    <h2>${msg === "notfound" ? `No data for “${symbol}”` : "Couldn't load gene"}</h2>
    <p>${msg === "notfound"
      ? "This gene isn't in the precomputed set — AlphaMissense doesn't cover every protein."
      : "Something went wrong fetching the gene bundle."}</p>
    <p><button class="chip" data-home>← Back to search</button></p>
  </div>`;
}

/** Route a free-text query to either a variant PM1 page or a gene map. */
function go(query: string): void {
  const q = query.trim();
  if (parseVariant(q)) showPM1(q);
  else if (looksLikeVariant(q)) showVariantError(q);
  else showGene(q);
}

function route(): void {
  const params = new URLSearchParams(location.search);
  const pm1 = params.get("pm1");
  const gene = params.get("gene");
  const page = params.get("page");
  if (pm1) showPM1(pm1, false);
  else if (page === "methods") showMethods(false);
  else if (page === "benchmark") showBenchmark(false);
  else if (gene) showGene(gene, false);
  else showLanding();
}

// ---- wire up ----
function init(): void {
  initTheme();

  const chips = $("example-chips");
  chips.innerHTML = EXAMPLES
    .map((g) => `<button class="chip" data-gene="${g}"><span class="mono">${g}</span></button>`)
    .join("")
    + EXAMPLE_VARIANTS.map((v) =>
      `<button class="chip chip-variant" data-variant="${v}" title="Predict PM1 for this variant">`
      + `<span class="mono">${v}</span></button>`).join("");

  attachAutocomplete($("search-input"), $("suggest") as HTMLUListElement, () => INDEX, go);
  attachAutocomplete($("search-input-top"), $("suggest-top") as HTMLUListElement, () => INDEX, go);

  document.addEventListener("click", (e) => {
    const t = e.target as HTMLElement;
    const gchip = t.closest<HTMLElement>(".chip[data-gene]");
    if (gchip) { showGene(gchip.dataset.gene!); return; }
    const vchip = t.closest<HTMLElement>("[data-variant]");
    if (vchip) { showPM1(vchip.dataset.variant!); return; }
    const glink = t.closest<HTMLElement>("[data-gene-link]");
    if (glink) { e.preventDefault(); showGene(glink.dataset.geneLink!); return; }
    const plink = t.closest<HTMLElement>("[data-pm1-link]");
    if (plink) { e.preventDefault(); showPM1(plink.dataset.pm1Link!); return; }
    if (t.closest("[data-home]")) { e.preventDefault(); showLanding(); history.pushState({}, "", "/"); return; }
    if (t.closest("[data-methods]")) { e.preventDefault(); showMethods(); return; }
    if (t.closest("[data-benchmark]")) { e.preventDefault(); showBenchmark(); return; }
    if (t.closest(".theme-toggle")) { toggleTheme(); return; }
  });

  window.addEventListener("popstate", route);
  loadIndex();
  route();
}

init();
