import "./styles.css";
import { attachAutocomplete, type GeneEntry } from "./search";
import { renderViewer, type Bundle } from "./viewer";
import { methodsHTML } from "./methods";

const EXAMPLES = ["KIF1A", "TP53", "BRCA1", "SCN2A", "PTEN", "SCN1A"];

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

async function showGene(symbol: string, push = true): Promise<void> {
  symbol = symbol.toUpperCase().trim();
  if (!symbol) return;
  if (push) history.pushState({ gene: symbol }, "", `?gene=${symbol}`);
  show("viewer");
  document.title = `${symbol} · Archipelago`;
  viewerRoot.innerHTML = `<div class="state"><div class="spinner"></div><p>Loading ${symbol}…</p></div>`;
  try {
    const bundle = await fetchBundle(symbol);
    renderViewer(viewerRoot, bundle);
  } catch (err) {
    const msg = (err as Error).message;
    viewerRoot.innerHTML = `<div class="state">
      <h2>${msg === "notfound" ? `No data for “${symbol}”` : "Couldn't load gene"}</h2>
      <p>${msg === "notfound"
        ? "This gene isn't in the precomputed set — it may have no AlphaMissense transcript (AlphaMissense doesn't cover every protein)."
        : "Something went wrong fetching the gene bundle."}</p>
      <p><button class="chip" data-home>← Back to search</button></p>
    </div>`;
  }
}

function route(): void {
  const params = new URLSearchParams(location.search);
  const gene = params.get("gene");
  const page = params.get("page");
  if (page === "methods") showMethods(false);
  else if (gene) showGene(gene, false);
  else showLanding();
}

// ---- wire up ----
function init(): void {
  initTheme();

  const chips = $("example-chips");
  chips.innerHTML = EXAMPLES
    .map((g) => `<button class="chip" data-gene="${g}"><span class="mono">${g}</span></button>`)
    .join("");

  attachAutocomplete($("search-input"), $("suggest") as HTMLUListElement, () => INDEX, (s) => showGene(s));
  attachAutocomplete($("search-input-top"), $("suggest-top") as HTMLUListElement, () => INDEX, (s) => showGene(s));

  // event delegation for nav + chips + theme
  document.addEventListener("click", (e) => {
    const t = e.target as HTMLElement;
    const chip = t.closest<HTMLElement>(".chip[data-gene]");
    if (chip) { showGene(chip.dataset.gene!); return; }
    if (t.closest("[data-home]")) { e.preventDefault(); showLanding(); history.pushState({}, "", "/"); return; }
    if (t.closest("[data-methods]")) { e.preventDefault(); showMethods(); return; }
    if (t.closest(".theme-toggle")) { toggleTheme(); return; }
  });

  window.addEventListener("popstate", route);
  loadIndex();
  route();
}

init();
