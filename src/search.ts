/** Reusable gene autocomplete bound to an <input> + suggestion <ul>. */

export interface GeneEntry { symbol: string; uniprot: string; len: number; }

export function attachAutocomplete(
  input: HTMLInputElement,
  listEl: HTMLUListElement,
  getIndex: () => GeneEntry[],
  onSelect: (symbol: string) => void,
): void {
  let matches: GeneEntry[] = [];
  let active = -1;

  const upper = (s: string) => s.toUpperCase();

  function rank(q: string): GeneEntry[] {
    const idx = getIndex();
    if (!q) return [];
    const Q = upper(q);
    const starts: GeneEntry[] = [];
    const contains: GeneEntry[] = [];
    for (const g of idx) {
      const s = g.symbol;
      if (s === Q) starts.unshift(g);
      else if (s.startsWith(Q)) starts.push(g);
      else if (s.includes(Q)) contains.push(g);
      if (starts.length + contains.length > 400) break;
    }
    return [...starts, ...contains].slice(0, 30);
  }

  function render(): void {
    if (!matches.length) {
      if (input.value.trim() && getIndex().length) {
        listEl.innerHTML = `<li class="empty">No gene matching “${input.value.trim()}”</li>`;
        listEl.hidden = false;
      } else {
        listEl.hidden = true;
      }
      return;
    }
    listEl.innerHTML = matches
      .map((g, i) =>
        `<li role="option" data-sym="${g.symbol}" aria-selected="${i === active}">
           <span class="sym">${g.symbol}</span>
           <span class="up mono">${g.uniprot}</span>
           <span class="len mono">${g.len} aa</span>
         </li>`)
      .join("");
    listEl.hidden = false;
  }

  function close(): void { listEl.hidden = true; active = -1; }

  input.addEventListener("input", () => {
    matches = rank(input.value.trim());
    active = matches.length ? 0 : -1;
    render();
  });

  input.addEventListener("keydown", (e) => {
    if (listEl.hidden && e.key !== "Enter") return;
    if (e.key === "ArrowDown") { e.preventDefault(); active = Math.min(active + 1, matches.length - 1); render(); }
    else if (e.key === "ArrowUp") { e.preventDefault(); active = Math.max(active - 1, 0); render(); }
    else if (e.key === "Enter") {
      e.preventDefault();
      const pick = matches[active];
      if (pick) { onSelect(pick.symbol); close(); }
      else if (input.value.trim()) { onSelect(upper(input.value.trim())); close(); }
    } else if (e.key === "Escape") { close(); }
  });

  listEl.addEventListener("mousedown", (e) => {
    const li = (e.target as HTMLElement).closest("li[data-sym]") as HTMLElement | null;
    if (li) { e.preventDefault(); onSelect(li.dataset.sym!); close(); }
  });

  input.addEventListener("blur", () => setTimeout(close, 120));
}
