/** Downloads page: the genome-wide AlphaMissense-islands BED resource. */

interface Build { label: string; file: string; size: string; islands: string; }

const BUILDS: Build[] = [
  { label: "GRCh38 / hg38", file: "am_islands_GRCh38.bed.gz", size: "296 KB", islands: "13,391" },
  { label: "GRCh37 / hg19", file: "am_islands_GRCh37.bed.gz", size: "291 KB", islands: "13,183" },
];

const COLUMNS: [string, string][] = [
  ["chrom", "Chromosome (e.g. chr1)"],
  ["start / end", "0-based, half-open genomic span of the island (BED convention)"],
  ["transcript_id", "Ensembl transcript the island is defined on (AlphaMissense isoform)"],
  ["uniprot_id", "UniProt accession of the protein"],
  ["island_mean_score", "Mean AlphaMissense pathogenicity across the island's residues"],
  ["island_length", "Island length in residues"],
  ["residue_start / residue_end", "1-based protein residue range of the island"],
];

export function downloadsHTML(): string {
  const cards = BUILDS.map((b) => `
    <div class="dlcard">
      <div class="dlcard-h">
        <span class="dlcard-build">${b.label}</span>
        <span class="dlcard-meta">${b.islands} islands · BED · gzip · ${b.size}</span>
      </div>
      <a class="dlbtn" href="/downloads/${b.file}" download>
        <svg viewBox="0 0 24 24" width="18" height="18" aria-hidden="true" fill="none"
             stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
          <path d="M12 3v12"/><path d="M7 12l5 5 5-5"/><path d="M4 21h16"/></svg>
        Download <span class="mono">${b.file}</span>
      </a>
    </div>`).join("");

  const cols = COLUMNS.map(([c, d]) =>
    `<tr><td class="mono">${c}</td><td>${d}</td></tr>`).join("");

  return `
    <div class="page">
      <p class="eyebrow">Community resource</p>
      <h1>Download the <em>AlphaMissense islands</em></h1>
      <p class="page-lead">Genome-wide regional missense hotspots derived from AlphaMissense, as a
        tab-delimited BED file for both human assemblies. Small, dependency-free, and ready to
        intersect with a VCF (<span class="mono">bedtools intersect</span>) or load in a genome
        browser.</p>

      <div class="dlgrid">${cards}</div>

      <h2>What's inside</h2>
      <p>One row per island; nine tab-separated columns with a <span class="mono">#</span> header line.</p>
      <div class="dltable-wrap">
        <table class="dltable"><thead><tr><th>Column</th><th>Description</th></tr></thead>
          <tbody>${cols}</tbody></table>
      </div>

      <h2>How the islands were built</h2>
      <p>Per transcript, the per-residue mean AlphaMissense score is median-smoothed (window 5);
        contiguous runs above the AlphaMissense "likely-pathogenic" boundary (&gt; 0.564) of at
        least 35 residues are kept as islands. The generation script
        (<span class="mono">computeAlphaMissenseIslands.py</span>) and the full method are
        described on the <a href="?page=methods" data-methods class="ilink">Methods</a> page.</p>

      <h2>License &amp; attribution</h2>
      <p>These islands are a derived work of <b>AlphaMissense</b> (DeepMind), distributed under
        <b>CC BY 4.0</b>; the BED files are released under the same terms. If you use them, please
        cite AlphaMissense (Cheng <em>et al.</em>, <em>Science</em> 2023) and this resource:</p>
      <pre class="dlcite">AlphaMissense Islands — a genome-wide regional missense-hotspot resource
for ACMG PM1 (SeqOne / Archipelago). Derived from AlphaMissense
(Cheng et al., Science 2023, CC BY 4.0). archipelago2.vercel.app</pre>
      <p class="page-note">A versioned, citable release (Zenodo DOI) is planned; until then please
        note the download date. Found an issue with the data?
        <a href="?page=feedback" data-feedback class="ilink">Send us feedback →</a></p>
    </div>`;
}
