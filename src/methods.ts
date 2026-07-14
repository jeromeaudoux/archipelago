/** Methods / "How the islands work" page content. */

export function methodsHTML(): string {
  return `
  <article class="methods">
    <p class="eyebrow">Methods</p>
    <h1>How AlphaMissense <em>islands</em> work</h1>
    <p class="lead">
      An <b>AlphaMissense island</b> is a contiguous stretch of protein residues that
      AlphaMissense predicts to be consistently intolerant to missense change. Because
      AlphaMissense is derived from a protein language model and structural context — not
      from clinical databases — islands are an <b>orthogonal, non-circular</b> signal for
      the ACMG <b>PM1</b> criterion (“located in a mutational hot spot / critical
      functional domain”), unlike the traditional approach of counting ClinVar variants in
      a fixed window.
    </p>

    <h2>How an island is defined</h2>
    <ol class="steps">
      <li>
        <span class="step-k">Per-residue score</span>
        For each residue, take the <b>mean AlphaMissense pathogenicity</b> across all
        possible substitutions at that position.
      </li>
      <li>
        <span class="step-k">Disorder mask <span class="opt">optional</span></span>
        Residues with AlphaFold <b>pLDDT &lt; 50</b> are set to 0, so islands terminate in
        disordered regions where the structural context — and thus the prediction — is
        unreliable.
      </li>
      <li>
        <span class="step-k">Smoothing</span>
        A <b>median filter (window = 5 residues)</b> removes isolated spikes/dips, yielding
        longer, more robust segments.
      </li>
      <li>
        <span class="step-k">Segmentation</span>
        Keep contiguous runs where the smoothed score exceeds <b>0.564</b> (the
        AlphaMissense “ambiguous” boundary). A run breaks at the threshold or at any gap in
        residue numbering, and runs shorter than <b>35 residues</b> are discarded.
      </li>
    </ol>

    <figure class="svgfig">
      <svg viewBox="0 0 720 250" role="img" aria-label="Island creation: per-residue AM score is smoothed, thresholded at 0.564, and contiguous runs of at least 35 residues become islands">
        <!-- detected island bands -->
        <rect x="118" y="24" width="140" height="150" fill="var(--isl-path)" opacity="0.16" rx="3" />
        <rect x="438" y="24" width="140" height="150" fill="var(--isl-path)" opacity="0.16" rx="3" />
        <!-- rejected short run -->
        <rect x="338" y="24" width="26" height="150" fill="none" stroke="var(--muted)" stroke-width="1" stroke-dasharray="3 3" rx="3" />
        <!-- axes -->
        <line x1="55" y1="174" x2="705" y2="174" stroke="var(--border)" stroke-width="1" />
        <!-- threshold -->
        <line x1="55" y1="94" x2="705" y2="94" stroke="var(--muted)" stroke-width="1" stroke-dasharray="5 4" />
        <text x="60" y="88" fill="var(--muted)" font-size="11" font-family="var(--font-mono)">threshold 0.564</text>
        <!-- raw (faint) -->
        <path d="M60 150 L82 128 L98 150 L120 84 L138 66 L160 52 L182 70 L210 58 L232 78 L255 88 L272 150 L285 150 L310 150 L330 150 L345 66 L360 74 L375 150 L400 150 L420 150 L440 84 L462 58 L490 46 L515 66 L540 54 L560 82 L575 92 L598 150 L610 150 L700 150"
          fill="none" stroke="var(--muted)" stroke-width="1" opacity="0.5" />
        <!-- smoothed -->
        <path d="M60 152 C80 150 96 140 120 92 S180 50 210 58 S262 84 285 150 L310 150 L332 150 L345 78 L360 80 L375 150 L400 150 L422 150 C446 150 462 66 490 52 S556 66 575 92 S610 150 700 150"
          fill="none" stroke="var(--accent)" stroke-width="2.4" stroke-linejoin="round" />
        <!-- labels -->
        <text x="188" y="192" fill="var(--isl-path)" font-size="11.5" font-weight="600" text-anchor="middle">island ✓ ≥35 aa</text>
        <text x="351" y="192" fill="var(--muted)" font-size="11" text-anchor="middle">✗ &lt;35 aa</text>
        <text x="508" y="192" fill="var(--isl-path)" font-size="11.5" font-weight="600" text-anchor="middle">island ✓ ≥35 aa</text>
        <text x="60" y="240" fill="var(--muted)" font-size="11">mean AM per residue &nbsp;→&nbsp; median-smooth &nbsp;→&nbsp; keep runs above threshold &nbsp;→&nbsp; length filter</text>
        <text x="700" y="240" fill="var(--muted)" font-size="11" text-anchor="end">residue →</text>
      </svg>
      <figcaption>Per-residue mean AlphaMissense (faint) is median-smoothed (blue); contiguous
        residues above 0.564 form islands, and runs shorter than 35 aa are discarded.</figcaption>
    </figure>

    <p class="src">Reference implementation: <code>computeAlphaMissenseIslands.py</code>
      (acmgscore-v2). Each island is emitted with its residue span, genomic coordinates and
      mean AM score — exactly the <code>islands</code> track you see in the viewer.</p>

    <h2>From islands to a PM1 prediction</h2>
    <p>
      At classification time a variant's residue is intersected with the island set. A
      variant that falls <b>inside an island</b> satisfies the AM-defined hot-spot condition;
      combining that with ClinVar pathogenic enrichment <i>within the island boundaries</i>
      gives a PM1 call that is both structurally grounded and clinically anchored. In this
      viewer the same islands are additionally <b>colored by their ClinVar overlap</b>
      (red = pathogenic-heavy, blue = benign-heavy), so you can see at a glance which
      predicted hot-spots are corroborated by real variants.
    </p>

    <h2>Grading PM1 by cumulative evidence</h2>
    <p>
      PM1 need not be all-or-nothing. Two independent lines of evidence can be combined and
      the criterion's <b>strength dialed</b> on the ACMG points scale (Supporting = 1,
      Moderate = 2, Strong = 4): the <b>AM island</b> (a structural/predictive hot-spot) and
      the <b>ClinVar pathogenic density inside that island</b> (empirical corroboration).
      The more the two agree — and the cleaner the island is of benign variation — the
      stronger the call.
    </p>
    <div class="grades">
      <table class="bench-table grades-table">
        <thead><tr><th>Evidence at the residue</th><th>PM1 strength</th></tr></thead>
        <tbody>
          <tr><td>Not in an AM island</td><td><span class="pill none">Not met</span></td></tr>
          <tr><td>In an AM island, little/no ClinVar corroboration</td><td><span class="pill sup">Supporting</span> <span class="pts">+1</span></td></tr>
          <tr><td>In an AM island <b>and</b> the island is pathogenic-enriched (several P/LP, no benign)</td><td><span class="pill mod">Moderate</span> <span class="pts">+2</span></td></tr>
          <tr><td>In a dense, well-established hot-spot island (many P/LP, zero benign)</td><td><span class="pill str">Strong</span> <span class="pts">+4</span></td></tr>
          <tr><td>Island contradicted by benign variation (B/LB-heavy)</td><td><span class="pill none">Downweight / not met</span></td></tr>
        </tbody>
      </table>
    </div>
    <p class="src">
      This is the framework Archipelago's data supports — the island track supplies the
      hot-spot boundary, and each island's red/blue ClinVar colouring supplies the
      corroboration. The exact thresholds and the maximum allowed strength are ultimately
      gene- and VCEP-specific (per ClinGen SVI recommendations), so a lab tunes the cut-offs;
      the viewer gives the two signals needed to make that call.
    </p>

    <h2>Benchmark — PM1</h2>
    <p>
      On an independent ACMG benchmark (support = 1,270 variants), the AlphaMissense-island
      PM1 is dramatically <b>more precise</b> than ClinVar-window approaches — ~3.6× the
      precision of InterVar and ~2× Franklin — and achieves the <b>best overall F1</b>.
    </p>
    <div class="bench">
      <table class="bench-table">
        <thead>
          <tr><th>Tool</th><th>Precision</th><th>Recall</th><th>F1</th></tr>
        </thead>
        <tbody>
          <tr class="win"><td>SeqOne <span class="tag">AM islands</span></td><td>0.70</td><td>0.50</td><td>0.59</td></tr>
          <tr><td>Franklin</td><td>0.33</td><td>0.86</td><td>0.47</td></tr>
          <tr><td>InterVar</td><td>0.19</td><td>0.84</td><td>0.31</td></tr>
        </tbody>
      </table>
      <figure class="bench-fig">
        <img src="/assets/benchmark_pm1.png" alt="PM1 precision, recall and F1 for SeqOne, InterVar and Franklin" />
        <figcaption>PM1 — mutational hotspot / functional domain. Precision / Recall / F1.</figcaption>
      </figure>
    </div>
    <p class="src">
      The high recall of window-based tools comes at a steep precision cost (they flag broad
      regions as hot-spots); the island approach trades a little recall for far fewer false
      positives — the right trade-off for an ACMG <i>moderate</i> criterion.
    </p>
    <p class="bench-cta"><a href="?page=benchmark" data-benchmark>Explore every PM1 call (TP / FP / FN) →</a></p>

    <p class="src">Cancer-hotspot track: cancerhotspots.org (Chang et al., 2018).</p>

    <p class="back"><a href="?page=methods" data-home>← Back to search</a></p>
  </article>`;
}
