/** Feedback page: an in-page form that POSTs to Formspree (no backend). */

const FORMSPREE = "https://formspree.io/f/xojgvlpe";

export function renderFeedback(root: HTMLElement): void {
  root.innerHTML = `
    <div class="page page-narrow">
      <p class="eyebrow">We'd love your input</p>
      <h1>Send <em>feedback</em></h1>
      <p class="page-lead">Bug, missing gene, a wrong-looking island, or an idea for the tool —
        tell us. A few words is plenty; an email is optional if you'd like a reply.</p>

      <form id="fb-form" class="fbform" novalidate>
        <label class="fb-field">
          <span class="fb-label">Topic</span>
          <select name="topic" class="fb-input">
            <option>General feedback</option>
            <option>Bug / something broken</option>
            <option>Feature request</option>
            <option>Data issue (island / ClinVar / gene)</option>
            <option>Scientific question</option>
          </select>
        </label>
        <label class="fb-field">
          <span class="fb-label">Your email <span class="fb-opt">(optional)</span></span>
          <input type="email" name="email" class="fb-input" placeholder="you@lab.org"
            autocomplete="email" spellcheck="false" />
        </label>
        <label class="fb-field">
          <span class="fb-label">Message</span>
          <textarea name="message" class="fb-input fb-textarea" rows="6" required
            placeholder="What happened, what you expected, or what you'd like to see…"></textarea>
        </label>
        <input type="text" name="_gotcha" class="fb-hp" tabindex="-1" autocomplete="off" aria-hidden="true" />
        <input type="hidden" name="page" />
        <div class="fb-actions">
          <button type="submit" class="fb-submit">Send feedback</button>
          <span id="fb-status" class="fb-status" role="status" aria-live="polite"></span>
        </div>
      </form>
      <p class="page-note">Submissions go to the Archipelago team via Formspree. We only read what
        you type here — no page or usage data is attached.</p>
    </div>`;

  const form = root.querySelector<HTMLFormElement>("#fb-form")!;
  const status = root.querySelector<HTMLElement>("#fb-status")!;
  const submit = form.querySelector<HTMLButtonElement>(".fb-submit")!;
  (form.querySelector('input[name="page"]') as HTMLInputElement).value =
    document.referrer || location.href;

  form.addEventListener("submit", async (e) => {
    e.preventDefault();
    const msg = (form.querySelector('[name="message"]') as HTMLTextAreaElement).value.trim();
    if (!msg) { setStatus("err", "Please add a message before sending."); return; }
    submit.disabled = true;
    setStatus("pending", "Sending…");
    try {
      const res = await fetch(FORMSPREE, {
        method: "POST",
        headers: { Accept: "application/json" },
        body: new FormData(form),
      });
      if (res.ok) {
        form.reset();
        setStatus("ok", "Thanks — your feedback is on its way. 🌱");
      } else {
        const data = await res.json().catch(() => null);
        const detail = data?.errors?.map((x: { message: string }) => x.message).join(", ");
        setStatus("err", detail || "Something went wrong. Please try again.");
        submit.disabled = false;
      }
    } catch {
      setStatus("err", "Network error — please check your connection and retry.");
      submit.disabled = false;
    }
  });

  function setStatus(kind: "pending" | "ok" | "err", text: string): void {
    status.className = `fb-status fb-${kind}`;
    status.textContent = text;
  }
}
