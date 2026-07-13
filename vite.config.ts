import { defineConfig } from "vite";

// Static SPA. Precomputed gene bundles live in public/data and are copied
// verbatim into the build output; they are served as static files on Vercel.
export default defineConfig({
  build: {
    target: "es2021",
    chunkSizeWarningLimit: 1200,
  },
});
