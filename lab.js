(() => {
  "use strict";

  const MU_MIN = 0.0;
  const MU_MAX = 0.36;
  const MU_DEFAULT = 0.048;

  const A1 = 1.40, B1 = 1.56, C1 = 1.40, D1 = -6.56;
  const A2 = 1.33, B2 = 1.49, C2 = 1.34, D2 = -6.38;
  const THETA = 0.22;
  const CT = Math.cos(THETA);
  const ST = Math.sin(THETA);

  const X_MIN = -8.35;
  const X_MAX = 8.35;
  const Y_MIN = -3.65;
  const Y_MAX = 3.65;

  const reducedMotionQuery = window.matchMedia("(prefers-reduced-motion: reduce)");

  function buildColorLUT() {
    const lut = new Uint32Array(256 * 256);
    const probe = new Uint8Array(new Uint32Array([0x11223344]).buffer);
    const littleEndian = probe[0] === 0x44;

    const obsR = 12.0, obsG = 12.0, obsB = 11.0;
    const umbR = 94.0, umbG = 60.0, umbB = 44.0;
    const copR = 223.0, copG = 159.0, copB = 100.0;
    const bonR = 236.0, bonG = 232.0, bonB = 225.0;

    for (let j = 0; j < 256; j++) {
      const v2 = j / 255.0;
      for (let i = 0; i < 256; i++) {
        const v1 = i / 255.0;
        const tot = Math.min(1.0, Math.max(0.0, 0.62 * v1 + 0.56 * v2));
        const sync = Math.min(1.0, Math.sqrt(v1 * v2) * 1.15);

        const sUmber = Math.pow(Math.min(1.0, tot * 1.48), 1.18);
        const sCopper = Math.pow(Math.min(1.0, Math.max(0.0, (v1 * 0.86 + 0.34 * v2 - 0.19) / 0.50)), 1.30);
        const maxPart = Math.max(0.0, (Math.max(v1, 0.85 * v2) - 0.55) / 0.40) * 0.80;
        const syncPart = Math.max(0.0, sync - 0.42) * 0.56;
        const sBone = Math.pow(Math.min(1.0, maxPart + syncPart), 1.50);

        let r = obsR * (1.0 - sUmber) + umbR * sUmber;
        let g = obsG * (1.0 - sUmber) + umbG * sUmber;
        let b = obsB * (1.0 - sUmber) + umbB * sUmber;

        r = r * (1.0 - sCopper) + copR * sCopper;
        g = g * (1.0 - sCopper) + copG * sCopper;
        b = b * (1.0 - sCopper) + copB * sCopper;

        r = r * (1.0 - sBone) + bonR * sBone;
        g = g * (1.0 - sBone) + bonG * sBone;
        b = b * (1.0 - sBone) + bonB * sBone;

        const ri = Math.min(255, Math.max(0, r + 0.5 | 0));
        const gi = Math.min(255, Math.max(0, g + 0.5 | 0));
        const bi = Math.min(255, Math.max(0, b + 0.5 | 0));

        lut[(j << 8) | i] = littleEndian
          ? (0xff000000 | (bi << 16) | (gi << 8) | ri) >>> 0
          : ((ri << 24) | (gi << 16) | (bi << 8) | 0xff) >>> 0;
      }
    }
    return lut;
  }

  function estimateDynamics(mu) {
    let x1 = 0.42, y1 = -0.31, x2 = 0.47, y2 = -0.26;
    for (let i = 0; i < 64; i++) {
      const nx1 = D1 * Math.sin(A1 * x1) - Math.sin(B1 * y1) + mu * Math.sin(y2 - y1);
      const ny1 = C1 * Math.cos(A1 * x1) + Math.cos(B1 * y1) + mu * Math.cos(x2 - x1);
      const nx2 = D2 * Math.sin(A2 * x2) - Math.sin(B2 * y2) + mu * Math.sin(y1 - y2);
      const ny2 = C2 * Math.cos(A2 * x2) + Math.cos(B2 * y2) + mu * Math.cos(x1 - x2);
      x1 = nx1; y1 = ny1; x2 = nx2; y2 = ny2;
    }

    const eps0 = 1e-6;
    let qx1 = x1 + eps0, qy1 = y1, qx2 = x2, qy2 = y2;
    let lyapSum = 0.0;
    let syncSum = 0.0;
    const steps = 140;

    for (let i = 0; i < steps; i++) {
      const nx1 = D1 * Math.sin(A1 * x1) - Math.sin(B1 * y1) + mu * Math.sin(y2 - y1);
      const ny1 = C1 * Math.cos(A1 * x1) + Math.cos(B1 * y1) + mu * Math.cos(x2 - x1);
      const nx2 = D2 * Math.sin(A2 * x2) - Math.sin(B2 * y2) + mu * Math.sin(y1 - y2);
      const ny2 = C2 * Math.cos(A2 * x2) + Math.cos(B2 * y2) + mu * Math.cos(x1 - x2);

      const pqx1 = D1 * Math.sin(A1 * qx1) - Math.sin(B1 * qy1) + mu * Math.sin(qy2 - qy1);
      const pqy1 = C1 * Math.cos(A1 * qx1) + Math.cos(B1 * qy1) + mu * Math.cos(qx2 - qx1);
      const pqx2 = D2 * Math.sin(A2 * qx2) - Math.sin(B2 * qy2) + mu * Math.sin(qy1 - qy2);
      const pqy2 = C2 * Math.cos(A2 * qx2) + Math.cos(B2 * qy2) + mu * Math.cos(qx1 - qx2);

      x1 = nx1; y1 = ny1; x2 = nx2; y2 = ny2;
      const dx1 = pqx1 - x1, dy1 = pqy1 - y1, dx2 = pqx2 - x2, dy2 = pqy2 - y2;
      const dist = Math.sqrt(dx1 * dx1 + dy1 * dy1 + dx2 * dx2 + dy2 * dy2) + 1e-15;
      lyapSum += Math.log(dist / eps0);

      const inv = eps0 / dist;
      qx1 = x1 + dx1 * inv;
      qy1 = y1 + dy1 * inv;
      qx2 = x2 + dx2 * inv;
      qy2 = y2 + dy2 * inv;

      const phaseSep = Math.hypot(x1 - x2, y1 - y2);
      syncSum += Math.exp(-0.35 * phaseSep);
    }

    return {
      lyap: lyapSum / steps,
      sync: syncSum / steps
    };
  }

  function smooth3Tap(src, tmp, out, w, h) {
    for (let y = 0; y < h; y++) {
      const row = y * w;
      tmp[row] = 0.78 * src[row] + 0.22 * src[row + 1];
      const end = row + w - 1;
      for (let i = row + 1; i < end; i++) {
        tmp[i] = 0.20 * src[i - 1] + 0.60 * src[i] + 0.20 * src[i + 1];
      }
      tmp[end] = 0.22 * src[end - 1] + 0.78 * src[end];
    }
    let peak = 1e-3;
    for (let y = 1; y < h - 1; y++) {
      const row = y * w;
      for (let x = 0; x < w; x++) {
        const i = row + x;
        const v = 0.68 * (0.20 * tmp[i - w] + 0.60 * tmp[i] + 0.20 * tmp[i + w]) + 0.32 * src[i];
        out[i] = v;
        if (v > peak) peak = v;
      }
    }
    return peak;
  }

  function initAttractorStage() {
    const canvas = document.getElementById("attractor-stage");
    if (!canvas) return;
    const ctx = canvas.getContext("2d", { alpha: false });
    if (!ctx) return;

    const tickEl = document.getElementById("stage-tick");
    const muEl = document.getElementById("tel-mu");
    const lyapEl = document.getElementById("tel-lyap");
    const syncEl = document.getElementById("tel-sync");
    const hintEl = document.getElementById("stage-hint");

    const W = canvas.width;
    const H = canvas.height;
    const N_PIXELS = W * H;
    const N_ORBITS = 3072;

    const colorLUT = buildColorLUT();
    const imgData = ctx.createImageData(W, H);
    const pixels32 = new Uint32Array(imgData.data.buffer);

    const d1 = new Float32Array(N_PIXELS);
    const d2 = new Float32Array(N_PIXELS);
    const tmpBuf = new Float32Array(N_PIXELS);
    const s1 = new Float32Array(N_PIXELS);
    const s2 = new Float32Array(N_PIXELS);

    const ox1 = new Float32Array(N_ORBITS);
    const oy1 = new Float32Array(N_ORBITS);
    const ox2 = new Float32Array(N_ORBITS);
    const oy2 = new Float32Array(N_ORBITS);

    let mu = MU_DEFAULT;
    let frameBudget = 0;
    let rafId = 0;
    let isDragging = false;

    function seedOrbits(currentMu) {
      let state = 0x10804208 >>> 0;
      const rand = () => {
        state ^= state << 13;
        state ^= state >>> 17;
        state ^= state << 5;
        return ((state >>> 0) / 4294967296.0);
      };

      for (let i = 0; i < N_ORBITS; i++) {
        const rx = (rand() * 3.0) - 1.5;
        const ry = (rand() * 3.0) - 1.5;
        ox1[i] = rx;
        oy1[i] = ry;
        ox2[i] = rx + (rand() - 0.5) * 0.14;
        oy2[i] = ry + (rand() - 0.5) * 0.14;
      }

      for (let s = 0; s < 32; s++) {
        for (let i = 0; i < N_ORBITS; i++) {
          const x1 = ox1[i], y1 = oy1[i], x2 = ox2[i], y2 = oy2[i];
          ox1[i] = D1 * Math.sin(A1 * x1) - Math.sin(B1 * y1) + currentMu * Math.sin(y2 - y1);
          oy1[i] = C1 * Math.cos(A1 * x1) + Math.cos(B1 * y1) + currentMu * Math.cos(x2 - x1);
          ox2[i] = D2 * Math.sin(A2 * x2) - Math.sin(B2 * y2) + currentMu * Math.sin(y1 - y2);
          oy2[i] = C2 * Math.cos(A2 * x2) + Math.cos(B2 * y2) + currentMu * Math.cos(x1 - x2);
        }
      }
    }

    function splatBatch(steps, currentMu) {
      const scaleX = (W - 1) / (X_MAX - X_MIN);
      const scaleY = (H - 1) / (Y_MAX - Y_MIN);
      const wLimit = W - 1.001;
      const hLimit = H - 1.001;

      for (let s = 0; s < steps; s++) {
        for (let i = 0; i < N_ORBITS; i++) {
          const x1 = ox1[i], y1 = oy1[i], x2 = ox2[i], y2 = oy2[i];
          const nx1 = D1 * Math.sin(A1 * x1) - Math.sin(B1 * y1) + currentMu * Math.sin(y2 - y1);
          const ny1 = C1 * Math.cos(A1 * x1) + Math.cos(B1 * y1) + currentMu * Math.cos(x2 - x1);
          const nx2 = D2 * Math.sin(A2 * x2) - Math.sin(B2 * y2) + currentMu * Math.sin(y1 - y2);
          const ny2 = C2 * Math.cos(A2 * x2) + Math.cos(B2 * y2) + currentMu * Math.cos(x1 - x2);
          ox1[i] = nx1; oy1[i] = ny1; ox2[i] = nx2; oy2[i] = ny2;

          const u1 = nx1 * CT - ny1 * ST;
          const v1 = nx1 * ST + ny1 * CT;
          const fx1 = (u1 - X_MIN) * scaleX;
          const fy1 = (Y_MAX - v1) * scaleY;
          if (fx1 >= 1.0 && fx1 < wLimit && fy1 >= 1.0 && fy1 < hLimit) {
            const ix = fx1 | 0;
            const iy = fy1 | 0;
            const dx = fx1 - ix;
            const dy = fy1 - iy;
            const idx = iy * W + ix;
            const w0 = 1.0 - dx;
            const h0 = 1.0 - dy;
            d1[idx] += w0 * h0;
            d1[idx + 1] += dx * h0;
            d1[idx + W] += w0 * dy;
            d1[idx + W + 1] += dx * dy;
          }

          const u2 = nx2 * CT - ny2 * ST;
          const v2 = nx2 * ST + ny2 * CT;
          const fx2 = (u2 - X_MIN) * scaleX;
          const fy2 = (Y_MAX - v2) * scaleY;
          if (fx2 >= 1.0 && fx2 < wLimit && fy2 >= 1.0 && fy2 < hLimit) {
            const ix = fx2 | 0;
            const iy = fy2 | 0;
            const dx = fx2 - ix;
            const dy = fy2 - iy;
            const idx = iy * W + ix;
            const w0 = 1.0 - dx;
            const h0 = 1.0 - dy;
            d2[idx] += w0 * h0;
            d2[idx + 1] += dx * h0;
            d2[idx + W] += w0 * dy;
            d2[idx + W + 1] += dx * dy;
          }
        }
      }
    }

    function paintBuffer() {
      const p1 = smooth3Tap(d1, tmpBuf, s1, W, H);
      const p2 = smooth3Tap(d2, tmpBuf, s2, W, H);
      const k1 = 1.0 / (0.015 * p1 + 1e-5);
      const k2 = 1.0 / (0.015 * p2 + 1e-5);
      const norm1 = 255.0 / Math.log1p(0.46 * p1 * k1 + 1.0);
      const norm2 = 255.0 / Math.log1p(0.46 * p2 * k2 + 1.0);
      const bgPixel = colorLUT[0];

      for (let idx = 0; idx < N_PIXELS; idx++) {
        const a = s1[idx];
        const b = s2[idx];
        if (a <= 0.0005 && b <= 0.0005) {
          pixels32[idx] = bgPixel;
          continue;
        }
        const i1 = a > 0.0 ? Math.min(255, (Math.log1p(a * k1) * norm1) | 0) : 0;
        const i2 = b > 0.0 ? Math.min(255, (Math.log1p(b * k2) * norm2) | 0) : 0;
        pixels32[idx] = colorLUT[(i2 << 8) | i1];
      }
      ctx.putImageData(imgData, 0, 0);
    }

    function updateTelemetry(currentMu) {
      const pct = ((currentMu - MU_MIN) / (MU_MAX - MU_MIN)) * 100.0;
      if (tickEl) tickEl.style.left = `${pct.toFixed(2)}%`;
      if (muEl) muEl.textContent = currentMu.toFixed(3);
      canvas.setAttribute("aria-valuenow", currentMu.toFixed(3));

      const dyn = estimateDynamics(currentMu);
      if (lyapEl) {
        const sign = dyn.lyap >= 0 ? "+" : "";
        lyapEl.textContent = `${sign}${dyn.lyap.toFixed(3)}`;
      }
      if (syncEl) {
        syncEl.textContent = dyn.sync.toFixed(2);
      }
    }

    function resetBuffers(currentMu) {
      d1.fill(0.0);
      d2.fill(0.0);
      s1.fill(0.0);
      s2.fill(0.0);
      seedOrbits(currentMu);
    }

    function renderStaticEquilibrium() {
      if (rafId) {
        cancelAnimationFrame(rafId);
        rafId = 0;
      }
      resetBuffers(mu);
      splatBatch(320, mu);
      paintBuffer();
      updateTelemetry(mu);
    }

    function tickProgressive() {
      rafId = 0;
      if (frameBudget <= 0) return;
      splatBatch(isDragging ? 32 : 42, mu);
      paintBuffer();
      frameBudget -= 1;
      if (frameBudget > 0) {
        rafId = requestAnimationFrame(tickProgressive);
      }
    }

    function triggerRender(initial) {
      updateTelemetry(mu);
      if (reducedMotionQuery.matches) {
        if (hintEl) hintEl.textContent = "scrub stage or press arrow keys · static equilibrium";
        renderStaticEquilibrium();
        return;
      }
      resetBuffers(mu);
      splatBatch(initial ? 140 : 44, mu);
      paintBuffer();
      frameBudget = isDragging ? 4 : 8;
      if (!rafId) {
        rafId = requestAnimationFrame(tickProgressive);
      }
    }

    function setMuFromClientX(clientX) {
      const rect = canvas.getBoundingClientRect();
      if (rect.width <= 0) return;
      const ratio = Math.min(1.0, Math.max(0.0, (clientX - rect.left) / rect.width));
      const nextMu = MU_MIN + ratio * (MU_MAX - MU_MIN);
      if (Math.abs(nextMu - mu) < 0.0004) return;
      mu = nextMu;
      triggerRender(false);
    }

    canvas.addEventListener("pointerdown", (ev) => {
      if (ev.button !== undefined && ev.button !== 0) return;
      isDragging = true;
      try {
        canvas.setPointerCapture(ev.pointerId);
      } catch (_) {}
      setMuFromClientX(ev.clientX);
    });

    canvas.addEventListener("pointermove", (ev) => {
      if (!isDragging) return;
      setMuFromClientX(ev.clientX);
    });

    const endDrag = (ev) => {
      if (!isDragging) return;
      isDragging = false;
      try {
        canvas.releasePointerCapture(ev.pointerId);
      } catch (_) {}
      if (!reducedMotionQuery.matches) {
        frameBudget = 8;
        if (!rafId) rafId = requestAnimationFrame(tickProgressive);
      }
    };

    canvas.addEventListener("pointerup", endDrag);
    canvas.addEventListener("pointercancel", endDrag);

    canvas.addEventListener("keydown", (ev) => {
      let step = 0.0;
      if (ev.key === "ArrowRight" || ev.key === "ArrowUp") step = 0.012;
      else if (ev.key === "ArrowLeft" || ev.key === "ArrowDown") step = -0.012;
      else if (ev.key === "Home") {
        mu = MU_MIN;
        ev.preventDefault();
        triggerRender(false);
        return;
      } else if (ev.key === "End") {
        mu = MU_MAX;
        ev.preventDefault();
        triggerRender(false);
        return;
      } else {
        return;
      }
      ev.preventDefault();
      mu = Math.min(MU_MAX, Math.max(MU_MIN, mu + step));
      triggerRender(false);
    });

    reducedMotionQuery.addEventListener("change", () => {
      if (hintEl) {
        hintEl.textContent = reducedMotionQuery.matches
          ? "scrub stage or press arrow keys · static equilibrium"
          : "scrub stage horizontally to perturb coupling";
      }
      triggerRender(true);
    });

    triggerRender(true);
  }

  function initContextScrubber() {
    const PROFILES = {
      c11: {
        catPerTurn: 9129.25,
        slicePerTurn: 346.0,
        savingPct: "-96.2%",
        latencyMs: "1.7 ms",
        symbol: "cmd_sandbox_impl",
        armor: "Landlock 0% CPU"
      },
      python: {
        catPerTurn: 7453.0,
        slicePerTurn: 484.0,
        savingPct: "-93.5%",
        latencyMs: "1.3 ms",
        symbol: "validate_pinned_image_ref",
        armor: "Pickle VM static scan"
      },
      rust: {
        catPerTurn: 9840.0,
        slicePerTurn: 442.0,
        savingPct: "-95.5%",
        latencyMs: "1.5 ms",
        symbol: "enforce_landlock_abi",
        armor: "ABI v1..v6 sandbox"
      }
    };

    const tabs = Array.from(document.querySelectorAll('.choice-word[role="tab"]'));
    const scrubber = document.getElementById("turn-scrubber");
    const fillEl = document.getElementById("turn-fill");
    const needleEl = document.getElementById("turn-needle");
    const turnReadout = document.getElementById("turn-readout");
    const savingEl = document.getElementById("scrubber-saving");
    const captionEl = document.getElementById("scrubber-caption");
    const valCat = document.getElementById("val-cat");
    const lblSlice = document.getElementById("lbl-slice");
    const valSlice = document.getElementById("val-slice");

    if (!scrubber) return;

    let activeLang = "c11";
    let activeTurn = 4;
    let dragging = false;

    const fmtInt = (n) => Math.round(n).toLocaleString("en-US");

    function renderState() {
      const p = PROFILES[activeLang] || PROFILES.c11;
      const catTok = Math.round(p.catPerTurn * activeTurn);
      const sliceTok = Math.round(p.slicePerTurn * activeTurn);
      const savedTok = Math.max(0, catTok - sliceTok);
      const catWinPct = ((catTok / 200000.0) * 100.0).toFixed(1);
      const sliceWinPct = ((sliceTok / 200000.0) * 100.0).toFixed(1);
      const turnPad = String(activeTurn).padStart(2, "0");
      const ratioPct = ((activeTurn - 1) / 11.0) * 100.0;

      if (fillEl) fillEl.style.width = `${ratioPct.toFixed(2)}%`;
      if (needleEl) needleEl.style.left = `${ratioPct.toFixed(2)}%`;
      if (turnReadout) turnReadout.textContent = `${turnPad} / 12 turns`;
      scrubber.setAttribute("aria-valuenow", String(activeTurn));

      if (savingEl) savingEl.textContent = p.savingPct;
      if (captionEl) {
        captionEl.textContent = `${fmtInt(savedTok)} tokens saved at turn ${turnPad} · ${p.latencyMs} mmap · ${p.armor}`;
      }
      if (valCat) valCat.textContent = `${fmtInt(catTok)} tok (${catWinPct}% of 200k window)`;
      if (lblSlice) lblSlice.textContent = `hslice mmap AST slice (${p.symbol})`;
      if (valSlice) valSlice.textContent = `${fmtInt(sliceTok)} tok (${sliceWinPct}% of 200k window)`;
    }

    function setTurnFromClientX(clientX) {
      const rect = scrubber.getBoundingClientRect();
      if (rect.width <= 0) return;
      const ratio = Math.min(1.0, Math.max(0.0, (clientX - rect.left) / rect.width));
      const nextTurn = Math.min(12, Math.max(1, Math.round(1 + ratio * 11)));
      if (nextTurn !== activeTurn) {
        activeTurn = nextTurn;
        renderState();
      }
    }

    tabs.forEach((tab, idx) => {
      tab.addEventListener("click", () => {
        const lang = tab.getAttribute("data-lang");
        if (!lang || !PROFILES[lang]) return;
        activeLang = lang;
        tabs.forEach((t) => t.setAttribute("aria-selected", t === tab ? "true" : "false"));
        renderState();
      });

      tab.addEventListener("keydown", (ev) => {
        if (ev.key !== "ArrowRight" && ev.key !== "ArrowLeft") return;
        ev.preventDefault();
        const nextIdx = (idx + (ev.key === "ArrowRight" ? 1 : tabs.length - 1)) % tabs.length;
        tabs[nextIdx].focus();
        tabs[nextIdx].click();
      });
    });

    scrubber.addEventListener("pointerdown", (ev) => {
      if (ev.button !== undefined && ev.button !== 0) return;
      dragging = true;
      try {
        scrubber.setPointerCapture(ev.pointerId);
      } catch (_) {}
      setTurnFromClientX(ev.clientX);
    });

    scrubber.addEventListener("pointermove", (ev) => {
      if (!dragging) return;
      setTurnFromClientX(ev.clientX);
    });

    const stopScrub = (ev) => {
      if (!dragging) return;
      dragging = false;
      try {
        scrubber.releasePointerCapture(ev.pointerId);
      } catch (_) {}
    };

    scrubber.addEventListener("pointerup", stopScrub);
    scrubber.addEventListener("pointercancel", stopScrub);

    scrubber.addEventListener("keydown", (ev) => {
      if (ev.key === "ArrowRight" || ev.key === "ArrowUp") {
        ev.preventDefault();
        activeTurn = Math.min(12, activeTurn + 1);
        renderState();
      } else if (ev.key === "ArrowLeft" || ev.key === "ArrowDown") {
        ev.preventDefault();
        activeTurn = Math.max(1, activeTurn - 1);
        renderState();
      } else if (ev.key === "Home") {
        ev.preventDefault();
        activeTurn = 1;
        renderState();
      } else if (ev.key === "End") {
        ev.preventDefault();
        activeTurn = 12;
        renderState();
      }
    });

    renderState();
  }

  function initCopyLines() {
    const lines = document.querySelectorAll(".cmd-line[data-copy]");
    lines.forEach((btn) => {
      let resetTimer = 0;
      btn.addEventListener("click", async () => {
        const text = btn.getAttribute("data-copy") || "";
        if (!text) return;
        const actionEl = btn.querySelector(".cmd-action");
        let ok = false;
        if (navigator.clipboard && navigator.clipboard.writeText) {
          try {
            await navigator.clipboard.writeText(text);
            ok = true;
          } catch (_) {}
        }
        if (!ok) {
          const area = document.createElement("textarea");
          area.value = text;
          area.setAttribute("readonly", "");
          area.style.position = "fixed";
          area.style.opacity = "0";
          document.body.appendChild(area);
          area.select();
          try {
            ok = document.execCommand("copy");
          } catch (_) {}
          document.body.removeChild(area);
        }
        if (actionEl && ok) {
          btn.setAttribute("data-copied", "true");
          actionEl.textContent = "copied";
          if (resetTimer) clearTimeout(resetTimer);
          resetTimer = window.setTimeout(() => {
            btn.removeAttribute("data-copied");
            actionEl.textContent = "copy";
          }, 1200);
        }
      });
    });
  }

  initAttractorStage();
  initContextScrubber();
  initCopyLines();
})();
