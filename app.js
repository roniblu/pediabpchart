/* =========================================================================
   PediBP — Pediatric Blood Pressure Percentile Calculator
   2017 American Academy of Pediatrics Clinical Practice Guideline
   ========================================================================= */

/* ----- 1. CDC 2000 Stature-for-age LMS data ----------------------------- */
/* Loaded from lms-data.js (sets LMS_BOYS and LMS_GIRLS) */

/* ----- 2. Math helpers --------------------------------------------------- */

/** Standard normal CDF via Abramowitz & Stegun 7.1.26 */
function normalCDF(z) {
  const sign = z < 0 ? -1 : 1;
  const x = Math.abs(z) / Math.SQRT2;
  // erf approximation
  const t = 1 / (1 + 0.3275911 * x);
  const y =
    1 -
    (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t +
      0.254829592) *
      t *
      Math.exp(-x * x);
  return 0.5 * (1 + sign * y);
}

/** Linearly interpolate LMS values for a fractional age in months. */
function interpolateLMS(ageMonths, table) {
  if (ageMonths <= table[0][0]) {
    const r = table[0];
    return { L: r[1], M: r[2], S: r[3] };
  }
  if (ageMonths >= table[table.length - 1][0]) {
    const r = table[table.length - 1];
    return { L: r[1], M: r[2], S: r[3] };
  }
  for (let i = 0; i < table.length - 1; i++) {
    const a = table[i],
      b = table[i + 1];
    if (ageMonths >= a[0] && ageMonths <= b[0]) {
      const f = (ageMonths - a[0]) / (b[0] - a[0]);
      return {
        L: a[1] + (b[1] - a[1]) * f,
        M: a[2] + (b[2] - a[2]) * f,
        S: a[3] + (b[3] - a[3]) * f,
      };
    }
  }
  return null;
}

/** Height Z from height in cm using LMS Box-Cox. */
function heightZFromCm(heightCm, ageMonths, sex) {
  const table = sex === 'M' ? LMS_BOYS : LMS_GIRLS;
  const lms = interpolateLMS(ageMonths, table);
  if (!lms) return null;
  const { L, M, S } = lms;
  if (Math.abs(L) < 1e-8) return Math.log(heightCm / M) / S;
  return (Math.pow(heightCm / M, L) - 1) / (L * S);
}

/** Height in cm from height percentile (1-99). */
function heightCmFromPercentile(percentile, ageMonths, sex) {
  const table = sex === 'M' ? LMS_BOYS : LMS_GIRLS;
  const lms = interpolateLMS(ageMonths, table);
  if (!lms) return null;
  const { L, M, S } = lms;
  // Inverse normal via Beasley-Springer-Moro
  const z = invNormal(percentile / 100);
  if (Math.abs(L) < 1e-8) return M * Math.exp(S * z);
  return M * Math.pow(1 + L * S * z, 1 / L);
}

/** Inverse normal CDF (Acklam's algorithm — accurate to ~1e-9). */
function invNormal(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  const a = [-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2, 1.38357751867269e2, -3.066479806614716e1, 2.506628277459239];
  const b = [-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2, 6.680131188771972e1, -1.328068155288572e1];
  const c = [-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
  const d = [7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996, 3.754408661907416];
  const pLow = 0.02425;
  const pHigh = 1 - pLow;
  let q, r;
  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
      ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else if (p <= pHigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
      (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
      ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
}

/* ----- 3. AAP 2017 / 4th Report BP regression --------------------------- */
/* Polynomial regression of conditional mean BP on age and height-Z.
   Coefficients from AAP 2017 appendix (carried forward from the 2004 4th Report).
   Z-score = (BP - mean) / sigma, percentile = Phi(Z) * 100. */

const BP_COEF = {
  SBP_M: {
    intercept: 102.19768,
    age: [1.82416, 0.12776, 0.00249, -0.00135],
    htZ: [2.73157, -0.19618, -0.04659, 0.00947],
    sigma: 10.7128,
  },
  SBP_F: {
    intercept: 102.01027,
    age: [1.94397, 0.00598, -0.00789, -0.00059],
    htZ: [2.03526, 0.02534, -0.01884, 0.00121],
    sigma: 10.4855,
  },
  DBP_M: {
    intercept: 61.01217,
    age: [0.68314, -0.09835, 0.01711, 0.00045],
    htZ: [1.46993, -0.07849, -0.03144, 0.00967],
    sigma: 11.6032,
  },
  DBP_F: {
    intercept: 60.5051,
    age: [1.01301, 0.01157, 0.00424, -0.00137],
    htZ: [1.16641, 0.12795, -0.03869, -0.00079],
    sigma: 10.9573,
  },
};

function bpMean(coef, ageMonths, htZ) {
  const a = ageMonths / 12 - 10;
  const a2 = a * a, a3 = a2 * a, a4 = a3 * a;
  const h = htZ;
  const h2 = h * h, h3 = h2 * h, h4 = h3 * h;
  return (
    coef.intercept +
    coef.age[0] * a + coef.age[1] * a2 + coef.age[2] * a3 + coef.age[3] * a4 +
    coef.htZ[0] * h + coef.htZ[1] * h2 + coef.htZ[2] * h3 + coef.htZ[3] * h4
  );
}

function bpZ(measuredBP, coef, ageMonths, htZ) {
  const mean = bpMean(coef, ageMonths, htZ);
  return { z: (measuredBP - mean) / coef.sigma, mean, sigma: coef.sigma };
}

function bpAtPercentile(coef, ageMonths, htZ, percentile) {
  const mean = bpMean(coef, ageMonths, htZ);
  const z = invNormal(percentile / 100);
  return mean + z * coef.sigma;
}

/* ----- 4. AAP 2017 classification --------------------------------------- */
/* For 1 ≤ age < 13: lower of percentile-based and absolute-mmHg cutoffs.
   For age ≥ 13: absolute adult thresholds. */

const CATEGORIES = ['normal', 'elevated', 'stage1', 'stage2'];
const CATEGORY_LABELS = {
  normal: 'Normal',
  elevated: 'Elevated',
  stage1: 'Stage 1 hypertension',
  stage2: 'Stage 2 hypertension',
};

function classifyChild(sbp, dbp, sbp_p90, sbp_p95, dbp_p90, dbp_p95) {
  // Returns the higher of SBP and DBP categories.
  function cat(value, p90, p95) {
    // percentile-based thresholds
    const p95plus12 = p95 + 12;
    let pctCat;
    if (value < p90) pctCat = 'normal';
    else if (value < p95) pctCat = 'elevated';
    else if (value < p95plus12) pctCat = 'stage1';
    else pctCat = 'stage2';
    return pctCat;
  }
  function catAbs(sbp, dbp) {
    // For children the absolute mmHg thresholds: 120/80 elevated, 130/80 stage1, 140/90 stage2
    if (sbp >= 140 || dbp >= 90) return 'stage2';
    if (sbp >= 130 || dbp >= 80) return 'stage1';
    if (sbp >= 120) return 'elevated';
    return 'normal';
  }
  const sbpCat = cat(sbp, sbp_p90, sbp_p95);
  const dbpCat = cat(dbp, dbp_p90, dbp_p95);
  const absCat = catAbs(sbp, dbp);
  // Lowest category between percentile-based-individual and absolute-cap
  // AAP rule says use the LOWER (stricter) numeric cutoff, which means whichever
  // method gives the HIGHER category for that BP value.
  // Wait: "whichever is lower" refers to the cutoff threshold, not the category.
  // E.g., if 95th percentile is 132 mmHg but 130 mmHg absolute is lower, the 130 mmHg
  // threshold defines stage 1 - which means a child with SBP 131 falls into stage 1
  // even though they're below the 95th percentile.
  // So: take the higher of the percentile-based and absolute-based categories.
  const sbpFinal = higherCat(sbpCat, catAbsSBP(sbp));
  const dbpFinal = higherCat(dbpCat, catAbsDBP(dbp));
  return {
    sbpCategory: sbpFinal,
    dbpCategory: dbpFinal,
    overall: higherCat(sbpFinal, dbpFinal),
  };
}

function catAbsSBP(sbp) {
  if (sbp >= 140) return 'stage2';
  if (sbp >= 130) return 'stage1';
  if (sbp >= 120) return 'elevated';
  return 'normal';
}
function catAbsDBP(dbp) {
  if (dbp >= 90) return 'stage2';
  if (dbp >= 80) return 'stage1';
  return 'normal';
}
function higherCat(a, b) {
  return CATEGORIES.indexOf(a) > CATEGORIES.indexOf(b) ? a : b;
}

function classifyAdolescent(sbp, dbp) {
  // ≥13 y: adult-style thresholds
  let cat = 'normal';
  if (sbp >= 140 || dbp >= 90) cat = 'stage2';
  else if (sbp >= 130 || dbp >= 80) cat = 'stage1';
  else if (sbp >= 120 && dbp < 80) cat = 'elevated';
  // SBP categorization
  let sbpCat = 'normal';
  if (sbp >= 140) sbpCat = 'stage2';
  else if (sbp >= 130) sbpCat = 'stage1';
  else if (sbp >= 120) sbpCat = 'elevated';
  let dbpCat = 'normal';
  if (dbp >= 90) dbpCat = 'stage2';
  else if (dbp >= 80) dbpCat = 'stage1';
  return { sbpCategory: sbpCat, dbpCategory: dbpCat, overall: cat };
}

/* ----- 5. End-to-end calculation ---------------------------------------- */

function calculate(input) {
  // input: { ageYears, sex, heightCm OR heightPct, sbp, dbp }
  const ageMonths = input.ageYears * 12;
  if (ageMonths < 12 || ageMonths > 215.99) {
    return { error: 'Age must be between 1.0 and 17.99 years.' };
  }

  // Derive height-Z and height in cm
  let heightCm, heightZ, heightPct;
  if (input.heightCm != null) {
    heightCm = input.heightCm;
    heightZ = heightZFromCm(heightCm, ageMonths, input.sex);
    if (heightZ == null) return { error: 'Could not compute height percentile.' };
    heightPct = normalCDF(heightZ) * 100;
  } else if (input.heightPct != null) {
    heightPct = input.heightPct;
    heightZ = invNormal(heightPct / 100);
    heightCm = heightCmFromPercentile(heightPct, ageMonths, input.sex);
  } else {
    return { error: 'Provide height in cm or as a percentile.' };
  }

  // Clamp height-Z to a reasonable range (formula was fit on the bulk of data)
  const clampedHtZ = Math.max(-3, Math.min(3, heightZ));

  // Compute BP percentiles
  const sbpKey = input.sex === 'M' ? 'SBP_M' : 'SBP_F';
  const dbpKey = input.sex === 'M' ? 'DBP_M' : 'DBP_F';
  const sbpRes = bpZ(input.sbp, BP_COEF[sbpKey], ageMonths, clampedHtZ);
  const dbpRes = bpZ(input.dbp, BP_COEF[dbpKey], ageMonths, clampedHtZ);
  const sbpPct = normalCDF(sbpRes.z) * 100;
  const dbpPct = normalCDF(dbpRes.z) * 100;
  const sbp_p50 = bpAtPercentile(BP_COEF[sbpKey], ageMonths, clampedHtZ, 50);
  const sbp_p90 = bpAtPercentile(BP_COEF[sbpKey], ageMonths, clampedHtZ, 90);
  const sbp_p95 = bpAtPercentile(BP_COEF[sbpKey], ageMonths, clampedHtZ, 95);
  const dbp_p50 = bpAtPercentile(BP_COEF[dbpKey], ageMonths, clampedHtZ, 50);
  const dbp_p90 = bpAtPercentile(BP_COEF[dbpKey], ageMonths, clampedHtZ, 90);
  const dbp_p95 = bpAtPercentile(BP_COEF[dbpKey], ageMonths, clampedHtZ, 95);

  // Classify
  let classification;
  if (input.ageYears >= 13) {
    classification = classifyAdolescent(input.sbp, input.dbp);
  } else {
    classification = classifyChild(input.sbp, input.dbp, sbp_p90, sbp_p95, dbp_p90, dbp_p95);
  }

  return {
    ageMonths,
    sex: input.sex,
    heightCm,
    heightZ,
    heightPct,
    sbp: input.sbp,
    sbpZ: sbpRes.z,
    sbpPercentile: sbpPct,
    sbpMean: sbpRes.mean,
    sbpSigma: sbpRes.sigma,
    sbp_p50, sbp_p90, sbp_p95, sbp_p95plus12: sbp_p95 + 12,
    dbp: input.dbp,
    dbpZ: dbpRes.z,
    dbpPercentile: dbpPct,
    dbpMean: dbpRes.mean,
    dbpSigma: dbpRes.sigma,
    dbp_p50, dbp_p90, dbp_p95, dbp_p95plus12: dbp_p95 + 12,
    classification,
  };
}

/* ----- 6. UI rendering -------------------------------------------------- */

const RISK_META = {
  normal: {
    color: 'var(--risk-normal)',
    soft: 'var(--risk-normal-soft)',
    strong: 'var(--risk-normal-strong)',
    icon: '<svg viewBox="0 0 24 24" width="22" height="22" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><polyline points="20 6 9 17 4 12"/></svg>',
    detail: 'BP is below the 90th percentile (or below 120/80 in adolescents).',
    action:
      'Repeat at the next routine well-child visit. Reinforce healthy lifestyle (DASH-type diet, ≥60 min daily activity, normal BMI, sleep, sodium &lt; 2.3 g/day).',
  },
  elevated: {
    color: 'var(--risk-elevated)',
    soft: 'var(--risk-elevated-soft)',
    strong: 'var(--risk-elevated-strong)',
    icon: '<svg viewBox="0 0 24 24" width="22" height="22" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="9"/><line x1="12" y1="8" x2="12" y2="13"/><circle cx="12" cy="16.5" r="0.6" fill="currentColor"/></svg>',
    detail:
      'BP is between the 90th and 95th percentile (or 120/&lt;80 to 129/&lt;80 in adolescents).',
    action:
      'Lifestyle counselling. Recheck BP in 6 months. If still elevated, repeat in another 6 months — three elevated readings overall warrant ABPM and workup.',
  },
  stage1: {
    color: 'var(--risk-stage1)',
    soft: 'var(--risk-stage1-soft)',
    strong: 'var(--risk-stage1-strong)',
    icon: '<svg viewBox="0 0 24 24" width="22" height="22" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><path d="M10.29 3.86 1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"/><line x1="12" y1="9" x2="12" y2="13"/><circle cx="12" cy="17" r="0.6" fill="currentColor"/></svg>',
    detail:
      'BP at or above the 95th percentile but &lt; 95th + 12 mmHg (or 130/80 to 139/89 in adolescents).',
    action:
      'Recheck BP and average over 1–2 weeks. If sustained at 3 visits, diagnose HTN — order ABPM, basic labs (CBC, BUN/Cr, lytes, UA, renal US), and lipid screening. Begin lifestyle therapy.',
  },
  stage2: {
    color: 'var(--risk-stage2)',
    soft: 'var(--risk-stage2-soft)',
    strong: 'var(--risk-stage2-strong)',
    icon: '<svg viewBox="0 0 24 24" width="22" height="22" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="9"/><line x1="15" y1="9" x2="9" y2="15"/><line x1="9" y1="9" x2="15" y2="15"/></svg>',
    detail:
      'BP at or above the 95th percentile + 12 mmHg (or ≥ 140/90 in adolescents).',
    action:
      'Confirm with repeat BP at the same visit. If still ≥ stage 2, refer within 1 week (or evaluate immediately if symptomatic). ABPM, labs, echocardiogram, and consideration of pharmacotherapy.',
  },
};

function fmtPct(p) {
  if (p < 0.1) return '<0.1';
  if (p > 99.9) return '>99.9';
  return p.toFixed(1);
}
function fmtZ(z) {
  return (z >= 0 ? '+' : '') + z.toFixed(2);
}

function renderResults(r, formInput) {
  const root = document.getElementById('results');
  if (r.error) {
    root.innerHTML = `<div class="results-empty"><h3>Cannot compute</h3><p>${r.error}</p></div>`;
    return;
  }
  const meta = RISK_META[r.classification.overall];

  // Position bar markers
  const sbpBarPos = Math.max(0, Math.min(100, r.sbpPercentile));
  const dbpBarPos = Math.max(0, Math.min(100, r.dbpPercentile));

  const heightPctText = fmtPct(r.heightPct);
  const heightCmText = r.heightCm.toFixed(1);

  const ageY = (r.ageMonths / 12).toFixed(1);
  const sexLabel = r.sex === 'M' ? 'Male' : 'Female';

  const sbpRiskMeta = RISK_META[r.classification.sbpCategory];
  const dbpRiskMeta = RISK_META[r.classification.dbpCategory];

  root.innerHTML = `
    <article class="result-card" style="--risk-color: ${meta.color}; --risk-soft: ${meta.soft}; --risk-strong: ${meta.strong};">
      <header class="verdict">
        <div class="verdict-icon" aria-hidden="true">${meta.icon}</div>
        <div class="verdict-body">
          <p class="verdict-label">2017 AAP classification</p>
          <h3 class="verdict-title">${CATEGORY_LABELS[r.classification.overall]}</h3>
          <p class="verdict-detail">${meta.detail}</p>
        </div>
      </header>

      <div class="metrics">
        <div class="metric" style="--metric-color: ${sbpRiskMeta.color};">
          <div class="metric-label"><span class="pip"></span>Systolic — ${CATEGORY_LABELS[r.classification.sbpCategory]}</div>
          <div class="metric-value">${r.sbp}<span class="unit">mmHg</span></div>
          <div class="metric-rows">
            <div class="metric-row"><span>Percentile</span><strong>${fmtPct(r.sbpPercentile)}</strong></div>
            <div class="metric-row"><span>Z-score</span><strong>${fmtZ(r.sbpZ)}</strong></div>
            <div class="metric-row"><span>50th / 90th / 95th</span><strong>${r.sbp_p50.toFixed(0)} / ${r.sbp_p90.toFixed(0)} / ${r.sbp_p95.toFixed(0)}</strong></div>
          </div>
          <div class="bar" aria-hidden="true">
            <div class="bar-marker" style="left: ${sbpBarPos}%;"></div>
          </div>
          <div class="bar-axis"><span>0</span><span>50</span><span>90</span><span>95</span><span>100</span></div>
        </div>

        <div class="metric" style="--metric-color: ${dbpRiskMeta.color};">
          <div class="metric-label"><span class="pip"></span>Diastolic — ${CATEGORY_LABELS[r.classification.dbpCategory]}</div>
          <div class="metric-value">${r.dbp}<span class="unit">mmHg</span></div>
          <div class="metric-rows">
            <div class="metric-row"><span>Percentile</span><strong>${fmtPct(r.dbpPercentile)}</strong></div>
            <div class="metric-row"><span>Z-score</span><strong>${fmtZ(r.dbpZ)}</strong></div>
            <div class="metric-row"><span>50th / 90th / 95th</span><strong>${r.dbp_p50.toFixed(0)} / ${r.dbp_p90.toFixed(0)} / ${r.dbp_p95.toFixed(0)}</strong></div>
          </div>
          <div class="bar" aria-hidden="true">
            <div class="bar-marker" style="left: ${dbpBarPos}%;"></div>
          </div>
          <div class="bar-axis"><span>0</span><span>50</span><span>90</span><span>95</span><span>100</span></div>
        </div>
      </div>

      <div class="aux">
        <dl>
          <dt>Patient</dt><dd>${ageY} y · ${sexLabel}</dd>
          <dt>Height</dt><dd>${heightCmText} cm · ${heightPctText}th %ile (Z ${fmtZ(r.heightZ)})</dd>
          <dt>SBP / DBP</dt><dd>${r.sbp} / ${r.dbp} mmHg</dd>
          <dt>SBP percentile / Z</dt><dd>${fmtPct(r.sbpPercentile)}th · Z ${fmtZ(r.sbpZ)}</dd>
          <dt>DBP percentile / Z</dt><dd>${fmtPct(r.dbpPercentile)}th · Z ${fmtZ(r.dbpZ)}</dd>
          <dt>SBP cutoffs (90 / 95 / 95+12)</dt><dd>${r.sbp_p90.toFixed(0)} / ${r.sbp_p95.toFixed(0)} / ${r.sbp_p95plus12.toFixed(0)} mmHg</dd>
          <dt>DBP cutoffs (90 / 95 / 95+12)</dt><dd>${r.dbp_p90.toFixed(0)} / ${r.dbp_p95.toFixed(0)} / ${r.dbp_p95plus12.toFixed(0)} mmHg</dd>
        </dl>
        <div class="action-tip">
          <svg viewBox="0 0 24 24" width="20" height="20" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M12 2 L2 7 L12 12 L22 7 Z"/><path d="M2 17 L12 22 L22 17"/><path d="M2 12 L12 17 L22 12"/></svg>
          <span><strong>Recommended action — </strong>${meta.action}</span>
        </div>
      </div>
    </article>
  `;
}

/* ----- 7. Form wiring --------------------------------------------------- */

document.addEventListener('DOMContentLoaded', () => {
  // Theme toggle
  (function () {
    const t = document.querySelector('[data-theme-toggle]');
    const r = document.documentElement;
    let d = matchMedia('(prefers-color-scheme:dark)').matches ? 'dark' : 'light';
    r.setAttribute('data-theme', d);
    function paint() {
      t.innerHTML =
        d === 'dark'
          ? '<svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.93 4.93l1.41 1.41M17.66 17.66l1.41 1.41M2 12h2M20 12h2M4.93 19.07l1.41-1.41M17.66 6.34l1.41-1.41"/></svg>'
          : '<svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/></svg>';
      t.setAttribute('aria-label', 'Switch to ' + (d === 'dark' ? 'light' : 'dark') + ' mode');
    }
    paint();
    t.addEventListener('click', () => {
      d = d === 'dark' ? 'light' : 'dark';
      r.setAttribute('data-theme', d);
      paint();
    });
  })();

  // Height-mode segmented control
  const cmField = document.getElementById('height-cm-field');
  const pctField = document.getElementById('height-pct-field');
  document.querySelectorAll('input[name="htmode"]').forEach((el) => {
    el.addEventListener('change', () => {
      const mode = document.querySelector('input[name="htmode"]:checked').value;
      if (mode === 'cm') {
        cmField.classList.remove('hidden');
        pctField.classList.add('hidden');
      } else {
        cmField.classList.add('hidden');
        pctField.classList.remove('hidden');
      }
    });
  });

  // Submit handler
  const form = document.getElementById('bp-form');
  const errorEl = document.getElementById('form-error');
  form.addEventListener('submit', (e) => {
    e.preventDefault();
    errorEl.hidden = true;
    errorEl.textContent = '';

    const ageYears = parseFloat(document.getElementById('age').value);
    const sex = document.querySelector('input[name="sex"]:checked').value;
    const sbp = parseFloat(document.getElementById('sbp').value);
    const dbp = parseFloat(document.getElementById('dbp').value);
    const mode = document.querySelector('input[name="htmode"]:checked').value;

    const errs = [];
    if (!isFinite(ageYears) || ageYears < 1 || ageYears >= 18)
      errs.push('Enter age between 1 and 17.99 years.');
    if (!isFinite(sbp) || sbp < 50 || sbp > 220)
      errs.push('Enter a plausible systolic BP (50–220 mmHg).');
    if (!isFinite(dbp) || dbp < 20 || dbp > 140)
      errs.push('Enter a plausible diastolic BP (20–140 mmHg).');
    if (sbp && dbp && dbp >= sbp) errs.push('Diastolic BP must be lower than systolic BP.');

    let input = { ageYears, sex, sbp, dbp };
    if (mode === 'cm') {
      const heightCm = parseFloat(document.getElementById('height').value);
      if (!isFinite(heightCm) || heightCm < 60 || heightCm > 210)
        errs.push('Enter a height between 60 and 210 cm.');
      input.heightCm = heightCm;
    } else {
      const heightPct = parseFloat(document.getElementById('height-pct').value);
      if (!isFinite(heightPct) || heightPct < 0.1 || heightPct > 99.9)
        errs.push('Enter a height percentile between 0.1 and 99.9.');
      input.heightPct = heightPct;
    }

    if (errs.length) {
      errorEl.textContent = errs.join(' ');
      errorEl.hidden = false;
      return;
    }

    const result = calculate(input);
    renderResults(result, input);
    // Smooth scroll on small screens so user sees the result
    if (window.innerWidth < 880) {
      document.getElementById('results').scrollIntoView({ behavior: 'smooth', block: 'start' });
    }
  });

  // Reset
  document.getElementById('reset-btn').addEventListener('click', () => {
    form.reset();
    errorEl.hidden = true;
    document.getElementById('results').innerHTML = `
      <div class="results-empty">
        <div class="empty-icon" aria-hidden="true">
          <svg width="44" height="44" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round">
            <path d="M3 12h4l3-9 4 18 3-9h4"/>
          </svg>
        </div>
        <h3>Awaiting input</h3>
        <p>Fill the patient fields, then press Calculate to see percentiles, Z-scores, and the AAP risk classification.</p>
      </div>`;
    // re-show cm field
    cmField.classList.remove('hidden');
    pctField.classList.add('hidden');
    document.getElementById('ht-cm').checked = true;
  });

  // Example loader
  document.getElementById('example-btn').addEventListener('click', () => {
    document.getElementById('age').value = '8';
    document.getElementById('sex-m').checked = true;
    document.getElementById('ht-cm').checked = true;
    cmField.classList.remove('hidden');
    pctField.classList.add('hidden');
    document.getElementById('height').value = '128';
    document.getElementById('sbp').value = '118';
    document.getElementById('dbp').value = '75';
    form.requestSubmit();
  });
});
