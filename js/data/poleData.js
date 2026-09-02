/**
 * poleData.js
 * ---------------------------------------------------------------------------
 * Single-source-of-truth accessor for the 3D-POLE-Modeler viewers.
 *
 * Loads data/pole_canonical_data.json once, caches it, and exposes a small
 * read API so every page (index.html, mutation.html, variant.html, ...) reads
 * domains / residues / variants from ONE place instead of re-declaring them.
 *
 * Design notes
 *  - ES module, zero dependencies. Pages must use <script type="module">.
 *  - Async: fetch() needs an HTTP origin. Works on GitHub Pages and with
 *    `npx serve .`; does NOT work from a file:// URL (the project already
 *    runs behind a server, so this is fine).
 *  - Every record in the JSON carries a `verification_status`. Anything that
 *    is not 'verified' / 'verified_as_disputed' is treated as PROVISIONAL so
 *    the renderer can style it distinctly (see isProvisional / buildSceneModel).
 *  - All residue positions are 1-based on UniProt Q07864 (2286 aa).
 *
 * Quick start
 *   import { loadPoleData, buildSceneModel } from './js/data/poleData.js';
 *   const data  = await loadPoleData();
 *   const model = buildSceneModel(data);   // normalized shape for the scene
 * ---------------------------------------------------------------------------
 */

/** Default location of the canonical JSON, resolved relative to THIS module. */
const DEFAULT_DATA_URL = new URL('../../data/pole_canonical_data.json', import.meta.url);

/** Statuses that count as fully verified. Everything else is provisional. */
const VERIFIED_STATUSES = new Set(['verified', 'verified_as_disputed']);

let _data = null;      // resolved data object once loaded
let _loading = null;   // in-flight promise (dedupes concurrent callers)

/* ===========================================================================
 * Loading
 * ======================================================================== */

/**
 * Load (and cache) the canonical data. Safe to call from every page; the
 * fetch happens at most once per session.
 * @param {string|URL} [url] override the data location
 * @returns {Promise<object>} the parsed canonical data
 */
export async function loadPoleData(url = DEFAULT_DATA_URL) {
  if (_data) return _data;
  if (_loading) return _loading;
  _loading = (async () => {
    const res = await fetch(url, { cache: 'no-cache' });
    if (!res.ok) {
      throw new Error(`[poleData] failed to load ${url}: ${res.status} ${res.statusText}`);
    }
    const json = await res.json();
    if (!json || !Array.isArray(json.domains) || !Array.isArray(json.variants)) {
      throw new Error('[poleData] loaded JSON is missing expected top-level arrays.');
    }
    _data = json;
    _loading = null;
    return json;
  })();
  return _loading;
}

/** Return the already-loaded data, or throw if loadPoleData() has not resolved. */
export function getLoadedData() {
  if (!_data) {
    throw new Error('[poleData] data not loaded yet — `await loadPoleData()` first.');
  }
  return _data;
}

/** Testing/util: clear the cache. */
function _resetCache() { _data = null; _loading = null; }

/* ===========================================================================
 * Provenance / verification status
 * ======================================================================== */

/** True if a record (domain/residue/variant/metal site) is NOT fully verified. */
export function isProvisional(rec) {
  return !VERIFIED_STATUSES.has(rec && rec.verification_status);
}

/** A short tag for tooltips: 'verified' or the specific provisional status. */
export function provenanceTag(rec) {
  return isProvisional(rec) ? (rec && rec.verification_status) || 'unverified' : 'verified';
}

/**
 * Resolve a record's `sources` keys into full reference objects.
 * @returns {Array<object>} each entry is the reference plus its `key`
 */
export function resolveSources(rec, data = getLoadedData()) {
  const refs = data.references || {};
  return ((rec && rec.sources) || []).map((k) => ({ key: k, ...(refs[k] || { _missing: true }) }));
}

/* ===========================================================================
 * Meta / references / notes
 * ======================================================================== */

export const getMeta = (data = getLoadedData()) => data.meta || {};
export const getReferences = (data = getLoadedData()) => data.references || {};
export const getReferenceById = (key, data = getLoadedData()) => (data.references || {})[key] || null;

/**
 * A compact provenance string for a UI footer, e.g.
 *   "POLE Q07864 · NM_006231.4 · GRCh38 · data schema 0.1.0-draft · generated 2026-05-29"
 */
export function provenanceFooter(data = getLoadedData()) {
  const m = getMeta(data);
  return [
    m.gene && `${m.gene} ${m.uniprot || ''}`.trim(),
    m.refseq_transcript_canonical,
    m.genome_build,
    m.schema_version && `data schema ${m.schema_version}`,
    m.generated && `generated ${m.generated}`,
  ].filter(Boolean).join(' · ');
}

/* ===========================================================================
 * Domains
 * ======================================================================== */

/** All domains (optionally excluding provisional ones). */
export function getDomains({ includeProvisional = true } = {}, data = getLoadedData()) {
  const ds = data.domains || [];
  return includeProvisional ? ds : ds.filter((d) => !isProvisional(d));
}

/** Every domain whose [start, end] span contains `pos` (null bounds skipped). */
export function domainsForResidue(pos, data = getLoadedData()) {
  return (data.domains || []).filter(
    (d) => Number.isFinite(d.start) && Number.isFinite(d.end) && pos >= d.start && pos <= d.end
  );
}

/**
 * The most specific domain for a residue = the smallest containing span.
 * This makes the P-domain (an insertion within the palm) win over the palm.
 */
export function primaryDomainForResidue(pos, data = getLoadedData()) {
  const hits = domainsForResidue(pos, data);
  if (!hits.length) return null;
  return hits.reduce((a, b) => (b.end - b.start < a.end - a.start ? b : a));
}

/* ===========================================================================
 * Residues (active sites, metal-coordinating, motif residues)
 * ======================================================================== */

/** Filter residues by (substring) role, (substring) motif, and/or domain id. */
export function getResidues({ role, motif, domainId } = {}, data = getLoadedData()) {
  let rs = data.residues || [];
  if (role) rs = rs.filter((r) => (r.role || '').toLowerCase().includes(role.toLowerCase()));
  if (motif) rs = rs.filter((r) => (r.motif || '').toLowerCase().includes(motif.toLowerCase()));
  if (domainId) rs = rs.filter((r) => r.domain_id === domainId);
  return rs;
}


/* ===========================================================================
 * Metal sites ([4Fe-4S], two-metal pol site, exo site, Zn fingers)
 * ======================================================================== */

export const getMetalSites = (data = getLoadedData()) => data.metal_sites || [];

/* ===========================================================================
 * Variants
 * ======================================================================== */

/**
 * Filter variants.
 * @param {object} [f]
 * @param {string}  [f.domainId]
 * @param {string}  [f.status]                  exact verification_status
 * @param {string}  [f.germlineClassification]  exact ClinVar germline class
 * @param {boolean} [f.hasClinVarId]            true=only with a VCV id, false=only without
 * @param {boolean} [f.somaticOnly]             only variants with a somatic_role
 */
export function getVariants(f = {}, data = getLoadedData()) {
  const { domainId, status, germlineClassification, hasClinVarId, somaticOnly } = f;
  let vs = data.variants || [];
  if (domainId) vs = vs.filter((v) => v.domain_id === domainId);
  if (status) vs = vs.filter((v) => v.verification_status === status);
  if (germlineClassification)
    vs = vs.filter((v) => (v.clinvar_germline_classification || '') === germlineClassification);
  if (hasClinVarId === true) vs = vs.filter((v) => v.clinvar_variation_id != null);
  if (hasClinVarId === false) vs = vs.filter((v) => v.clinvar_variation_id == null);
  if (somaticOnly) vs = vs.filter((v) => !!v.somatic_role);
  return vs;
}

export const getVariantById = (id, data = getLoadedData()) =>
  (data.variants || []).find((v) => v.id === id) || null;

/** Look up by the legacy label, e.g. 'P286R', 'L424V'. */
export const getVariantByLegacy = (legacy, data = getLoadedData()) =>
  (data.variants || []).find((v) => v.legacy === legacy) || null;

export const getVariantsByCodon = (codon, data = getLoadedData()) =>
  (data.variants || []).filter((v) => v.codon === codon);

/** True for truncating / loss-of-function protein consequences (frameshift, nonsense). */
function isTruncating(v) {
  const p = (v && v.hgvs_p) || '';
  return /fs|Ter|\*\d|\*\)/.test(p);
}

/** Canonical ClinVar URL for a variant, or null if no VCV id is recorded yet. */
export const clinVarUrl = (v) =>
  v && v.clinvar_variation_id != null
    ? `https://www.ncbi.nlm.nih.gov/clinvar/variation/${v.clinvar_variation_id}/`
    : null;

/* ===========================================================================
 * Color helpers for the viewer's color modes
 *   - germline pathogenicity is one track (ClinVar germline classification)
 *   - somatic/oncogenic hotspot is a SEPARATE track. Never merge them:
 *     P286R is germline VUS but a somatic oncogenic hotspot.
 * ======================================================================== */

export const CLASSIFICATION_COLORS = {
  default: {
    Pathogenic: '#c0392b',
    'Likely pathogenic': '#e07a3f',
    'Uncertain significance': '#cfcfcf',
    'Likely benign': '#5aa0c4',
    Benign: '#2e7d4f',
    confirm_live: '#8a8a8a',
    unknown: '#8a8a8a',
  },
  // Okabe-Ito-derived, colorblind-friendly
  accessible: {
    Pathogenic: '#d55e00',
    'Likely pathogenic': '#e69f00',
    'Uncertain significance': '#999999',
    'Likely benign': '#56b4e9',
    Benign: '#009e73',
    confirm_live: '#777777',
    unknown: '#777777',
  },
};

/** Color for the germline "Pathogenicity" color mode. */
export function germlineClassColor(classification, scheme = 'default') {
  const table = CLASSIFICATION_COLORS[scheme] || CLASSIFICATION_COLORS.default;
  return table[classification] || table.unknown;
}

/** Color for the separate somatic/oncogenic-hotspot track. */
export const SOMATIC_HOTSPOT_COLOR = '#7b1fa2';

/* ===========================================================================
 * Adapter - normalized "scene model"
 *   Consume THIS from your Three.js code. It flattens the canonical records
 *   into stable, camelCase shapes and precomputes `provisional`. If your
 *   existing render functions expect different field names, map them here in
 *   ONE place rather than scattering field access through the scene code.
 * ======================================================================== */

export function buildSceneModel(data = getLoadedData()) {
  return {
    meta: getMeta(data),
    references: getReferences(data),
    provenanceFooter: provenanceFooter(data),

    domains: getDomains({}, data).map((d) => ({
      id: d.id,
      label: d.label,
      code: d.short_code,
      start: d.start,
      end: d.end,
      insertionWithin: d.is_insertion_within || null,
      provisional: isProvisional(d),
      status: d.verification_status,
      sources: (d.sources || []).slice(),
    })),

    residues: (data.residues || []).map((r) => ({
      id: r.id,
      position: r.position,
      aa: r.aa,
      role: r.role,
      motif: r.motif || null,
      domainId: r.domain_id || null,
      yeastEquivalent: r.yeast_equivalent || null,
      provisional: isProvisional(r),
      status: r.verification_status,
    })),

    metalSites: (data.metal_sites || []).map((m) => ({
      id: m.id,
      label: m.label,
      metals: m.metals || (m.metal ? [m.metal] : []),
      coordinatingResidues: m.coordinating_residues || [],
      disputed: !!m._DISPUTED_LOCATION,
      provisional: isProvisional(m),
      status: m.verification_status,
      raw: m, // full record (e.g. the dual [4Fe-4S] annotations array)
    })),

    variants: (data.variants || []).map((v) => ({
      id: v.id,
      label: v.legacy,
      codon: v.codon,
      domainId: v.domain_id || null,
      hgvsC: v.hgvs_c || null,
      hgvsP: v.hgvs_p || null,
      germline: v.clinvar_germline_classification || null,
      somatic: v.somatic_role || null,
      isSomaticHotspot: !!v.somatic_role,
      isTruncating: isTruncating(v),
      mechanism: v.mechanism || null,
      conditions: v.conditions || [],
      clinvarId: v.clinvar_variation_id != null ? v.clinvar_variation_id : null,
      clinvarUrl: clinVarUrl(v),
      investigationStatus: v.investigation_status || null,
      candidateMechanisms: v.candidate_mechanisms || [],
      clinicalObservation: v.clinical_observation || null,
      clinicalLabClassification: v.clinical_lab_classification || null,
      variantType: v.variant_type || null,
      provisional: isProvisional(v),
      status: v.verification_status,
      sources: (v.sources || []).slice(),
    })),
  };
}

/** Convenience: load + adapt in one call. */
export async function loadSceneModel(url) {
  const data = await loadPoleData(url);
  return buildSceneModel(data);
}
