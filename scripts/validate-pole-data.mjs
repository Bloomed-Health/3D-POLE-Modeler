#!/usr/bin/env node
/**
 * validate-pole-data.mjs
 * ---------------------------------------------------------------------------
 * Integrity / cross-reference validator for the POLE canonical data file.
 * Zero dependencies - runs anywhere Node >= 16 is available, offline.
 *
 * This complements JSON-Schema validation (which checks shape/types/enums).
 * Here we check the things a schema can't easily express:
 *   - unique ids within domains / residues / metal_sites / variants
 *   - every `sources[]` key resolves to an entry in `references`
 *   - every variant.domain_id and residue.domain_id resolves to a domain id
 *   - every metal_sites[].coordinating_residues id resolves to a residue id
 *   - verification_status is in the allowed set
 *   - domains with numeric bounds satisfy start <= end
 *   - clinvar_variation_id, when present, is a positive integer
 *   - hgvs_c, when present, is prefixed with the canonical transcript (WARN)
 *
 * Usage:
 *   node scripts/validate-pole-data.mjs [path/to/pole_canonical_data.json]
 * Defaults to data/pole_canonical_data.json relative to the current dir.
 *
 * Exit code 0 = OK (warnings allowed), 1 = errors found / file unreadable.
 * ---------------------------------------------------------------------------
 */

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';

const ALLOWED_STATUS = new Set([
  'verified',
  'verified_as_disputed',
  'estimated',
  'schematic',
  'todo',
  'verified_hgvs_pending_vcv',
  'hgvs_c_pending_verification',
  'hgvs_pending_verification',
]);

const errors = [];
const warnings = [];
const err = (m) => errors.push(m);
const warn = (m) => warnings.push(m);

// ---- load ----------------------------------------------------------------
const path = resolve(process.argv[2] || 'data/pole_canonical_data.json');
let data;
try {
  data = JSON.parse(readFileSync(path, 'utf8'));
} catch (e) {
  console.error(`Could not read/parse ${path}\n  ${e.message}`);
  process.exit(1);
}

// ---- helpers ---------------------------------------------------------
function checkUniqueIds(arr, name) {
  const seen = new Set();
  for (const item of arr || []) {
    const id = item && item.id;
    if (id == null) { err(`${name}: an item is missing an "id"`); continue; }
    if (seen.has(id)) err(`${name}: duplicate id "${id}"`);
    seen.add(id);
  }
  return seen;
}
function checkStatus(item, name) {
  if (!ALLOWED_STATUS.has(item.verification_status)) {
    err(`${name} "${item.id}": invalid verification_status "${item.verification_status}"`);
  }
}
function checkSources(item, name, refKeys) {
  for (const key of item.sources || []) {
    if (!refKeys.has(key)) err(`${name} "${item.id}": sources -> unknown reference key "${key}"`);
  }
}

// ---- top-level presence --------------------------------------------------
for (const key of ['meta', 'references', 'domains', 'residues', 'metal_sites', 'variants']) {
  if (data[key] == null) err(`missing top-level key "${key}"`);
}

const refKeys = new Set(Object.keys(data.references || {}));
const domainIds = checkUniqueIds(data.domains, 'domains');
const residueIds = checkUniqueIds(data.residues, 'residues');
checkUniqueIds(data.metal_sites, 'metal_sites');
checkUniqueIds(data.variants, 'variants');

// ---- meta sanity ---------------------------------------------------------
const m = data.meta || {};
if (m.refseq_transcript_canonical && !/^NM_\d+\.\d+$/.test(m.refseq_transcript_canonical)) {
  err(`meta.refseq_transcript_canonical "${m.refseq_transcript_canonical}" is not a versioned RefSeq id`);
}
const transcript = m.refseq_transcript_canonical || '';

// ---- domains -------------------------------------------------------------
for (const d of data.domains || []) {
  checkStatus(d, 'domains');
  checkSources(d, 'domains', refKeys);
  if (Number.isFinite(d.start) && Number.isFinite(d.end) && d.start > d.end) {
    err(`domains "${d.id}": start (${d.start}) > end (${d.end})`);
  }
  if (d.is_insertion_within && !domainIds.has(d.is_insertion_within)) {
    err(`domains "${d.id}": is_insertion_within -> unknown domain "${d.is_insertion_within}"`);
  }
}

// ---- residues ------------------------------------------------------------
for (const r of data.residues || []) {
  checkStatus(r, 'residues');
  checkSources(r, 'residues', refKeys);
  if (r.domain_id && !domainIds.has(r.domain_id)) {
    err(`residues "${r.id}": domain_id -> unknown domain "${r.domain_id}"`);
  }
}

// ---- metal sites ---------------------------------------------------------
for (const ms of data.metal_sites || []) {
  checkStatus(ms, 'metal_sites');
  checkSources(ms, 'metal_sites', refKeys);
  for (const rid of ms.coordinating_residues || []) {
    if (!residueIds.has(rid)) {
      err(`metal_sites "${ms.id}": coordinating_residues -> unknown residue "${rid}"`);
    }
  }
  if (ms.domain_id && !domainIds.has(ms.domain_id)) {
    err(`metal_sites "${ms.id}": domain_id -> unknown domain "${ms.domain_id}"`);
  }
}

// ---- variants ------------------------------------------------------------
let nullVcv = 0;
for (const v of data.variants || []) {
  checkStatus(v, 'variants');
  checkSources(v, 'variants', refKeys);
  if (v.domain_id && !domainIds.has(v.domain_id)) {
    err(`variants "${v.id}": domain_id -> unknown domain "${v.domain_id}"`);
  }
  if (v.clinvar_variation_id != null) {
    if (!Number.isInteger(v.clinvar_variation_id) || v.clinvar_variation_id <= 0) {
      err(`variants "${v.id}": clinvar_variation_id must be a positive integer`);
    }
  } else {
    nullVcv += 1;
  }
  if (typeof v.hgvs_c === 'string' && transcript && !v.hgvs_c.startsWith(transcript + ':')) {
    warn(`variants "${v.id}": hgvs_c "${v.hgvs_c}" is not prefixed with ${transcript}:`);
  }
  // Verified variants should carry a ClinVar id (soft expectation).
  if (v.verification_status === 'verified' && v.clinvar_variation_id == null) {
    warn(`variants "${v.id}": status "verified" but no clinvar_variation_id`);
  }
  // Citation-required: variants with mechanism/conditions but empty sources
  if ((v.mechanism || (v.conditions && v.conditions.length)) && (!v.sources || !v.sources.length)) {
    warn(`variants "${v.id}": has mechanism or conditions but empty sources array`);
  }
  // Citation-required: has germline classification but no ClinVar id
  if (v.clinvar_germline_classification && v.clinvar_germline_classification !== 'confirm_live' && v.clinvar_variation_id == null) {
    warn(`variants "${v.id}": has clinvar_germline_classification "${v.clinvar_germline_classification}" but no clinvar_variation_id`);
  }
  // Investigation status consistency: if set, candidate_mechanisms should be non-empty
  if (v.investigation_status && (!v.candidate_mechanisms || !v.candidate_mechanisms.length)) {
    warn(`variants "${v.id}": investigation_status "${v.investigation_status}" but no candidate_mechanisms`);
  }
}

// ---- report --------------------------------------------------------------
const counts = {
  references: refKeys.size,
  domains: (data.domains || []).length,
  residues: (data.residues || []).length,
  metal_sites: (data.metal_sites || []).length,
  variants: (data.variants || []).length,
};
console.log(`POLE data integrity check - ${path}`);
console.log(
  `  records: ${counts.domains} domains, ${counts.residues} residues, ` +
    `${counts.metal_sites} metal sites, ${counts.variants} variants, ${counts.references} references`
);
console.log(`  variants still awaiting a ClinVar VCV id: ${nullVcv}`);

if (warnings.length) {
  console.log(`\n${warnings.length} warning(s):`);
  for (const w of warnings) console.log(`  - ${w}`);
}
if (errors.length) {
  console.error(`\n${errors.length} error(s):`);
  for (const e of errors) console.error(`  - ${e}`);
  process.exit(1);
}
console.log('\nNo integrity errors.');
