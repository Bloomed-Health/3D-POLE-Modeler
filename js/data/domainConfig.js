// ==========================================================================
// Shared POLE Domain Configuration
//
// Single source of truth for domain boundaries, colors, ordering, and
// canonical mutation hotspots. Consumed by pole-viewer.js, molstar-viewer.js,
// mutation-viewer.js, and the molstar.html inline script.
//
// DOMAIN BOUNDARY NOTE: These ranges are approximate and have not been
// derived from DSSP analysis of deposited coordinates.
// See: data/pole_canonical_data.json for the canonical values.
// ==========================================================================

export const DOMAIN_COLORS = {
  ntd:      '#c4a47a',
  exo:      '#b85c5c',
  palm:     '#6b8e6b',
  pdomain:  '#4a8a7a',
  fingers:  '#5b7f99',
  thumb:    '#8b7ba8',
  inactpol: '#7a6a8a',
  ctd:      '#b8a04a',
};

export const DOMAINS = {
  ntd:      { range: [1, 280],    color: '#c4a47a', name: 'N-terminal domain',     abbrev: 'NTD' },
  exo:      { range: [268, 471],  color: '#b85c5c', name: "3'\u21925' Exonuclease", abbrev: 'EXO' },
  palm:     { range: [472, 699],  color: '#6b8e6b', name: 'Palm subdomain',        abbrev: 'PALM' },
  pdomain:  { range: [700, 760],  color: '#4a8a7a', name: 'Processivity domain',   abbrev: 'P-DOM' },
  fingers:  { range: [761, 1020], color: '#5b7f99', name: 'Fingers subdomain',     abbrev: 'FNG' },
  thumb:    { range: [1021,1190], color: '#8b7ba8', name: 'Thumb subdomain',       abbrev: 'THM' },
  inactpol: { range: [1191,1900], color: '#7a6a8a', name: 'Inactive Pol module',   abbrev: 'INACT' },
  ctd:      { range: [1901,2286], color: '#b8a04a', name: 'C-terminal domain',     abbrev: 'CTD' },
};

export const DOMAIN_ORDER = ['ntd','exo','palm','pdomain','fingers','thumb','inactpol','ctd'];

export const MUTATIONS = [
  { id: 'P286R', residue: 286, domain: 'exo', label: 'Pro286\u2192Arg',
    detail: 'Somatic ultra-mutator / adjacent to ExoI / CRC & endometrial / SBS10a/b' },
  { id: 'V411L', residue: 411, domain: 'exo', label: 'Val411\u2192Leu',
    detail: 'Somatic proofreading-deficient / TMB >100 mut/Mb / SBS10a' },
  { id: 'S297F', residue: 297, domain: 'exo', label: 'Ser297\u2192Phe',
    detail: 'Germline ExoI / PPAP predisposition / SBS28' },
  { id: 'L424V', residue: 424, domain: 'exo', label: 'Leu424\u2192Val',
    detail: 'Germline ExoIII / PPAP-associated / Mur 2023' },
  { id: 'D287E', residue: 287, domain: 'exo', label: 'Asp287\u2192Glu',
    detail: 'Somatic ExoII adjacent / CRC / SBS10b' },
  { id: 'P436R', residue: 436, domain: 'exo', label: 'Pro436\u2192Arg',
    detail: 'Somatic ExoIII / endometrial / SBS10a' },
  { id: 'M444K', residue: 444, domain: 'exo', label: 'Met444\u2192Lys',
    detail: 'Somatic / ExoIII motif / CRC / SBS10b' },
  { id: 'S459F', residue: 459, domain: 'exo', label: 'Ser459\u2192Phe',
    detail: 'Germline / PPAP-associated / Valle 2023' },
  { id: 'F367S', residue: 367, domain: 'exo', label: 'Phe367\u2192Ser',
    detail: 'Germline PPAP / adjacent to ExoIII catalytic triad' },
];
