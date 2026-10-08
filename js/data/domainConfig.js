// ==========================================================================
// Shared POLE Domain Configuration (synchronous)
//
// Synchronous domain/mutation config for viewers that need static imports at
// module load time. Values MUST match data/pole_canonical_data.json.
// For async metadata (provenance, verification status, references), use
// poleData.js which loads from the canonical JSON.
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

// Domain boundary convention: Atlas of Genetics / UniProt (broad) ranges.
// Matches pole_canonical_data.json and pipeline/src/config.jl.
// Alternative convention (Church et al. 2013, narrow EXO 268-471) is not used here.
export const DOMAINS = {
  ntd:      { range: [1, 222],    color: '#c4a47a', name: 'N-terminal domain',     abbrev: 'NTD',
    tag: 'scaffold',
    desc: 'Structural scaffolding region preceding the catalytic core. Contains residues affected by the ultra-rare c.138del germline truncation.' },
  exo:      { range: [223, 517],  color: '#b85c5c', name: "3'\u21925' Exonuclease", abbrev: 'EXO',
    tag: 'proofreading',
    desc: 'Proofreading domain that excises misincorporated nucleotides. Two divalent metal ions are coordinated by D275 and E277 (ExoI) and D462 (ExoIII; verified in the 9F6L mismatch-excision structure). D368 (ExoII) is structural / DNA-binding rather than metal-coordinating. Primary site of cancer-driving POLE mutations.' },
  palm:     { range: [518, 870],  color: '#6b8e6b', name: 'Palm subdomain',        abbrev: 'PALM',
    tag: 'catalytic core',
    desc: 'Houses the polymerase active site. Conserved aspartate carboxylates (D640, D642 motif A; D860 motif C) coordinate two Mg\u00B2\u207A ions for phosphodiester bond formation via the two-metal-ion catalytic mechanism.' },
  pdomain:  { range: [640, 760],  color: '#4a8a7a', name: 'Processivity domain',   abbrev: 'P-DOM',
    tag: '[4Fe-4S] \u00B7 PCNA',
    desc: 'Pol \u03B5-specific insertion between Palm and Fingers that distinguishes it from all other B-family polymerases. Houses the [4Fe-4S] iron-sulfur cluster (CysX motif) essential for polymerase activity, contacts the nascent duplex, and forms one of three PCNA contact surfaces.' },
  fingers:  { range: [870, 1020], color: '#5b7f99', name: 'Fingers subdomain',     abbrev: 'FNG',
    tag: 'fidelity gate',
    desc: 'Predominantly \u03B1-helical region. The O-helix undergoes a 28\u00B0/19 \u00C5 conformational closure upon dNTP binding (open \u2192 ajar \u2192 closed; PDB 9F6D/9F6E/9F6F), ensuring Watson\u2013Crick base-pair geometry. Rate-limiting step for nucleotide incorporation and the structural basis of replicative fidelity.' },
  thumb:    { range: [1020,1189], color: '#8b7ba8', name: 'Thumb subdomain',       abbrev: 'THM',
    tag: 'duplex grip',
    desc: 'Helical bundle that grips the minor groove of the nascent duplex DNA product. Contacts extend ~8 bp into the product duplex. The Thumb insertion (res. 1102\u20131122) forms a second PCNA contact, contributing to the tripartite interface that makes Pol \u03B5\u2013PCNA resistant to displacement by Pol \u03B4 and FEN1.' },
  inactpol: { range: [1190,1900], color: '#7a6a8a', name: 'Inactive Pol module',   abbrev: 'INACT',
    tag: 'subunit scaffold',
    desc: 'Non-catalytic B-family polymerase fold (the \u201CPol2 module\u201D of Yuan et al. 2020) that retains the canonical Palm/Fingers/Thumb architecture but has lost all enzymatic activity. Serves as the structural platform for POLE2/3/4 accessory subunit assembly and bridges the catalytic N-terminal half to the CTD zinc-finger appendix.' },
  ctd:      { range: [1900,2286], color: '#b8a04a', name: 'C-terminal domain',     abbrev: 'CTD',
    tag: 'Zn-finger',
    desc: 'Contains CysA and CysB zinc-finger motifs that serve as the structural platform for binding accessory subunits POLE2, POLE3, and POLE4, completing the heterotetrameric holoenzyme.' },
};

export const DOMAIN_ORDER = ['ntd','exo','palm','pdomain','fingers','thumb','inactpol','ctd'];

/** Array form of DOMAINS for consumers that need iteration with key access. */
export const DOMAINS_LIST = DOMAIN_ORDER.map(key => ({ key, ...DOMAINS[key] }));

export const MUTATIONS = [
  { id: 'P286R', residue: 286, domain: 'exo', label: 'Pro286\u2192Arg',
    origins: ['Somatic'], loc: 'EXO \u00B7 ultra-mutator',
    detail: 'Somatic ultra-mutator / adjacent to ExoI / CRC & endometrial / SBS10a/b',
    text: 'Most frequent cancer-associated POLE mutation. Pro\u2192Arg substitution adjacent to ExoI disrupts the proofreading active site, causing an ultra-mutation phenotype (TMB >100 mut/Mb) in colorectal & endometrial carcinoma.' },
  { id: 'V411L', residue: 411, domain: 'exo', label: 'Val411\u2192Leu',
    origins: ['Somatic'], loc: 'EXO \u00B7 proofreading-deficient',
    detail: 'Somatic proofreading-deficient / TMB >100 mut/Mb / SBS10a',
    text: 'Recurrent exonuclease-domain mutation that abolishes editing capacity. Associated with microsatellite-stable tumours with very high mutational burden, often exceeding 100 mut/Mb.' },
  { id: 'S297F', residue: 297, domain: 'exo', label: 'Ser297\u2192Phe',
    origins: ['Somatic', 'Germline'], loc: 'EXO \u00B7 germline PPAP',
    detail: 'Germline ExoI / PPAP predisposition / SBS28',
    text: 'ExoI-region variant; originally reported as somatic endometrial variant (Church et al. 2013). Proofreading-impairing; germline PPAP association requires further pedigree evidence.' },
  { id: 'L424V', residue: 424, domain: 'exo', label: 'Leu424\u2192Val',
    origins: ['Germline'], loc: 'EXO \u00B7 germline PPAP founder',
    detail: 'Germline PPAP founder variant (Palles 2013) / ~50% PPAP families / Note: somatic L424I is a distinct variant',
    text: 'Germline PPAP founder variant (Palles 2013) / ~50% PPAP families / Note: somatic L424I is a distinct variant.' },
  { id: 'D287E', residue: 287, domain: 'exo', label: 'Asp287\u2192Glu',
    origins: ['Somatic'], loc: 'EXO \u00B7 somatic recurrence',
    detail: 'Somatic ExoII adjacent / CRC / SBS10b',
    text: 'Somatic ExoII-adjacent recurrence \u00B7 CRC \u00B7 SBS10b. Not in the Leon-Castillo 2020 pathogenic set.' },
  { id: 'P436R', residue: 436, domain: 'exo', label: 'Pro436\u2192Arg',
    origins: ['Somatic'], loc: 'EXO \u00B7 somatic ExoIII',
    detail: 'Somatic ExoIII / endometrial / SBS10a',
    text: 'Somatic ExoIII / endometrial / SBS10a.' },
  { id: 'M444K', residue: 444, domain: 'exo', label: 'Met444\u2192Lys',
    origins: ['Germline'], loc: 'EXO \u00B7 germline PPAP',
    detail: 'Germline PPAP / exonuclease domain',
    text: 'Germline PPAP / exonuclease domain.' },
  { id: 'S459F', residue: 459, domain: 'exo', label: 'Ser459\u2192Phe',
    origins: ['Somatic', 'Germline'], loc: 'EXO \u00B7 somatic EDM + germline',
    detail: 'Germline / PPAP-associated / Valle 2023',
    text: 'Somatic EDM in the Leon-Castillo 2020 pathogenic set; also reported as a germline PPAP-associated variant (Valle 2023).' },
  { id: 'F367S', residue: 367, domain: 'exo', label: 'Phe367\u2192Ser',
    origins: ['Germline'], loc: 'EXO \u00B7 germline PPAP',
    detail: 'Germline PPAP / adjacent to ExoII structural residue D368',
    text: 'Adjacent to ExoII structural residue D368. Germline variant associated with polymerase proofreading-associated polyposis.' },
  { id: 'A456P', residue: 456, domain: 'exo', label: 'Ala456\u2192Pro',
    origins: ['Somatic'], loc: 'EXO \u00B7 somatic EDM',
    detail: 'Somatic EDM / Leon-Castillo 2020 pathogenic set / endometrial',
    text: 'Somatic EDM / Leon-Castillo 2020 pathogenic set / endometrial.' },
];
