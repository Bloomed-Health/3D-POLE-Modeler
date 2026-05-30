// ==========================================================================
// POLE Viewer Palette & Color Schemes
//
// Shared THREE.Color palette and color scheme definitions for the
// procedural ribbon viewers (pole-viewer.js, mutation-viewer.js).
// ==========================================================================

import * as THREE from 'three';
import { DOMAIN_COLORS } from './data/domainConfig.js';

// ==================== PALETTE (muted, publication-quality) ====================

export const PAL = {
  // Domain colors derived from shared config
  ntd:      new THREE.Color(DOMAIN_COLORS.ntd),
  exo:      new THREE.Color(DOMAIN_COLORS.exo),
  palm:     new THREE.Color(DOMAIN_COLORS.palm),
  pdomain:  new THREE.Color(DOMAIN_COLORS.pdomain),
  fingers:  new THREE.Color(DOMAIN_COLORS.fingers),
  thumb:    new THREE.Color(DOMAIN_COLORS.thumb),
  inactpol: new THREE.Color(DOMAIN_COLORS.inactpol),
  ctd:      new THREE.Color(DOMAIN_COLORS.ctd),
  linker:   new THREE.Color(0x464040),

  // DNA
  dnaTemplate: new THREE.Color(0x8c9aa8),
  dnaPrimer:   new THREE.Color(0x6a9e96),
  dnaBase:     new THREE.Color(0x5a6a78),

  // Atomic element colors
  carbon:   new THREE.Color(0x707070),
  nitrogen: new THREE.Color(0x3050c0),
  oxygen:   new THREE.Color(0xc03030),
  phosphorus: new THREE.Color(0xd07020),
  magnesium:  new THREE.Color(0x60a830),

  // UI
  mutMarker: 0xc45c5c,
  accent:    0x7abba8,
  bg:        0x080b11,
  annotate:  new THREE.Color(0x506070),
};

// ==================== COLOR SCHEMES ====================

export const COLOR_SCHEMES = {
  default:    { ntd:0xc4a47a, exo:0xb85c5c, palm:0x6b8e6b, pdomain:0x4a8a7a, fingers:0x5b7f99, thumb:0x8b7ba8, inactpol:0x7a6a8a, ctd:0xb8a04a },
  pymol:      { ntd:0x33cc33, exo:0xff4444, palm:0x3366ff, pdomain:0x00aa88, fingers:0xff66cc, thumb:0x00cccc, inactpol:0xaa66ff, ctd:0xffcc00 },
  colorblind: { ntd:0xe69f00, exo:0x0072b2, palm:0x009e73, pdomain:0x44aa99, fingers:0xcc79a7, thumb:0x56b4e9, inactpol:0xaa4499, ctd:0xf0e442 },
};
