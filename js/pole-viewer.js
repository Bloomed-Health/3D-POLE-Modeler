// ==========================================================================
// POLE p261 -3D Structural Viewer
// Cartoon Ribbon Representation
//
// The protein fold is rendered as cartoon ribbons -helices as wide flat
// spirals, strands as flat arrows with arrowheads, loops as thin tubes.
// DNA template/primer rendered as backbone tubes with base-pair rungs.
//
// References:
//   UniProt Q07864 (POLE_HUMAN)
//   PDB: 9F6D / 9F6E / 9F6F (Roske & Yeeles 2024), 9B8S (He et al. 2024),
//        4M8O (Hogg et al. 2014), 6WJV (Lancey et al. 2020, Pol delta for comparison)
//   Roske & Yeeles, Nat Struct Mol Biol 31:1921 (2024) — human Pol ε structures
//   He, Wang, Yao, O'Donnell & Li, Nat Commun 15:7847 (2024) — human Pol ε–PCNA
//   Mur et al., Genome Med 15:85 (2023) — POLE/POLD1 ACMG classification
//   Yuan, Georgescu, Schauer, O'Donnell & Li, Nat Commun 11:3156 (2020) — yeast Pol ε
//   Robinson, Coorens et al., Nat Genet 53:1434 (2021) — normal-tissue mutation rates
//   Hogg et al., Nat Struct Mol Biol 21:49 (2014) — S. cerevisiae Pol2 NTD (4M8O)
//
// NOTE: This viewer renders a procedural cartoon-ribbon schematic, not atomic
// coordinates from PDB depositions. Color modes (B-factor, conservation, charge,
// pathogenicity) use approximate domain-level or synthetic per-residue values.
// For coordinate-derived analysis, load PDB 9F6D into Mol* or PyMOL.
// ==========================================================================

import * as THREE from 'three';
import { buildSSElement, segmentBySSType } from './ribbon.js';


// ==================== PALETTE (muted, publication-quality) ====================

const PAL = {
  ntd:      new THREE.Color(0xc4a47a),
  exo:      new THREE.Color(0xb85c5c),
  palm:     new THREE.Color(0x6b8e6b),
  pdomain:  new THREE.Color(0x4a8a7a),
  fingers:  new THREE.Color(0x5b7f99),
  thumb:    new THREE.Color(0x8b7ba8),
  inactpol: new THREE.Color(0x7a6a8a),
  ctd:      new THREE.Color(0xb8a04a),
  linker:   new THREE.Color(0x464040),

  dnaTemplate: new THREE.Color(0x8c9aa8),
  dnaPrimer:   new THREE.Color(0x6a9e96),
  dnaBase:     new THREE.Color(0x5a6a78),

  carbon:   new THREE.Color(0x707070),
  nitrogen: new THREE.Color(0x3050c0),
  oxygen:   new THREE.Color(0xc03030),
  phosphorus: new THREE.Color(0xd07020),
  magnesium:  new THREE.Color(0x60a830),

  mutMarker: 0xc45c5c,
  accent:    0x7abba8,
  bg:        0x080b11,
  annotate:  new THREE.Color(0x506070),
};

// ==================== COLOR SCHEMES ====================

const COLOR_SCHEMES = {
  default:    { ntd:0xc4a47a, exo:0xb85c5c, palm:0x6b8e6b, pdomain:0x4a8a7a, fingers:0x5b7f99, thumb:0x8b7ba8, inactpol:0x7a6a8a, ctd:0xb8a04a },
  pymol:      { ntd:0x33cc33, exo:0xff4444, palm:0x3366ff, pdomain:0x00aa88, fingers:0xff66cc, thumb:0x00cccc, inactpol:0xaa66ff, ctd:0xffcc00 },
  colorblind: { ntd:0xe69f00, exo:0x0072b2, palm:0x009e73, pdomain:0x44aa99, fingers:0xcc79a7, thumb:0x56b4e9, inactpol:0xaa4499, ctd:0xf0e442 },
};

// ==================== PIPELINE DATA (defaults, overridden by JSON) ====================

// B-factor pseudo-values by SS type (0=rigid, 1=flexible)
let BFACTOR_MAP = { H: 0.2, E: 0.5, L: 0.9 };

// Conservation pseudo-scores by domain (0=variable, 1=conserved)
let CONSERVATION = { ntd:0.3, exo:0.85, palm:0.95, pdomain:0.8, fingers:0.7, thumb:0.6, inactpol:0.35, ctd:0.4 };

// Variant density pseudo-scores by domain (COSMIC/ClinVar)
let VARIANT_DENSITY = { ntd:0.1, exo:0.95, palm:0.3, pdomain:0.15, fingers:0.15, thumb:0.1, inactpol:0.05, ctd:0.05 };

// Per-residue scoring arrays (populated from pipeline JSON when available)
let PER_RESIDUE_SCORES = null;  // { bfactor:[], conservation:[], variant_density:[], pathogenicity:[], electrostatic_charge:[] }

// Mutation data with pipeline annotations (populated from JSON when available)
let PIPELINE_MUTATIONS = null;  // Full mutation data with signatures, classification, DDG

// Pathogenicity spectrum per-domain: residue-range → score (0=benign, 1=pathogenic)
let PATHOGENICITY = {
  ntd: [ { start:1, end:280, score:0.15 } ],
  exo: [
    { start:268, end:280, score:0.3 },
    { start:281, end:285, score:0.6 },
    { start:286, end:286, score:1.0 },   // P286R hotspot
    { start:287, end:296, score:0.5 },
    { start:297, end:297, score:0.9 },   // S297F hotspot
    { start:298, end:367, score:0.4 },
    { start:368, end:370, score:0.75 },  // EXO catalytic D368/D370
    { start:371, end:410, score:0.35 },
    { start:411, end:411, score:0.95 },  // V411L hotspot
    { start:412, end:471, score:0.3 },
  ],
  palm: [
    { start:600, end:639, score:0.2 },
    { start:640, end:642, score:0.85 },  // catalytic D640/D642
    { start:643, end:859, score:0.25 },
    { start:860, end:860, score:0.8 },   // catalytic D860
    { start:861, end:870, score:0.2 },
  ],
  pdomain: [ { start:700, end:760, score:0.2 } ],
  fingers: [ { start:870, end:1020, score:0.15 } ],
  thumb:   [ { start:1020, end:1190, score:0.1 } ],
  inactpol:[ { start:1198, end:1900, score:0.05 } ],
  ctd:     [ { start:1900, end:2286, score:0.05 } ],
};

function getPathogenicityScore(domain, residueIdx, totalResidues) {
  const dom = DOMAINS[domain];
  if (!dom) return 0;
  const residue = dom.range[0] + Math.round((residueIdx / Math.max(1, totalResidues - 1)) * (dom.range[1] - dom.range[0]));
  const ranges = PATHOGENICITY[domain];
  if (!ranges) return 0;
  for (const r of ranges) {
    if (residue >= r.start && residue <= r.end) return r.score;
  }
  return 0;
}

// ==================== UTILITIES ====================

function makeRng(seed) {
  let s = seed;
  return () => { s = (s * 16807) % 2147483647; return (s & 0x7fffffff) / 2147483647; };
}

// ==================== DOMAIN ARCHITECTURE ====================

const SCALE = 0.5;

const DOMAINS = {
  ntd: { center:[-18,8,-5], range:[1,280], name:'N-terminal domain', abbrev:'NTD', color:'ntd',
    header:'// NTD (1-280) // Structural scaffold',
    ss:[ {t:'H',n:14},{t:'L',n:6},{t:'E',n:7},{t:'L',n:4},{t:'E',n:6},{t:'L',n:5},{t:'H',n:18},{t:'L',n:5},{t:'E',n:8},{t:'L',n:4},{t:'E',n:7},{t:'L',n:6},{t:'H',n:12},{t:'L',n:8},{t:'H',n:10},{t:'L',n:5},{t:'E',n:6},{t:'L',n:4},{t:'H',n:16},{t:'L',n:7},{t:'H',n:10},{t:'L',n:8},{t:'E',n:5},{t:'L',n:6},{t:'H',n:8},{t:'L',n:9} ] },
  exo: { center:[-10,2,7], range:[268,471], name:"3\u2032\u21925\u2032 Exonuclease", abbrev:'EXO', color:'exo',
    header:"// EXO (268-471) // 3'\u21925' proofreading / D275/E277/D368/D370",
    ss:[ {t:'E',n:6},{t:'L',n:4},{t:'H',n:14},{t:'L',n:5},{t:'E',n:7},{t:'L',n:3},{t:'E',n:6},{t:'H',n:12},{t:'L',n:6},{t:'H',n:16},{t:'L',n:4},{t:'E',n:8},{t:'L',n:5},{t:'H',n:10},{t:'L',n:6},{t:'E',n:7},{t:'L',n:3},{t:'H',n:14},{t:'L',n:5},{t:'H',n:8},{t:'L',n:6},{t:'E',n:6},{t:'L',n:4},{t:'H',n:10},{t:'L',n:4} ] },
  palm: { center:[0,-1,0], range:[600,870], name:'Palm subdomain', abbrev:'PALM', color:'palm',
    header:'// PALM (600-870) // Catalytic core / motifs A(D640xD642) C(D860)',
    ss:[ {t:'E',n:8},{t:'L',n:3},{t:'E',n:7},{t:'L',n:5},{t:'H',n:12},{t:'L',n:4},{t:'E',n:9},{t:'L',n:3},{t:'E',n:7},{t:'L',n:6},{t:'H',n:10},{t:'L',n:4},{t:'E',n:8},{t:'L',n:5},{t:'E',n:6},{t:'L',n:3},{t:'H',n:14},{t:'L',n:5},{t:'E',n:7},{t:'L',n:4},{t:'H',n:8},{t:'L',n:6},{t:'E',n:9},{t:'L',n:3},{t:'E',n:7},{t:'L',n:5},{t:'H',n:12},{t:'L',n:4} ] },
  pdomain: { center:[4,1,-3], range:[700,760], name:'Processivity (P) domain', abbrev:'P-DOM', color:'pdomain',
    header:'// P-DOM (700-760) // [4Fe-4S] cluster / PCNA contact',
    ss:[ {t:'H',n:12},{t:'L',n:4},{t:'H',n:10},{t:'L',n:5},{t:'E',n:6},{t:'L',n:3},{t:'H',n:8},{t:'L',n:4} ] },
  fingers: { center:[7,5,4], range:[870,1020], name:'Fingers subdomain', abbrev:'FNG', color:'fingers',
    header:'// FNG (870-1020) // O-helix / Watson-Crick fidelity gate',
    ss:[ {t:'H',n:20},{t:'L',n:4},{t:'H',n:16},{t:'L',n:6},{t:'H',n:22},{t:'L',n:5},{t:'H',n:14},{t:'L',n:4},{t:'H',n:12},{t:'L',n:5},{t:'H',n:10},{t:'L',n:6},{t:'H',n:8},{t:'L',n:4} ] },
  thumb: { center:[3,-9,7], range:[1020,1190], name:'Thumb subdomain', abbrev:'THM', color:'thumb',
    header:'// THM (1020-1190) // Helical bundle / grips duplex DNA',
    ss:[ {t:'H',n:18},{t:'L',n:4},{t:'H',n:22},{t:'L',n:5},{t:'H',n:16},{t:'L',n:6},{t:'H',n:20},{t:'L',n:4},{t:'H',n:14},{t:'L',n:5},{t:'H',n:12},{t:'L',n:6},{t:'H',n:10},{t:'L',n:4} ] },
  inactpol: { center:[10,-4,1], range:[1198,1900], name:'Inactive Pol module', abbrev:'INACT', color:'inactpol',
    header:'// INACT (1198-1900) // Non-catalytic Pol fold / subunit scaffold',
    ss:[ {t:'H',n:14},{t:'L',n:6},{t:'E',n:7},{t:'L',n:4},{t:'H',n:10},{t:'L',n:5},{t:'E',n:6},{t:'L',n:4},{t:'H',n:12},{t:'L',n:6},{t:'E',n:8},{t:'L',n:3},{t:'H',n:8},{t:'L',n:5},{t:'E',n:7},{t:'L',n:4},{t:'H',n:10},{t:'L',n:6},{t:'E',n:6},{t:'L',n:4},{t:'H',n:14},{t:'L',n:5},{t:'E',n:5},{t:'L',n:4},{t:'H',n:8},{t:'L',n:6} ] },
  ctd: { center:[16,1,-5], range:[1900,2286], name:'C-terminal domain', abbrev:'CTD', color:'ctd',
    header:'// CTD (1900-2286) // CysA/CysB Zn-finger / POLE2/3/4 anchor',
    ss:[ {t:'H',n:10},{t:'L',n:6},{t:'E',n:5},{t:'L',n:4},{t:'E',n:6},{t:'L',n:3},{t:'H',n:14},{t:'L',n:5},{t:'E',n:7},{t:'L',n:4},{t:'E',n:5},{t:'L',n:6},{t:'H',n:12},{t:'L',n:4},{t:'E',n:6},{t:'L',n:5},{t:'H',n:10},{t:'L',n:8},{t:'E',n:5},{t:'L',n:4},{t:'H',n:8},{t:'L',n:6},{t:'E',n:7},{t:'L',n:3},{t:'H',n:12},{t:'L',n:5},{t:'E',n:6},{t:'L',n:4},{t:'H',n:10},{t:'L',n:7},{t:'E',n:5},{t:'L',n:5},{t:'H',n:14},{t:'L',n:6} ] },
};

const DOMAIN_ORDER = ['ntd','exo','palm','pdomain','fingers','thumb','inactpol','ctd'];

let MUTATIONS = [
  { id:'P286R', residue:286, domain:'exo', offset:[2.5,1.2,2.0],
    label:'Pro286\u2192Arg', detail:'Somatic ultra-mutator / ExoII motif / CRC & endometrial / SBS10a/b' },
  { id:'V411L', residue:411, domain:'exo', offset:[-2.8,-0.5,1.8],
    label:'Val411\u2192Leu', detail:'Somatic proofreading-deficient / TMB >100 mut/Mb / SBS10a' },
  { id:'S297F', residue:297, domain:'exo', offset:[0.5,2.8,-1.2],
    label:'Ser297\u2192Phe', detail:'Germline ExoI / PPAP predisposition / SBS28' },
  { id:'L424V', residue:424, domain:'exo', offset:[-1.5,2.0,2.5],
    label:'Leu424\u2192Val', detail:'Germline ExoIII / PPAP-associated / Mur 2023' },
  { id:'D287E', residue:287, domain:'exo', offset:[3.0,-1.0,1.0],
    label:'Asp287\u2192Glu', detail:'Somatic ExoII adjacent / CRC / SBS10b' },
  { id:'P436R', residue:436, domain:'exo', offset:[-2.0,1.5,-1.5],
    label:'Pro436\u2192Arg', detail:'Somatic ExoIII / endometrial / SBS10a' },
  { id:'M444K', residue:444, domain:'exo', offset:[1.0,-2.0,2.5],
    label:'Met444\u2192Lys', detail:'Somatic / ExoIII motif / CRC / SBS10b' },
  { id:'S459F', residue:459, domain:'exo', offset:[-3.0,0.5,-0.5],
    label:'Ser459\u2192Phe', detail:'Germline / PPAP-associated / Valle 2023' },
  { id:'F367S', residue:367, domain:'exo', offset:[2.0,0.8,-2.0],
    label:'Phe367\u2192Ser', detail:'Germline PPAP / adjacent to ExoIII catalytic triad' },
];

// ==================== PIPELINE DATA LOADER ====================

/**
 * Load pipeline-generated JSON data. Falls back to hardcoded defaults on failure.
 * Call before instantiating POLEViewer for pipeline-enhanced data.
 */
export async function loadPipelineData(basePath = './data') {
  try {
    const [structure, scores, mutations] = await Promise.all([
      fetch(`${basePath}/pole_structure.json`).then(r => { if (!r.ok) throw new Error(r.status); return r.json(); }),
      fetch(`${basePath}/pole_scores.json`).then(r => { if (!r.ok) throw new Error(r.status); return r.json(); }),
      fetch(`${basePath}/pole_mutations.json`).then(r => { if (!r.ok) throw new Error(r.status); return r.json(); }),
    ]);

    // Apply structure data
    if (structure?.domains) {
      for (const [key, data] of Object.entries(structure.domains)) {
        if (DOMAINS[key]) {
          if (data.center) DOMAINS[key].center = data.center;
          if (data.range) DOMAINS[key].range = data.range;
          if (data.ss && data.ss.length > 0) DOMAINS[key].ss = data.ss;
        }
      }
    }
    if (structure?.bfactor_map) {
      BFACTOR_MAP = { ...BFACTOR_MAP, ...structure.bfactor_map };
    }

    // Apply per-residue scores
    if (scores?.per_residue) {
      PER_RESIDUE_SCORES = scores.per_residue;
    }
    if (scores?.per_domain) {
      if (scores.per_domain.conservation) {
        CONSERVATION = { ...CONSERVATION, ...scores.per_domain.conservation };
      }
      if (scores.per_domain.variant_density) {
        VARIANT_DENSITY = { ...VARIANT_DENSITY, ...scores.per_domain.variant_density };
      }
      if (scores.per_domain.electrostatic_charge) {
        // Store for electrostatic mode
        CONSERVATION._charge = scores.per_domain.electrostatic_charge;
      }
    }

    // Apply mutation data
    if (mutations?.mutations) {
      PIPELINE_MUTATIONS = mutations.mutations;
      // Enrich existing MUTATIONS array with pipeline data
      for (const pm of mutations.mutations) {
        const existing = MUTATIONS.find(m => m.id === pm.id);
        if (existing) {
          existing.pathogenicity_score = pm.pathogenicity_score;
          existing.ddg = pm.ddg_kcal_mol;
          existing.signature = pm.signature_attribution;
          existing.classification = pm.classification;
        }
      }
    }

    console.info('[POLEPipeline] Loaded pipeline data successfully');
    return { structure, scores, mutations };
  } catch (e) {
    console.warn('[POLEPipeline] Pipeline data not available, using defaults:', e.message);
    return null;
  }
}


// ==================== BACKBONE GENERATION ====================

function generateDomainBackbone(key) {
  const dom = DOMAINS[key];
  const center = new THREE.Vector3(...dom.center);
  const rng = makeRng(key.charCodeAt(0) * 137 + key.charCodeAt(1) * 31);
  const points = [];
  let cur = center.clone().add(new THREE.Vector3((rng()-0.5)*5,(rng()-0.5)*5,(rng()-0.5)*5));
  let dir = new THREE.Vector3(rng()-0.5,rng()-0.5,rng()-0.5).normalize();

  for (const el of dom.ss) {
    if (el.t === 'H') {
      const axis = dir.clone();
      let pA = new THREE.Vector3(1,0,0);
      if (Math.abs(axis.dot(pA)) > 0.85) pA.set(0,1,0);
      const pB = new THREE.Vector3().crossVectors(axis,pA).normalize();
      pA.crossVectors(pB,axis).normalize();
      // Handedness check: ensure right-handed winding
      const check = new THREE.Vector3().crossVectors(pA, pB);
      if (check.dot(axis) < 0) pB.negate();
      for (let i = 0; i < el.n; i++) {
        const ang = i * (2*Math.PI/3.6);
        const p = cur.clone()
          .add(axis.clone().multiplyScalar(i * 1.5 * SCALE))
          .add(pA.clone().multiplyScalar(Math.cos(ang) * 2.3 * SCALE))
          .add(pB.clone().multiplyScalar(Math.sin(ang) * 2.3 * SCALE));
        points.push({ pos:p, ss:'H' });
      }
    } else if (el.t === 'E') {
      // Proper 120° pleated-sheet rotation using perpendicular frame
      const sDir = dir.clone().normalize();
      let sA = new THREE.Vector3(1,0,0);
      if (Math.abs(sDir.dot(sA)) > 0.85) sA.set(0,1,0);
      const sB = new THREE.Vector3().crossVectors(sDir, sA).normalize();
      sA.crossVectors(sB, sDir).normalize();
      for (let i = 0; i < el.n; i++) {
        const ang = i * (2 * Math.PI / 3);
        const p = cur.clone()
          .add(sDir.clone().multiplyScalar(i * 3.3 * SCALE))
          .add(sA.clone().multiplyScalar(Math.cos(ang) * 0.7 * SCALE))
          .add(sB.clone().multiplyScalar(Math.sin(ang) * 0.7 * SCALE));
        points.push({ pos:p, ss:'E' });
      }
    } else {
      // 3.8 Å step chain with dihedral-inspired angular deflection
      const target = center.clone().add(new THREE.Vector3((rng()-0.5)*7,(rng()-0.5)*7,(rng()-0.5)*7));
      let loopCur = cur.clone();
      let loopDir = target.clone().sub(cur);
      if (loopDir.length() < 0.01) loopDir.set(rng()-0.5, rng()-0.5, rng()-0.5);
      loopDir.normalize();
      for (let i = 0; i < el.n; i++) {
        const toTarget = target.clone().sub(loopCur);
        if (toTarget.length() > 0.01) loopDir.lerp(toTarget.normalize(), 0.3).normalize();
        const phi = (rng() - 0.5) * Math.PI * 0.6;
        const psi = (rng() - 0.5) * Math.PI * 0.4;
        let perpSeed = new THREE.Vector3(1,0,0);
        if (Math.abs(loopDir.dot(perpSeed)) > 0.85) perpSeed.set(0,1,0);
        const perp1 = new THREE.Vector3().crossVectors(loopDir, perpSeed).normalize();
        const perp2 = new THREE.Vector3().crossVectors(loopDir, perp1).normalize();
        const step = loopDir.clone()
          .add(perp1.clone().multiplyScalar(Math.sin(phi) * 0.5))
          .add(perp2.clone().multiplyScalar(Math.sin(psi) * 0.3))
          .normalize().multiplyScalar(3.8 * SCALE);
        loopCur = loopCur.clone().add(step);
        points.push({ pos: loopCur.clone(), ss:'L' });
        loopDir = step.clone().normalize();
      }
    }
    cur = points[points.length-1].pos.clone();
    const toC = center.clone().sub(cur).normalize();
    dir = dir.clone().lerp(toC, 0.35+rng()*0.3).normalize();
    dir.add(new THREE.Vector3((rng()-0.5)*0.5,(rng()-0.5)*0.5,(rng()-0.5)*0.5)).normalize();
  }
  return points;
}

function generateLinker(fromPos, toPos, n) {
  const pts = [];
  const rng = makeRng(Math.floor(fromPos.x*100+toPos.y*100));
  for (let i = 0; i <= n; i++) {
    const t = i/n;
    const p = fromPos.clone().lerp(toPos,t);
    const w = 3*Math.sin(t*Math.PI);
    p.add(new THREE.Vector3((rng()-0.5)*w*2,(rng()-0.5)*w*2,(rng()-0.5)*w*1.5));
    pts.push({ pos:p, ss:'L' });
  }
  return pts;
}


// ==================== VIEWER CLASS ====================

export class POLEViewer {
  constructor(containerEl, tooltipEl) {
    this.container = containerEl;
    this.tooltipEl = tooltipEl;
    this.autoRotate = true;
    this.isDragging = false;
    this.prevMouse = { x:0, y:0 };
    this.cameraTheta = -0.3;
    this.cameraPhi = 0.15;
    this.cameraRadius = 95;
    this.hoverTarget = null;
    this.groups = {};
    this.domainMeshes = {};
    this.mutMarkers = {};
    this.activeSitePos = new THREE.Vector3();
    this.cameraTarget = new THREE.Vector3(0, 0, 0);

    // New state for features
    this.backboneCache = {};
    this.colorMode = 'default';
    this.colorScheme = 'default';
    this.measureMode = false;
    this.measurePoints = [];
    this.measureLine = null;
    this.measureLabel = null;
    this.exoPolT = 0;
    this.exoPolPlaying = false;
    this.fpsFrames = 0;
    this.fpsTime = 0;
    this.fpsCallback = null;

    this._init();
  }

  _init() {
    this._createScene();
    this._createLighting();
    this._buildProtein();
    this._buildDNA();
    this._buildActiveSite();
    this._buildZincMotifs();
    this._buildFeSCluster();
    this._buildExoSite();
    this._buildAccessorySubunits();
    this._buildPCNA();
    this._buildMutationMarkers();
    this._buildAnnotations();
    this._buildHydrogenBonds();
    this._buildBetaSheetHBonds();
    this._buildExoPolPath();
    this._buildParticles();
    this._setupControls();
    this._setupTouch();
    this._startRenderLoop();
  }

  _createScene() {
    this.scene = new THREE.Scene();
    this.scene.fog = new THREE.FogExp2(PAL.bg, 0.003);
    this.camera = new THREE.PerspectiveCamera(40,1,0.1,2000);
    this.renderer = new THREE.WebGLRenderer({ antialias:true, alpha:true, preserveDrawingBuffer:true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio,2));
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
    this.renderer.toneMappingExposure = 1.4;
    this.container.appendChild(this.renderer.domElement);
    const resize = () => {
      const w = this.container.clientWidth, h = this.container.clientHeight;
      this.renderer.setSize(w,h,false);
      this.camera.aspect = w/h;
      this.camera.updateProjectionMatrix();
    };
    new ResizeObserver(resize).observe(this.container);
    resize();
    this._updateCamera();
  }

  _createLighting() {
    this.scene.add(new THREE.AmbientLight(0x8090a0, 0.8));
    const key = new THREE.DirectionalLight(0xeee8dd, 1.3);
    key.position.set(18,25,20);
    this.scene.add(key);
    const fill = new THREE.DirectionalLight(0x6a7a8a, 0.4);
    fill.position.set(-10,-15,10);
    this.scene.add(fill);
    const rim = new THREE.DirectionalLight(0x7abba8, 0.5);
    rim.position.set(-20,8,-15);
    this.scene.add(rim);
    this.activeSiteLight = new THREE.PointLight(0xc45c8c, 0.6, 18);
    this.scene.add(this.activeSiteLight);
  }

  // ==================== PROTEIN: cartoon ribbon representation ====================

  _buildProtein() {
    this.groups.protein = new THREE.Group();
    this.groups.surface = new THREE.Group();
    this.groups.surface.visible = true;
    this.scene.add(this.groups.protein);
    this.groups.protein.add(this.groups.surface);

    const backbones = {};

    for (const key of DOMAIN_ORDER) {
      const dom = DOMAINS[key];
      const bb = generateDomainBackbone(key);
      backbones[key] = bb;
      this.backboneCache[key] = bb;
      const domColor = PAL[dom.color];

      // Build ribbon meshes per SS segment
      const segments = segmentBySSType(bb);
      const ribbons = [];
      for (let s = 0; s < segments.length; s++) {
        const seg = segments[s];
        // 1-point overlap at boundaries for smooth joins
        const lo = Math.max(0, seg.startIdx - (s > 0 ? 1 : 0));
        const hi = Math.min(bb.length - 1, seg.endIdx + (s < segments.length - 1 ? 1 : 0));
        const pts = [];
        for (let i = lo; i <= hi; i++) pts.push(bb[i].pos);
        if (pts.length >= 2) {
          const mesh = buildSSElement(pts, seg.type, domColor);
          mesh.userData.ssType = seg.type;
          mesh.userData.domainKey = key;
          mesh.userData.startIdx = seg.startIdx;
          mesh.userData.endIdx = seg.endIdx;
          this.groups.protein.add(mesh);
          ribbons.push(mesh);
        }
      }

      // Raycasting proxy
      const bbox = new THREE.Box3();
      for (const r of bb) bbox.expandByPoint(r.pos);
      const center = new THREE.Vector3();
      bbox.getCenter(center);
      const size = new THREE.Vector3();
      bbox.getSize(size);
      const proxy = new THREE.Mesh(
        new THREE.SphereGeometry(Math.max(size.x,size.y,size.z)*0.5, 8, 8),
        new THREE.MeshBasicMaterial({ visible:false })
      );
      proxy.position.copy(center);
      proxy.userData = { domain:key, name:dom.name, detail:`res ${dom.range[0]}\u2013${dom.range[1]}` };
      this.groups.protein.add(proxy);
      this.domainMeshes[key] = { proxy, ribbons };
    }

    // Inter-domain linkers as thin tubes
    for (let i = 0; i < DOMAIN_ORDER.length - 1; i++) {
      const kA = DOMAIN_ORDER[i], kB = DOMAIN_ORDER[i+1];
      const bbA = backbones[kA], bbB = backbones[kB];
      const linkerN = Math.max(8, Math.floor((DOMAINS[kB].range[0] - DOMAINS[kA].range[1])/5));
      const linkerBB = generateLinker(bbA[bbA.length-1].pos, bbB[0].pos, linkerN);
      const linkerPts = linkerBB.map(p => p.pos);
      if (linkerPts.length >= 2) {
        const mesh = buildSSElement(linkerPts, 'L', PAL.linker);
        this.groups.protein.add(mesh);
      }
    }

    // Surface envelope (toggle)
    for (const key of DOMAIN_ORDER) {
      const c = DOMAINS[key].center;
      const col = PAL[DOMAINS[key].color];
      const mesh = new THREE.Mesh(
        new THREE.IcosahedronGeometry(6,3),
        new THREE.MeshPhysicalMaterial({
          color:col, transparent:true, opacity:0.09, roughness:0.9,
          side:THREE.DoubleSide, depthWrite:false })
      );
      mesh.position.set(c[0],c[1],c[2]);
      mesh.userData.surfaceDomain = key;
      this.groups.surface.add(mesh);
    }
  }

  // ==================== DNA: backbone tubes + base-pair rungs ====================

  _buildDNA() {
    this.groups.dna = new THREE.Group();
    this.scene.add(this.groups.dna);

    const axisPts = [];
    for (let i = 0; i < 80; i++) {
      const t = i/79;
      axisPts.push(new THREE.Vector3(
        -8+t*26, -9+Math.sin(t*Math.PI)*5+t*4,
        8-t*10+Math.sin(t*Math.PI*1.3)*2.5));
    }
    const axisCurve = new THREE.CatmullRomCurve3(axisPts);

    const helixR = 10.0*SCALE*0.5;
    const turnsPerUnit = 6;
    const segments = 400;
    const primerCutoff = 0.52;
    const grooveOffset = Math.PI*(154/180);
    const strand1=[], strand2=[];
    const basePairs = [];

    for (let i = 0; i <= segments; i++) {
      const t = i/segments;
      const center = axisCurve.getPoint(t);
      const tangent = axisCurve.getTangent(t).normalize();
      let up = Math.abs(tangent.y)<0.95 ? new THREE.Vector3(0,1,0) : new THREE.Vector3(1,0,0);
      const normal = new THREE.Vector3().crossVectors(tangent,up).normalize();
      const binormal = new THREE.Vector3().crossVectors(tangent,normal).normalize();
      const angle = t*turnsPerUnit*Math.PI*2;

      const o1 = normal.clone().multiplyScalar(Math.cos(angle)*helixR)
        .add(binormal.clone().multiplyScalar(Math.sin(angle)*helixR));
      const o2 = normal.clone().multiplyScalar(Math.cos(angle+grooveOffset)*helixR)
        .add(binormal.clone().multiplyScalar(Math.sin(angle+grooveOffset)*helixR));

      strand1.push(center.clone().add(o1));
      strand2.push(center.clone().add(o2));

      if (i % 14 === 0 && t < primerCutoff) {
        basePairs.push({ p1: center.clone().add(o1), p2: center.clone().add(o2) });
      }
    }

    // Template backbone tube
    if (strand1.length > 2) {
      const curve = new THREE.CatmullRomCurve3(strand1);
      const geo = new THREE.TubeGeometry(curve, 200, 0.12, 8, false);
      const mat = new THREE.MeshPhysicalMaterial({
        color: PAL.dnaTemplate, roughness:0.5, metalness:0.05, clearcoat:0.15
      });
      this.groups.dna.add(new THREE.Mesh(geo,mat));
    }
    // Primer backbone tube
    const primerPts = strand2.slice(0, Math.floor(strand2.length*primerCutoff));
    if (primerPts.length > 2) {
      const curve = new THREE.CatmullRomCurve3(primerPts);
      const geo = new THREE.TubeGeometry(curve, 120, 0.12, 8, false);
      const mat = new THREE.MeshPhysicalMaterial({
        color: PAL.dnaPrimer, roughness:0.5, metalness:0.05, clearcoat:0.15
      });
      this.groups.dna.add(new THREE.Mesh(geo,mat));
    }

    // Base-pair rungs (cylinders)
    for (const bp of basePairs) {
      const dir = bp.p2.clone().sub(bp.p1);
      const len = dir.length();
      const rung = new THREE.Mesh(
        new THREE.CylinderGeometry(0.04, 0.04, len, 6),
        new THREE.MeshPhysicalMaterial({ color: PAL.dnaBase, roughness:0.5, metalness:0.05 })
      );
      rung.position.copy(bp.p1.clone().add(bp.p2).multiplyScalar(0.5));
      rung.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0), dir.normalize());
      this.groups.dna.add(rung);
    }

    const pidx = Math.floor(segments*primerCutoff);
    this.activeSitePos.copy(strand2[pidx] || new THREE.Vector3());
  }

  // ==================== ACTIVE SITE: ball-and-stick + annotation ====================

  _buildActiveSite() {
    this.groups.metals = new THREE.Group();
    this.groups.dntp = new THREE.Group();
    this.scene.add(this.groups.metals);
    this.scene.add(this.groups.dntp);

    const asp = this.activeSitePos;
    const mgApos = asp.clone().add(new THREE.Vector3(-0.5,0.25,0.1));
    const mgBpos = asp.clone().add(new THREE.Vector3(0.7,-0.15,0.35));

    for (const [pos, label] of [[mgApos,'A'],[mgBpos,'B']]) {
      const core = new THREE.Mesh(
        new THREE.SphereGeometry(0.35,16,16),
        new THREE.MeshPhysicalMaterial({
          color:0x60a830, roughness:0.2, metalness:0.7,
          emissive:0x60a830, emissiveIntensity:0.5 })
      );
      core.position.copy(pos);
      core.userData = { name:`Mg\u00b2\u207a (Metal ${label})`,
        detail: label==='A' ? 'Activates 3\u2032-OH \u00b7 D640/D642' : 'Coordinates \u03b2,\u03b3-phosphates \u00b7 D860' };
      this.groups.metals.add(core);
    }

    // Octahedral coordination: 6 ligand spheres per metal
    const octDirs = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]];
    const ligDist = 2.1 * SCALE;
    for (const pos of [mgApos,mgBpos]) {
      for (const d of octDirs) {
        const lPos = pos.clone().add(new THREE.Vector3(...d).multiplyScalar(ligDist));
        const lig = this._atom(0xc03030, 0.08);
        lig.position.copy(lPos);
        this.groups.metals.add(lig);
        this.groups.metals.add(this._dash(pos,lPos,0x60a830));
      }
    }
    this.groups.metals.add(this._dash(mgApos,mgBpos,0x60a830));

    // Mg-Mg distance label
    const mgMgMid = mgApos.clone().add(mgBpos).multiplyScalar(0.5);
    const mgMgDist = mgApos.distanceTo(mgBpos) / SCALE;
    const mgMgLabelPos = mgMgMid.clone().add(new THREE.Vector3(6, -8, 6));
    const mgMgLabel = this._spriteLabel(`Mg\u2013Mg: ${mgMgDist.toFixed(1)} \u00c5`, '#60a830', '#60a830', 220);
    mgMgLabel.scale.set(7, 2, 1);
    mgMgLabel.position.copy(mgMgLabelPos);
    mgMgLabel.userData.labelPriority = 1;
    this.groups.metals.add(mgMgLabel);
    this.groups.metals.add(this._leaderLine(mgMgLabelPos, mgMgMid, PAL.magnesium));

    const catRes = [
      { l:'D640', pos:asp.clone().add(new THREE.Vector3(-1.2,0.8,-0.5)), mg:mgApos },
      { l:'D642', pos:asp.clone().add(new THREE.Vector3(-0.8,-0.6,-0.8)), mg:mgApos },
      { l:'D860', pos:asp.clone().add(new THREE.Vector3(1.2,0.6,-0.4)), mg:mgBpos },
    ];
    for (const res of catRes) {
      const dir = res.pos.clone().sub(asp).normalize();
      const end = res.pos.clone().add(dir.multiplyScalar(0.8));
      this.groups.metals.add(this._stick(res.pos,end,0x707070,0.05));
      const perp = new THREE.Vector3(-dir.y,dir.x,dir.z*0.5).normalize();
      for (const s of [1,-1]) {
        const o = end.clone().add(perp.clone().multiplyScalar(s*0.35));
        this.groups.metals.add(this._atom(0xc03030,0.12));
        this.groups.metals.children[this.groups.metals.children.length-1].position.copy(o);
        this.groups.metals.add(this._stick(end,o,0xc03030,0.03));
        this.groups.metals.add(this._dash(o,res.mg,0x888888));
      }
      // Asp-Mg distance label -pushed away from center
      const distVal = res.pos.distanceTo(res.mg) / SCALE;
      const distMid = res.pos.clone().add(res.mg).multiplyScalar(0.5);
      const outDir = res.pos.clone().sub(asp).normalize();
      const dLabelPos = distMid.clone().add(outDir.multiplyScalar(6)).add(new THREE.Vector3(0, 4, 0));
      const dLabel = this._spriteLabel(`${res.l}: ${distVal.toFixed(1)} \u00c5`, '#aaaaaa', '#888888', 200);
      dLabel.scale.set(6, 1.8, 1);
      dLabel.position.copy(dLabelPos);
      dLabel.userData.labelPriority = 0;
      this.groups.metals.add(dLabel);
      this.groups.metals.add(this._leaderLine(dLabelPos, distMid, new THREE.Color(0x888888)));
    }

    const dBase = asp.clone().add(new THREE.Vector3(1.8,0.8,-0.6));
    const pA=dBase.clone(), pB=dBase.clone().add(new THREE.Vector3(0.65,0.15,0)),
          pG=dBase.clone().add(new THREE.Vector3(1.3,0.3,0.05));
    for (const pp of [pA,pB,pG]) {
      this.groups.dntp.add(this._atom(0xd07020,0.18));
      this.groups.dntp.children[this.groups.dntp.children.length-1].position.copy(pp);
    }
    for (const [a,b] of [[pA,pB],[pB,pG]]) {
      const mid=a.clone().add(b).multiplyScalar(0.5);
      this.groups.dntp.add(this._atom(0xc03030,0.09));
      this.groups.dntp.children[this.groups.dntp.children.length-1].position.copy(mid);
      this.groups.dntp.add(this._stick(a,mid,0xd07020,0.03));
      this.groups.dntp.add(this._stick(mid,b,0xd07020,0.03));
    }
    const sugarC = dBase.clone().add(new THREE.Vector3(-0.7,0,0));
    this.groups.dntp.add(this._atom(0x707070,0.12));
    this.groups.dntp.children[this.groups.dntp.children.length-1].position.copy(sugarC);
    this.groups.dntp.add(this._stick(sugarC,pA,0xd07020,0.03));
    const baseC = sugarC.clone().add(new THREE.Vector3(0,0.6,0));
    this.groups.dntp.add(this._atom(0x3050c0,0.1));
    this.groups.dntp.children[this.groups.dntp.children.length-1].position.copy(baseC);
    this.groups.dntp.add(this._stick(sugarC,baseC,0x3050c0,0.025));

    this.groups.dntp.add(this._dash(pB,mgBpos,0x60a830));
    this.groups.dntp.add(this._dash(pG,mgBpos,0x60a830));

    // Template nucleotide opposite the dNTP (Watson-Crick base pairing)
    const templateBase = baseC.clone().add(new THREE.Vector3(0, 1.4, 0));
    const templateAtom = this._atom(0x3050c0, 0.1);
    templateAtom.position.copy(templateBase);
    this.groups.dntp.add(templateAtom);
    const templateRing = templateBase.clone().add(new THREE.Vector3(0.3, 0.2, 0));
    const templateRingAtom = this._atom(0x3050c0, 0.08);
    templateRingAtom.position.copy(templateRing);
    this.groups.dntp.add(templateRingAtom);
    this.groups.dntp.add(this._stick(templateBase, templateRing, 0x3050c0, 0.025));

    // 3 Watson-Crick H-bond dashed lines
    const hbondColor = 0x88aacc;
    const hbondOffsets = [
      [0, 0, 0], [0.15, 0.05, 0.1], [-0.15, 0.05, -0.1]
    ];
    for (const off of hbondOffsets) {
      const from = baseC.clone().add(new THREE.Vector3(off[0], 0.1 + off[1], off[2]));
      const to = templateBase.clone().add(new THREE.Vector3(off[0], -0.1 + off[1], off[2]));
      const geo = new THREE.BufferGeometry().setFromPoints([from, to]);
      const mat = new THREE.LineDashedMaterial({
        color: hbondColor, dashSize: 0.08, gapSize: 0.06,
        transparent: true, opacity: 0.5
      });
      const line = new THREE.Line(geo, mat);
      line.computeLineDistances();
      this.groups.dntp.add(line);
    }

    // W-C base pair label -pushed out from model center
    const wcAnchor = baseC.clone().add(new THREE.Vector3(0, 0.7, 0));
    const wcLabelPos = wcAnchor.clone().add(new THREE.Vector3(8, 6, -6));
    const wcLabel = this._spriteLabel('W\u2013C base pair', '#88aacc', '#88aacc', 220);
    wcLabel.scale.set(7, 2, 1);
    wcLabel.position.copy(wcLabelPos);
    wcLabel.userData.labelPriority = 1;
    this.groups.dntp.add(wcLabel);
    this.groups.dntp.add(this._leaderLine(wcLabelPos, wcAnchor, new THREE.Color(0x88aacc)));

    this.groups.dntp.userData = { name:'Incoming dNTP',
      detail:'Templated nucleotide \u00b7 \u03b2,\u03b3-phosphates \u2192 Metal B' };
    this.activeSiteLight.position.copy(asp);
  }

  _atom(color,r) {
    return new THREE.Mesh(new THREE.SphereGeometry(r,12,12),
      new THREE.MeshPhysicalMaterial({ color, roughness:0.25, metalness:0.5,
        emissive:color, emissiveIntensity:0.25 }));
  }
  _stick(a,b,color,r) {
    const d=b.clone().sub(a);
    const m=new THREE.Mesh(new THREE.CylinderGeometry(r,r,d.length(),6),
      new THREE.MeshPhysicalMaterial({ color, roughness:0.4, metalness:0.3 }));
    m.position.copy(a.clone().add(b).multiplyScalar(0.5));
    m.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0),d.normalize());
    return m;
  }
  _dash(a,b,color) {
    const g=new THREE.BufferGeometry().setFromPoints([a,b]);
    const m=new THREE.LineDashedMaterial({ color, dashSize:0.2, gapSize:0.15, transparent:true, opacity:0.4 });
    const l=new THREE.Line(g,m); l.computeLineDistances(); return l;
  }

  // ==================== ZINC-BINDING MOTIFS (CTD) ====================

  _buildZincMotifs() {
    this.groups.zincMotifs = new THREE.Group();
    this.scene.add(this.groups.zincMotifs);

    const ctdCenter = new THREE.Vector3(...DOMAINS.ctd.center);
    const sites = [
      { name: 'CysA Zn-finger', offset: [-2.5, 1, -1.5] },
      { name: 'CysB Zn-finger', offset: [2, -1.5, 1] },
    ];
    const znPositions = [];

    for (const site of sites) {
      const znPos = ctdCenter.clone().add(new THREE.Vector3(...site.offset));
      znPositions.push(znPos);

      // Zn²⁺ ion
      const znSphere = new THREE.Mesh(
        new THREE.SphereGeometry(0.3, 16, 16),
        new THREE.MeshPhysicalMaterial({
          color: 0x7799aa, roughness: 0.2, metalness: 0.7,
          emissive: 0x7799aa, emissiveIntensity: 0.4
        })
      );
      znSphere.position.copy(znPos);
      znSphere.userData = { name: `Zn\u00b2\u207a (${site.name})`, detail: 'Tetrahedral Cys\u2084 coordination' };
      this.groups.zincMotifs.add(znSphere);

      // 4 coordinating cysteine sulfurs in tetrahedral geometry
      const tetDirs = [
        [1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]
      ];
      const bondLen = 2.3 * SCALE;
      for (const td of tetDirs) {
        const sDir = new THREE.Vector3(...td).normalize();
        const sPos = znPos.clone().add(sDir.clone().multiplyScalar(bondLen));

        // Yellow sulfur atom
        const sulfur = this._atom(0xcccc00, 0.12);
        sulfur.position.copy(sPos);
        this.groups.zincMotifs.add(sulfur);

        // Zn-S stick
        this.groups.zincMotifs.add(this._stick(znPos, sPos, 0x99aa88, 0.04));

        // C-alpha stub
        const caPos = sPos.clone().add(sDir.clone().multiplyScalar(0.5));
        const ca = this._atom(0x707070, 0.08);
        ca.position.copy(caPos);
        this.groups.zincMotifs.add(ca);
        this.groups.zincMotifs.add(this._stick(sPos, caPos, 0x707070, 0.025));
      }

      // Label -pushed away from model
      const labelPos = znPos.clone().add(new THREE.Vector3(
        site.offset[0] > 0 ? 4 : -4, 6, site.offset[2] > 0 ? 3 : -3));
      const label = this._spriteLabel(site.name, '#d0d8e0', '#7799aa', 260);
      label.scale.set(9, 2.5, 1);
      label.position.copy(labelPos);
      label.userData.labelPriority = 2;
      this.groups.zincMotifs.add(label);
      this.groups.zincMotifs.add(this._leaderLine(labelPos, znPos, new THREE.Color(0x7799aa)));
    }

    // Bridging dash between the two Zn sites
    if (znPositions.length === 2) {
      this.groups.zincMotifs.add(this._dash(znPositions[0], znPositions[1], 0x7799aa));
    }
  }

  // ==================== [4Fe-4S] CLUSTER (P DOMAIN) ====================

  _buildFeSCluster() {
    this.groups.fesCluster = new THREE.Group();
    this.scene.add(this.groups.fesCluster);

    const pCenter = new THREE.Vector3(...DOMAINS.pdomain.center);
    const clusterPos = pCenter.clone().add(new THREE.Vector3(0.5, 0.8, -0.5));

    // Cubane [4Fe-4S] core: alternating Fe and S at cube vertices
    const cubeR = 0.55 * SCALE;
    const fePositions = [
      clusterPos.clone().add(new THREE.Vector3( cubeR,  cubeR,  cubeR)),
      clusterPos.clone().add(new THREE.Vector3(-cubeR, -cubeR,  cubeR)),
      clusterPos.clone().add(new THREE.Vector3(-cubeR,  cubeR, -cubeR)),
      clusterPos.clone().add(new THREE.Vector3( cubeR, -cubeR, -cubeR)),
    ];
    const sPositions = [
      clusterPos.clone().add(new THREE.Vector3(-cubeR,  cubeR,  cubeR)),
      clusterPos.clone().add(new THREE.Vector3( cubeR, -cubeR,  cubeR)),
      clusterPos.clone().add(new THREE.Vector3( cubeR,  cubeR, -cubeR)),
      clusterPos.clone().add(new THREE.Vector3(-cubeR, -cubeR, -cubeR)),
    ];

    // Fe atoms (rust-orange)
    for (let i = 0; i < 4; i++) {
      const fe = new THREE.Mesh(
        new THREE.SphereGeometry(0.2, 12, 12),
        new THREE.MeshPhysicalMaterial({
          color: 0xcc6633, roughness: 0.3, metalness: 0.8,
          emissive: 0xcc6633, emissiveIntensity: 0.4
        })
      );
      fe.position.copy(fePositions[i]);
      fe.userData = { name: `Fe${i+1} ([4Fe-4S])`, detail: 'Iron-sulfur cluster / CysX motif' };
      this.groups.fesCluster.add(fe);
    }

    // Inorganic S atoms (yellow)
    for (let i = 0; i < 4; i++) {
      const s = new THREE.Mesh(
        new THREE.SphereGeometry(0.18, 12, 12),
        new THREE.MeshPhysicalMaterial({
          color: 0xcccc33, roughness: 0.4, metalness: 0.5,
          emissive: 0xcccc33, emissiveIntensity: 0.3
        })
      );
      s.position.copy(sPositions[i]);
      this.groups.fesCluster.add(s);
    }

    // Fe-S bonds within cubane
    for (const feP of fePositions) {
      for (const sP of sPositions) {
        if (feP.distanceTo(sP) < cubeR * 2.2) {
          this.groups.fesCluster.add(this._stick(feP, sP, 0x888866, 0.03));
        }
      }
    }

    // 4 coordinating cysteine sulfurs (CysX motif)
    const cysDirs = [[1.2,0.5,0],[-0.3,1.2,0.8],[0.5,-0.3,1.2],[-1.0,-0.8,-0.3]];
    for (let i = 0; i < 4; i++) {
      const cysS = fePositions[i].clone().add(
        new THREE.Vector3(cysDirs[i][0], cysDirs[i][1], cysDirs[i][2]).normalize().multiplyScalar(1.2 * SCALE)
      );
      const atom = this._atom(0xcccc33, 0.12);
      atom.position.copy(cysS);
      this.groups.fesCluster.add(atom);
      this.groups.fesCluster.add(this._stick(fePositions[i], cysS, 0xcccc33, 0.025));
    }

    // Label
    const labelPos = clusterPos.clone().add(new THREE.Vector3(6, 14, -6));
    const label = this._spriteLabel('[4Fe-4S] cluster', '#d0d8e0', '#cc6633', 240);
    label.scale.set(10, 2.8, 1);
    label.position.copy(labelPos);
    label.userData.labelPriority = 2;
    this.groups.fesCluster.add(label);
    this.groups.fesCluster.add(this._leaderLine(labelPos, clusterPos, new THREE.Color(0xcc6633)));
  }

  // ==================== EXONUCLEASE ACTIVE SITE ====================

  _buildExoSite() {
    this.groups.exoSite = new THREE.Group();
    this.scene.add(this.groups.exoSite);

    const exoCenter = new THREE.Vector3(...DOMAINS.exo.center);
    const mgOffsets = [[0.3, -0.5, 1.0], [-0.8, -0.3, 1.5]];
    const mgPositions = mgOffsets.map(o => exoCenter.clone().add(new THREE.Vector3(...o)));
    this.exoSitePos = exoCenter.clone();

    // Two Mg²⁺ ions with octahedral coordination
    for (let m = 0; m < 2; m++) {
      const mgPos = mgPositions[m];
      const mgSphere = new THREE.Mesh(
        new THREE.SphereGeometry(0.35, 16, 16),
        new THREE.MeshPhysicalMaterial({
          color: 0x60a830, roughness: 0.2, metalness: 0.7,
          emissive: 0x60a830, emissiveIntensity: 0.5
        })
      );
      mgSphere.position.copy(mgPos);
      mgSphere.userData = { name: `Mg\u00b2\u207a (EXO Metal ${m === 0 ? 'A' : 'B'})`,
        detail: 'Exonuclease two-metal-ion center' };
      this.groups.exoSite.add(mgSphere);

      // 6 octahedral ligand spheres
      const octDirs = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]];
      const ligDist = 2.1 * SCALE;
      for (const d of octDirs) {
        const lPos = mgPos.clone().add(new THREE.Vector3(...d).multiplyScalar(ligDist));
        const lig = this._atom(0xc03030, 0.08);
        lig.position.copy(lPos);
        this.groups.exoSite.add(lig);
        this.groups.exoSite.add(this._dash(mgPos, lPos, 0x60a830));
      }
    }

    // Mg-Mg bridging dash
    this.groups.exoSite.add(this._dash(mgPositions[0], mgPositions[1], 0x60a830));

    // Catalytic residues: D275, E277, D368, D370
    const catRes = [
      { l: 'D275', offset: [-1.5, 0.8, 0.5] },
      { l: 'E277', offset: [-1.2, -0.6, 1.8] },
      { l: 'D368', offset: [0.8, 0.5, 2.0] },
      { l: 'D370', offset: [1.2, -0.8, 0.8] },
    ];
    for (const res of catRes) {
      const resPos = exoCenter.clone().add(new THREE.Vector3(...res.offset));
      // Ball-and-stick
      const ca = this._atom(0x707070, 0.12);
      ca.position.copy(resPos);
      this.groups.exoSite.add(ca);

      const dir = resPos.clone().sub(exoCenter).normalize();
      const end = resPos.clone().add(dir.clone().multiplyScalar(0.6));
      this.groups.exoSite.add(this._stick(resPos, end, 0x707070, 0.04));

      // Carboxylate oxygens
      const perp = new THREE.Vector3(-dir.y, dir.x, dir.z * 0.5).normalize();
      for (const s of [1, -1]) {
        const oPos = end.clone().add(perp.clone().multiplyScalar(s * 0.3));
        const oxy = this._atom(0xc03030, 0.1);
        oxy.position.copy(oPos);
        this.groups.exoSite.add(oxy);
        this.groups.exoSite.add(this._stick(end, oPos, 0xc03030, 0.025));
        // Coordination dash to nearest metal
        const nearMg = mgPositions[0].distanceTo(oPos) < mgPositions[1].distanceTo(oPos)
          ? mgPositions[0] : mgPositions[1];
        this.groups.exoSite.add(this._dash(oPos, nearMg, 0x888888));
      }
    }

    // Label -pushed away from model
    const exoLabelPos = exoCenter.clone().add(new THREE.Vector3(-8, 16, 8));
    const label = this._spriteLabel('EXO SITE \u00b7 2\u00d7 Mg\u00b2\u207a', '#d0d8e0', '#b85c5c', 300);
    label.scale.set(11, 3, 1);
    label.position.copy(exoLabelPos);
    label.userData.labelPriority = 2;
    this.groups.exoSite.add(label);
    this.groups.exoSite.add(this._leaderLine(exoLabelPos, exoCenter, PAL.exo));
  }

  // ==================== ACCESSORY SUBUNITS ====================

  _buildAccessorySubunits() {
    this.groups.accessory = new THREE.Group();
    this.groups.accessory.visible = true;
    this.scene.add(this.groups.accessory);

    const ctdCenter = new THREE.Vector3(...DOMAINS.ctd.center);
    const subunits = [
      { name: 'POLE2 (p59)', offset: [8, 2, -3], radius: 5.0, color: 0x6699aa },
      { name: 'POLE3 (p17)', offset: [10, -3, 2], radius: 2.8, color: 0x99aa66 },
      { name: 'POLE4 (p12)', offset: [10, -3, -3], radius: 2.2, color: 0xaa8866 },
    ];

    for (const sub of subunits) {
      const pos = ctdCenter.clone().add(new THREE.Vector3(...sub.offset));

      // Ghost surface
      const ghost = new THREE.Mesh(
        new THREE.IcosahedronGeometry(sub.radius, 3),
        new THREE.MeshPhysicalMaterial({
          color: sub.color, transparent: true, opacity: 0.06,
          roughness: 0.9, side: THREE.DoubleSide, depthWrite: false
        })
      );
      ghost.position.copy(pos);
      ghost.userData = { name: sub.name, detail: 'Accessory subunit' };
      this.groups.accessory.add(ghost);

      // Wireframe outline
      const wire = new THREE.Mesh(
        new THREE.IcosahedronGeometry(sub.radius, 2),
        new THREE.MeshBasicMaterial({
          color: sub.color, wireframe: true, transparent: true, opacity: 0.08
        })
      );
      wire.position.copy(pos);
      this.groups.accessory.add(wire);

      // Label -pushed away from model
      const subLabelPos = pos.clone().add(new THREE.Vector3(
        sub.offset[0] > 9 ? 4 : -3, sub.radius + 5, sub.offset[2] > 0 ? 3 : -3));
      const label = this._spriteLabel(sub.name, '#d0d8e0', '#' + new THREE.Color(sub.color).getHexString(), 240);
      label.scale.set(9, 2.5, 1);
      label.position.copy(subLabelPos);
      label.userData.labelPriority = 2;
      this.groups.accessory.add(label);
      this.groups.accessory.add(this._leaderLine(subLabelPos, pos, new THREE.Color(sub.color)));

      // Dashed line from CTD center to subunit
      this.groups.accessory.add(this._dash(ctdCenter, pos, sub.color));
    }
  }

  // ==================== PCNA SLIDING CLAMP ====================

  _buildPCNA() {
    this.groups.pcna = new THREE.Group();
    this.groups.pcna.visible = true;
    this.scene.add(this.groups.pcna);

    const thumbCenter = new THREE.Vector3(...DOMAINS.thumb.center);
    const pcnaPos = thumbCenter.clone().add(new THREE.Vector3(5, -6, 10));

    // Torus geometry for PCNA ring
    const torusGeo = new THREE.TorusGeometry(4.0, 1.2, 16, 48);
    const ghostMat = new THREE.MeshPhysicalMaterial({
      color: 0x7abba8, transparent: true, opacity: 0.12,
      emissive: 0x7abba8, emissiveIntensity: 0.15,
      roughness: 0.8, side: THREE.DoubleSide, depthWrite: false
    });
    const ghostTorus = new THREE.Mesh(torusGeo, ghostMat);
    ghostTorus.position.copy(pcnaPos);
    ghostTorus.rotation.y = Math.PI / 2;
    ghostTorus.rotation.x = 0.2;
    ghostTorus.userData = { name: 'PCNA sliding clamp', detail: 'Homotrimeric processivity factor · Pol \u03b5 contacts all three protomers via PIP-box (Q1180), Thumb insertion (res. 1102\u20131122), and P domain \u2014 tripartite interface impervious to FEN1/Pol \u03b4 displacement (Roske & Yeeles 2024)' };
    this.groups.pcna.add(ghostTorus);

    // Wireframe torus
    const wireTorus = new THREE.Mesh(
      new THREE.TorusGeometry(4.0, 1.2, 8, 24),
      new THREE.MeshBasicMaterial({
        color: 0x7abba8, wireframe: true, transparent: true, opacity: 0.1
      })
    );
    wireTorus.position.copy(pcnaPos);
    wireTorus.rotation.y = Math.PI / 2;
    wireTorus.rotation.x = 0.2;
    this.groups.pcna.add(wireTorus);

    // 3 symmetry markers (PCNA homotrimer)
    for (let i = 0; i < 3; i++) {
      const ang = (i / 3) * Math.PI * 2;
      const markerPos = pcnaPos.clone().add(new THREE.Vector3(
        Math.cos(ang) * 4.0 * Math.cos(0.2),
        Math.sin(0.2) * Math.cos(ang) * 4.0,
        Math.sin(ang) * 4.0
      ));
      const marker = new THREE.Mesh(
        new THREE.SphereGeometry(0.25, 8, 8),
        new THREE.MeshPhysicalMaterial({
          color: 0x7abba8, emissive: 0x7abba8, emissiveIntensity: 0.5,
          roughness: 0.3
        })
      );
      marker.position.copy(markerPos);
      this.groups.pcna.add(marker);
    }

    // Tripartite interface labels
    const triLabels = [
      { text:'PIP-box (Q1180)', offset:[-7, -5, 3] },
      { text:'Thumb insert (1102\u20131122)', offset:[5, -4, -5] },
      { text:'P-domain contact', offset:[-5, 4, -6] },
    ];
    for (const tl of triLabels) {
      const tlPos = pcnaPos.clone().add(new THREE.Vector3(...tl.offset));
      const tlLabel = this._spriteLabel(tl.text, '#7abba8', '#7abba8', 240);
      tlLabel.scale.set(9, 2.5, 1);
      tlLabel.position.copy(tlPos);
      tlLabel.userData.labelPriority = 1;
      this.groups.pcna.add(tlLabel);
      this.groups.pcna.add(this._leaderLine(tlPos, pcnaPos, new THREE.Color(0x7abba8)));
    }

    // Main label
    const mainLabelPos = pcnaPos.clone().add(new THREE.Vector3(0, 10, 0));
    const mainLabel = this._spriteLabel('PCNA sliding clamp', '#d0d8e0', '#7abba8', 280);
    mainLabel.scale.set(11, 3, 1);
    mainLabel.position.copy(mainLabelPos);
    mainLabel.userData.labelPriority = 2;
    this.groups.pcna.add(mainLabel);
    this.groups.pcna.add(this._leaderLine(mainLabelPos, pcnaPos, new THREE.Color(0x7abba8)));

    // Dashed line from thumb center to PCNA
    this.groups.pcna.add(this._dash(thumbCenter, pcnaPos, 0x7abba8));
  }

  // ==================== MUTATION MARKERS ====================

  _buildMutationMarkers() {
    this.groups.mutations = new THREE.Group();
    this.scene.add(this.groups.mutations);

    for (const mut of MUTATIONS) {
      const domCenter = new THREE.Vector3(...DOMAINS[mut.domain].center);
      const pos = domCenter.clone().add(new THREE.Vector3(...mut.offset));
      const g = new THREE.Group();

      const ring = new THREE.Mesh(new THREE.RingGeometry(0.55,0.72,24),
        new THREE.MeshBasicMaterial({ color:PAL.mutMarker, side:THREE.DoubleSide,
          transparent:true, opacity:0.8 }));
      g.add(ring);
      g.add(new THREE.Mesh(new THREE.SphereGeometry(0.25,12,12),
        new THREE.MeshPhysicalMaterial({ color:PAL.mutMarker,
          emissive:PAL.mutMarker, emissiveIntensity:0.8, roughness:0.25, metalness:0.4 })));
      const halo = new THREE.Mesh(new THREE.SphereGeometry(1.0,12,12),
        new THREE.MeshBasicMaterial({ color:PAL.mutMarker, transparent:true, opacity:0.1, depthWrite:false }));
      g.add(halo);
      g.position.copy(pos);
      // Enrich tooltip detail with pipeline data (signatures, classification, DDG)
      let enrichedDetail = mut.detail;
      if (mut.classification) {
        enrichedDetail += ` / ${mut.classification.acmg_class} (P=${mut.classification.posterior})`;
      }
      if (mut.ddg !== undefined) {
        enrichedDetail += ` / \u0394\u0394G=${mut.ddg} kcal/mol`;
      }
      if (mut.signature) {
        const topSig = Object.entries(mut.signature).sort((a,b) => b[1]-a[1])[0];
        if (topSig) enrichedDetail += ` / ${topSig[0]}:${(topSig[1]*100).toFixed(0)}%`;
      }
      g.userData = { halo, ring, name:mut.label, detail:enrichedDetail, mut:mut.id, pulseBoost:0 };
      this.groups.mutations.add(g);
      this.mutMarkers[mut.id] = g;
    }
  }

  // ==================== ANNOTATIONS (domain headers, labels) ====================

  _buildAnnotations() {
    this.groups.labels = new THREE.Group();
    this.scene.add(this.groups.labels);

    // Stagger labels radially outward and high above the model
    const labelOffsets = {
      ntd:     { dx: -8, dy: 16, dz: -6 },
      exo:     { dx: -8, dy: 16, dz:  8 },
      palm:    { dx:-12, dy: 22, dz:-10 },
      pdomain: { dx:  6, dy: 18, dz: -8 },
      fingers: { dx: 10, dy: 20, dz:  8 },
      thumb:   { dx:-10, dy:-16, dz: 14 },
      inactpol:{ dx: 14, dy: 14, dz:  4 },
      ctd:     { dx:  8, dy: 16, dz: -6 },
    };

    for (const key of DOMAIN_ORDER) {
      const dom = DOMAINS[key];
      const col = PAL[dom.color];
      const hex = '#'+col.getHexString();
      const off = labelOffsets[key];
      const anchor = new THREE.Vector3(...dom.center);
      const labelPos = anchor.clone().add(new THREE.Vector3(off.dx, off.dy, off.dz));

      // Domain abbreviated label -high priority
      const label = this._spriteLabel(dom.abbrev, '#d0d8e0', hex);
      label.position.copy(labelPos);
      label.userData.labelPriority = 3;
      this.groups.labels.add(label);

      // Leader line from label down to domain center
      this.groups.labels.add(this._leaderLine(labelPos, anchor, col));

      // Domain header comment -positioned above the abbreviation
      const header = this._spriteLabel(
        `// ${dom.abbrev} (${dom.range[0]}-${dom.range[1]})`,
        hex, hex, 320
      );
      header.scale.set(10, 2.5, 1);
      header.position.copy(labelPos.clone().add(new THREE.Vector3(0, 3.2, 0)));
      header.userData.labelPriority = 1;
      this.groups.labels.add(header);
    }

    const asAnchor = this.activeSitePos.clone();
    const asPos = asAnchor.clone().add(new THREE.Vector3(-10, 20, -10));
    const asLabel = this._spriteLabel('ACTIVE SITE \u00b7 2\u00d7 Mg\u00b2\u207a', '#d0d8e0',
      '#'+PAL.magnesium.getHexString());
    asLabel.scale.set(12,3,1);
    asLabel.position.copy(asPos);
    asLabel.userData.labelPriority = 3;
    this.groups.labels.add(asLabel);
    this.groups.labels.add(this._leaderLine(asPos, asAnchor, PAL.magnesium));

    const proofAnchor = new THREE.Vector3(-10, 2, 7);
    const proofPos = proofAnchor.clone().add(new THREE.Vector3(-8, 18, 8));
    const proofLabel = this._spriteLabel('PROOFREADING \u00b7 ExoI/II/III', '#d0d8e0', '#b85c5c');
    proofLabel.scale.set(12,3,1);
    proofLabel.position.copy(proofPos);
    proofLabel.userData.labelPriority = 2;
    this.groups.labels.add(proofLabel);
    this.groups.labels.add(this._leaderLine(proofPos, proofAnchor, PAL.exo));

    const annotations = [
      { text:'D640xD642 (Motif A)', pos:[0,-1,0], offset:[-12,-14,10], color:'#6b8e6b' },
      { text:'D860 (Motif C)', pos:[0,-1,0], offset:[12,-14,-8], color:'#6b8e6b' },
      { text:'O-helix', pos:[7,5,4], offset:[8,16,6], color:'#5b7f99' },
    ];
    for (const a of annotations) {
      const anchor = new THREE.Vector3(...a.pos);
      const lPos = anchor.clone().add(new THREE.Vector3(...a.offset));
      const lbl = this._spriteLabel(a.text, a.color, a.color, 280);
      lbl.scale.set(9,2.5,1);
      lbl.position.copy(lPos);
      lbl.userData.labelPriority = 1;
      this.groups.labels.add(lbl);
      this.groups.labels.add(this._leaderLine(lPos, anchor, new THREE.Color(a.color)));
    }
  }

  _leaderLine(from, to, color) {
    const geo = new THREE.BufferGeometry().setFromPoints([from, to]);
    const mat = new THREE.LineDashedMaterial({
      color, dashSize: 0.3, gapSize: 0.2,
      transparent: true, opacity: 0.25
    });
    const line = new THREE.Line(geo, mat);
    line.computeLineDistances();
    return line;
  }

  _spriteLabel(text, textColor, accentColor, width) {
    const w = width || 320, h = 96;
    const dpr = 2; // HiDPI for crisp text
    const c = document.createElement('canvas');
    c.width = w * dpr; c.height = h * dpr;
    const ctx = c.getContext('2d');
    ctx.scale(dpr, dpr);
    ctx.clearRect(0,0,w,h);
    ctx.fillStyle = 'rgba(8,11,17,0.9)';
    ctx.fillRect(2,18,w-4,48);
    ctx.fillStyle = accentColor;
    ctx.fillRect(2,18,5,48);
    ctx.font = '600 24px "JetBrains Mono", monospace';
    ctx.fillStyle = textColor;
    ctx.textBaseline = 'middle';
    ctx.fillText(text,18,42);
    ctx.strokeStyle = accentColor;
    ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(0,18); ctx.lineTo(14,10); ctx.stroke();
    const tex = new THREE.CanvasTexture(c);
    tex.colorSpace = THREE.SRGBColorSpace;
    const sprite = new THREE.Sprite(
      new THREE.SpriteMaterial({ map:tex, transparent:true, depthTest:true, depthWrite:false }));
    sprite.scale.set(9,2.8,1);
    sprite.userData.labelPriority = 0;
    sprite.userData.baseOpacity = 1;
    return sprite;
  }

  // ==================== LABEL OVERLAP CULLING ====================

  _updateLabelVisibility() {
    if (!this.groups.labels || !this.groups.labels.visible) return;

    const labels = [];
    this.groups.labels.children.forEach(s => {
      if (s.isSprite) labels.push(s);
    });
    // Also include labels from other visible groups
    for (const gk of ['exoPolPath', 'zincMotifs', 'fesCluster', 'exoSite', 'accessory', 'pcna']) {
      if (this.groups[gk] && this.groups[gk].visible) {
        this.groups[gk].children.forEach(s => {
          if (s.isSprite) labels.push(s);
        });
      }
    }
    if (labels.length === 0) return;

    const w2 = this.renderer.domElement.width / 2;
    const h2 = this.renderer.domElement.height / 2;

    // Project each label to screen space
    const projected = [];
    const v = new THREE.Vector3();
    for (const s of labels) {
      v.setFromMatrixPosition(s.matrixWorld);
      v.project(this.camera);
      // Convert NDC (-1..1) to pixel coords
      projected.push({
        sprite: s,
        x: (v.x + 1) * w2,
        y: (-v.y + 1) * h2,
        z: v.z,
        priority: s.userData.labelPriority || 0,
        baseOpacity: s.userData.baseOpacity || 1,
      });
    }

    // Sort: highest priority first, then closest to camera (lower z)
    projected.sort((a, b) => b.priority - a.priority || a.z - b.z);

    // Check overlap -each label occupies a screen-space rect
    const placed = [];
    const LABEL_W = 130; // approximate half-width in pixels
    const LABEL_H = 35; // approximate half-height in pixels

    for (const p of projected) {
      // Behind camera -hide
      if (p.z > 1 || p.z < -1) {
        p.sprite.material.opacity = 0;
        p.sprite.visible = false;
        continue;
      }

      // Check overlap with already-placed labels
      let overlapping = false;
      for (const q of placed) {
        if (Math.abs(p.x - q.x) < LABEL_W && Math.abs(p.y - q.y) < LABEL_H) {
          overlapping = true;
          break;
        }
      }

      const targetOpacity = overlapping ? 0 : p.baseOpacity;
      // Smooth transition
      const current = p.sprite.material.opacity;
      const newOpacity = current + (targetOpacity - current) * 0.15;
      p.sprite.material.opacity = newOpacity;
      p.sprite.visible = newOpacity > 0.02;

      if (!overlapping) {
        placed.push(p);
      }
    }
  }

  // ==================== HYDROGEN BONDS (i→i+4 in helices) ====================

  _buildHydrogenBonds() {
    this.groups.hbonds = new THREE.Group();
    this.groups.hbonds.visible = true;
    this.scene.add(this.groups.hbonds);

    for (const key of DOMAIN_ORDER) {
      const bb = this.backboneCache[key];
      if (!bb) continue;
      const domColor = PAL[DOMAINS[key].color];

      for (let i = 0; i < bb.length - 4; i++) {
        if (bb[i].ss === 'H' && bb[i+4].ss === 'H') {
          const geo = new THREE.BufferGeometry().setFromPoints([bb[i].pos, bb[i+4].pos]);
          const mat = new THREE.LineDashedMaterial({
            color: domColor, dashSize:0.12, gapSize:0.1,
            transparent:true, opacity:0.35
          });
          const line = new THREE.Line(geo, mat);
          line.computeLineDistances();
          this.groups.hbonds.add(line);
        }
      }
    }
  }

  // ==================== β-SHEET HYDROGEN BONDS ====================

  _buildBetaSheetHBonds() {
    this.groups.sheetHbonds = new THREE.Group();
    this.groups.sheetHbonds.visible = true;
    this.scene.add(this.groups.sheetHbonds);

    for (const key of DOMAIN_ORDER) {
      const bb = this.backboneCache[key];
      if (!bb) continue;
      const segments = segmentBySSType(bb);
      const strands = segments.filter(s => s.type === 'E');

      // Find spatially proximal strand pairs
      for (let i = 0; i < strands.length; i++) {
        for (let j = i + 1; j < strands.length; j++) {
          const sA = strands[i], sB = strands[j];
          // Compute midpoints
          const midA = new THREE.Vector3();
          for (let k = sA.startIdx; k <= sA.endIdx; k++) midA.add(bb[k].pos);
          midA.divideScalar(sA.endIdx - sA.startIdx + 1);
          const midB = new THREE.Vector3();
          for (let k = sB.startIdx; k <= sB.endIdx; k++) midB.add(bb[k].pos);
          midB.divideScalar(sB.endIdx - sB.startIdx + 1);

          if (midA.distanceTo(midB) > 8.0) continue;

          // Draw inter-strand H-bonds every 2nd residue
          const lenA = sA.endIdx - sA.startIdx + 1;
          const lenB = sB.endIdx - sB.startIdx + 1;
          const nBonds = Math.min(lenA, lenB);
          for (let k = 0; k < nBonds; k += 2) {
            const idxA = sA.startIdx + k;
            const idxB = sB.startIdx + Math.min(k, lenB - 1);
            const pA = bb[idxA].pos;
            const pB = bb[idxB].pos;
            const geo = new THREE.BufferGeometry().setFromPoints([pA, pB]);
            const mat = new THREE.LineDashedMaterial({
              color: PAL[DOMAINS[key].color],
              dashSize: 0.12, gapSize: 0.1,
              transparent: true, opacity: 0.3
            });
            const line = new THREE.Line(geo, mat);
            line.computeLineDistances();
            this.groups.sheetHbonds.add(line);
          }
        }
      }
    }
  }

  // ==================== EXO↔POL SWITCHING PATH ====================

  _buildExoPolPath() {
    this.groups.exoPolPath = new THREE.Group();
    this.groups.exoPolPath.visible = true;
    this.scene.add(this.groups.exoPolPath);

    const exoCenter = new THREE.Vector3(...DOMAINS.exo.center);
    const palmCenter = this.activeSitePos.clone();
    const mid = exoCenter.clone().add(palmCenter).multiplyScalar(0.5);
    mid.y += 8;

    const curve = new THREE.CatmullRomCurve3([exoCenter, mid, palmCenter]);
    this.exoPolCurve = curve;

    // Path tube
    const tubeGeo = new THREE.TubeGeometry(curve, 60, 0.06, 6, false);
    const tubeMat = new THREE.MeshPhysicalMaterial({
      color:0x7abba8, transparent:true, opacity:0.3, roughness:0.5
    });
    this.groups.exoPolPath.add(new THREE.Mesh(tubeGeo, tubeMat));

    // Animated sphere
    this.exoPolSphere = new THREE.Mesh(
      new THREE.SphereGeometry(0.3, 12, 12),
      new THREE.MeshPhysicalMaterial({
        color:0x7abba8, emissive:0x7abba8, emissiveIntensity:0.6,
        roughness:0.2, metalness:0.4
      })
    );
    this.exoPolSphere.position.copy(exoCenter);
    this.groups.exoPolPath.add(this.exoPolSphere);

    // Endpoint labels
    const exoLabel = this._spriteLabel('EXO site', '#b85c5c', '#b85c5c', 160);
    exoLabel.scale.set(6, 1.5, 1);
    exoLabel.position.copy(exoCenter.clone().add(new THREE.Vector3(0, 2, 0)));
    this.groups.exoPolPath.add(exoLabel);

    const polLabel = this._spriteLabel('POL site', '#6b8e6b', '#6b8e6b', 160);
    polLabel.scale.set(6, 1.5, 1);
    polLabel.position.copy(palmCenter.clone().add(new THREE.Vector3(0, 2, 0)));
    this.groups.exoPolPath.add(polLabel);

    // Distance annotation (Roske & Yeeles 2024: "almost 40 A" between active sites)
    const distLabel = this._spriteLabel('~40 \u00c5 shuttling distance', '#7abba8', '#7abba8', 320);
    distLabel.scale.set(10, 2, 1);
    distLabel.position.copy(mid.clone().add(new THREE.Vector3(0, 2, 0)));
    this.groups.exoPolPath.add(distLabel);
  }

  // ==================== PARTICLES ====================

  _buildParticles() {
    const n = 150;
    const pos = new Float32Array(n*3);
    for (let i=0;i<n;i++) { pos[i*3]=(Math.random()-0.5)*100; pos[i*3+1]=(Math.random()-0.5)*70; pos[i*3+2]=(Math.random()-0.5)*100; }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(pos,3));
    this.scene.add(new THREE.Points(geo,
      new THREE.PointsMaterial({ color:PAL.accent, size:0.04,
        transparent:true, opacity:0.15, depthWrite:false })));
  }

  // ==================== CONTROLS ====================

  _setupControls() {
    const el = this.renderer.domElement;
    this._shiftDown = false;
    this._pointerDownPos = { x: 0, y: 0 };
    el.addEventListener('pointerdown', e => {
      if (this.measureMode) { this._handleMeasureClick(e); return; }
      this.isDragging=true;
      this._shiftDown = e.shiftKey || e.button === 1;
      this.prevMouse={x:e.clientX,y:e.clientY};
      this._pointerDownPos = { x: e.clientX, y: e.clientY };
    });
    window.addEventListener('pointerup', () => { this.isDragging=false; this._shiftDown=false; });
    el.addEventListener('click', e => {
      if (this.measureMode) return;
      // Only treat as click if pointer didn't move (not a drag)
      const dx = e.clientX - this._pointerDownPos.x;
      const dy = e.clientY - this._pointerDownPos.y;
      if (dx * dx + dy * dy > 25) return;
      this._handleDomainClick(e);
    });
    window.addEventListener('pointermove', e => {
      if (!this.isDragging) return;
      const dx = e.clientX-this.prevMouse.x;
      const dy = e.clientY-this.prevMouse.y;
      if (this._shiftDown || e.shiftKey) {
        // Pan: move look-at target along camera-plane axes
        const panScale = this.cameraRadius * 0.002;
        const rx = Math.cos(this.cameraTheta);
        const rz = -Math.sin(this.cameraTheta);
        this.cameraTarget.x -= dx * panScale * rx;
        this.cameraTarget.z -= dx * panScale * rz;
        this.cameraTarget.y += dy * panScale;
      } else {
        // Orbit
        this.cameraTheta -= dx*0.005;
        this.cameraPhi = Math.max(-1.2,Math.min(1.2, this.cameraPhi+dy*0.005));
      }
      this.prevMouse = {x:e.clientX,y:e.clientY};
      this._updateCamera();
    });
    el.addEventListener('wheel', e => {
      e.preventDefault();
      this.cameraRadius = Math.max(25,Math.min(200,this.cameraRadius+e.deltaY*0.04));
      this._updateCamera();
    }, {passive:false});

    this.raycaster = new THREE.Raycaster();
    this.mouseNDC = new THREE.Vector2();
    el.addEventListener('mousemove', e => this._hover(e));
    el.style.cursor = 'grab';
  }

  _setupTouch() {
    const el = this.renderer.domElement;
    let lastPinchDist = 0;
    let longPressTimer = null;

    el.addEventListener('touchstart', e => {
      e.preventDefault();
      clearTimeout(longPressTimer);
      if (e.touches.length === 1) {
        this.isDragging = true;
        this.prevMouse = { x:e.touches[0].clientX, y:e.touches[0].clientY };
        longPressTimer = setTimeout(() => {
          this._hover({ clientX:e.touches[0].clientX, clientY:e.touches[0].clientY });
        }, 500);
      } else if (e.touches.length === 2) {
        this.isDragging = false;
        lastPinchDist = Math.hypot(
          e.touches[0].clientX - e.touches[1].clientX,
          e.touches[0].clientY - e.touches[1].clientY);
      }
    }, {passive:false});

    let lastTwoFingerMid = null;
    el.addEventListener('touchmove', e => {
      e.preventDefault();
      clearTimeout(longPressTimer);
      if (e.touches.length === 1 && this.isDragging) {
        const dx = e.touches[0].clientX - this.prevMouse.x;
        const dy = e.touches[0].clientY - this.prevMouse.y;
        this.cameraTheta -= dx * 0.005;
        this.cameraPhi = Math.max(-1.2, Math.min(1.2, this.cameraPhi + dy * 0.005));
        this.prevMouse = { x:e.touches[0].clientX, y:e.touches[0].clientY };
        this._updateCamera();
      } else if (e.touches.length === 2) {
        const dist = Math.hypot(
          e.touches[0].clientX - e.touches[1].clientX,
          e.touches[0].clientY - e.touches[1].clientY);
        const delta = lastPinchDist - dist;
        this.cameraRadius = Math.max(25, Math.min(200, this.cameraRadius + delta * 0.1));
        lastPinchDist = dist;

        // Two-finger pan
        const midX = (e.touches[0].clientX + e.touches[1].clientX) / 2;
        const midY = (e.touches[0].clientY + e.touches[1].clientY) / 2;
        if (lastTwoFingerMid) {
          const panDx = midX - lastTwoFingerMid.x;
          const panDy = midY - lastTwoFingerMid.y;
          const panScale = this.cameraRadius * 0.002;
          const rx = Math.cos(this.cameraTheta);
          const rz = -Math.sin(this.cameraTheta);
          this.cameraTarget.x -= panDx * panScale * rx;
          this.cameraTarget.z -= panDx * panScale * rz;
          this.cameraTarget.y += panDy * panScale;
        }
        lastTwoFingerMid = { x: midX, y: midY };
        this._updateCamera();
      }
    }, {passive:false});

    el.addEventListener('touchend', () => {
      clearTimeout(longPressTimer);
      this.isDragging = false;
      lastTwoFingerMid = null;
      this.tooltipEl.classList.remove('show');
    });
  }

  _updateCamera() {
    const {cameraTheta:th,cameraPhi:ph,cameraRadius:r,cameraTarget:tgt} = this;
    this.camera.position.set(
      tgt.x + Math.cos(ph)*Math.sin(th)*r,
      tgt.y + Math.sin(ph)*r,
      tgt.z + Math.cos(ph)*Math.cos(th)*r);
    this.camera.lookAt(tgt);
  }

  _hover(e) {
    const rect = this.renderer.domElement.getBoundingClientRect();
    this.mouseNDC.set(((e.clientX-rect.left)/rect.width)*2-1,
      -((e.clientY-rect.top)/rect.height)*2+1);
    this.raycaster.setFromCamera(this.mouseNDC,this.camera);
    const targets = [];
    this.groups.protein.traverse(o => { if (o.isMesh&&o.userData.domain) targets.push(o); });
    this.groups.metals.traverse(o => { if (o.isMesh&&o.userData.name) targets.push(o); });
    this.groups.mutations.traverse(o => { if (o.isGroup&&o.userData.name) targets.push(o); });
    for (const gk of ['exoSite','zincMotifs','fesCluster','accessory','pcna']) {
      if (this.groups[gk] && this.groups[gk].visible) {
        this.groups[gk].traverse(o => { if (o.isMesh && o.userData.name) targets.push(o); });
      }
    }
    const hits = this.raycaster.intersectObjects(targets,true);
    if (hits.length) {
      let t=hits[0].object; while (t&&!t.userData.name) t=t.parent;
      if (t?.userData.name) {
        this.hoverTarget=t;
        this.tooltipEl.querySelector('.t-name').textContent=t.userData.name;
        this.tooltipEl.querySelector('.t-detail').textContent=t.userData.detail||'';
        this.tooltipEl.classList.add('show');
        this.tooltipEl.style.left=(e.clientX-rect.left+12)+'px';
        this.tooltipEl.style.top=(e.clientY-rect.top+12)+'px';
        this.renderer.domElement.style.cursor='pointer';
        return;
      }
    }
    this.hoverTarget=null;
    this.tooltipEl.classList.remove('show');
    this.renderer.domElement.style.cursor = this.measureMode ? 'crosshair' : (this.isDragging?'grabbing':'grab');
  }

  // ==================== CLICK-TO-SELECT DOMAIN ====================

  _handleDomainClick(e) {
    const rect = this.renderer.domElement.getBoundingClientRect();
    const ndc = new THREE.Vector2(
      ((e.clientX - rect.left) / rect.width) * 2 - 1,
      -((e.clientY - rect.top) / rect.height) * 2 + 1);
    this.raycaster.setFromCamera(ndc, this.camera);

    // Raycast against all ribbon meshes (they have domainKey in userData)
    const targets = [];
    this.groups.protein.traverse(o => {
      if (o.isMesh && (o.userData.domain || o.userData.domainKey)) targets.push(o);
    });
    const hits = this.raycaster.intersectObjects(targets);

    if (hits.length) {
      let obj = hits[0].object;
      const domainKey = obj.userData.domainKey || obj.userData.domain;
      if (domainKey) {
        const alreadySelected = (this._highlightedDomain === domainKey);
        if (alreadySelected) {
          // Deselect and zoom back out
          this.highlightDomain(domainKey, false);
          this.resetView();
        } else {
          // Select this domain and zoom in
          this.highlightDomain(domainKey, true);
          this.focusDomain(domainKey);
        }
        // Notify external listener (sidebar sync)
        if (this.onDomainSelect) {
          this.onDomainSelect(alreadySelected ? null : domainKey);
        }
        return;
      }
    }

    // Clicked empty space -deselect and zoom back out
    if (this._highlightedDomain) {
      this.highlightDomain(this._highlightedDomain, false);
      this.resetView();
      if (this.onDomainSelect) this.onDomainSelect(null);
    }
  }

  // ==================== DISTANCE MEASUREMENT ====================

  _handleMeasureClick(e) {
    const rect = this.renderer.domElement.getBoundingClientRect();
    const ndc = new THREE.Vector2(
      ((e.clientX-rect.left)/rect.width)*2-1,
      -((e.clientY-rect.top)/rect.height)*2+1);
    this.raycaster.setFromCamera(ndc, this.camera);

    const targets = [];
    this.groups.protein.traverse(o => { if (o.isMesh && o.visible) targets.push(o); });
    this.groups.dna.traverse(o => { if (o.isMesh) targets.push(o); });
    this.groups.metals.traverse(o => { if (o.isMesh) targets.push(o); });
    const hits = this.raycaster.intersectObjects(targets);

    if (!hits.length) return;

    this.measurePoints.push(hits[0].point.clone());

    // Show point marker
    const marker = new THREE.Mesh(
      new THREE.SphereGeometry(0.2, 8, 8),
      new THREE.MeshBasicMaterial({ color:0xffaa00 })
    );
    marker.position.copy(hits[0].point);
    marker.userData.isMeasure = true;
    this.scene.add(marker);

    if (this.measurePoints.length === 2) {
      const [p1, p2] = this.measurePoints;
      const dist = p1.distanceTo(p2) / SCALE;

      // Draw dashed line
      const geo = new THREE.BufferGeometry().setFromPoints([p1, p2]);
      const mat = new THREE.LineDashedMaterial({
        color:0xffaa00, dashSize:0.2, gapSize:0.1, linewidth:2
      });
      this.measureLine = new THREE.Line(geo, mat);
      this.measureLine.computeLineDistances();
      this.measureLine.userData.isMeasure = true;
      this.scene.add(this.measureLine);

      // Distance label
      const mid = p1.clone().add(p2).multiplyScalar(0.5);
      this.measureLabel = this._spriteLabel(
        `${dist.toFixed(1)} \u00c5`, '#ffaa00', '#ffaa00', 140);
      this.measureLabel.scale.set(5, 1.5, 1);
      this.measureLabel.position.copy(mid.clone().add(new THREE.Vector3(0, 1, 0)));
      this.measureLabel.userData.isMeasure = true;
      this.scene.add(this.measureLabel);

      this.measurePoints = [];
    }
  }

  toggleMeasure() {
    this.measureMode = !this.measureMode;
    this.measurePoints = [];
    this.renderer.domElement.style.cursor = this.measureMode ? 'crosshair' : 'grab';
    return this.measureMode;
  }

  clearMeasurements() {
    const toRemove = [];
    this.scene.traverse(o => { if (o.userData.isMeasure) toRemove.push(o); });
    for (const o of toRemove) this.scene.remove(o);
    this.measureLine = null;
    this.measureLabel = null;
    this.measurePoints = [];
  }

  // ==================== RENDER LOOP ====================

  _startRenderLoop() {
    this.clock = new THREE.Clock();
    const tick = () => {
      const dt=this.clock.getDelta(), t=this.clock.getElapsedTime();
      if (this.autoRotate) { this.cameraTheta+=dt*0.06; this._updateCamera(); }

      // Mutation marker animation
      this.groups.mutations.children.forEach((g,i) => {
        if (!g.userData?.halo) return;
        const phase=t*1.2+i*0.6;
        g.userData.halo.scale.setScalar(1+Math.sin(phase)*0.15+(g.userData.pulseBoost||0)*0.4);
        g.userData.halo.material.opacity=0.1+Math.abs(Math.sin(phase))*0.08+(g.userData.pulseBoost||0)*0.25;
        g.userData.ring.lookAt(this.camera.position);
        g.userData.ring.rotateZ(t*0.4);
        if (g.userData.pulseBoost>0) g.userData.pulseBoost=Math.max(0,g.userData.pulseBoost-dt*0.8);
      });

      this.activeSiteLight.intensity=0.5+Math.sin(t*2.5)*0.15;
      this.groups.protein.scale.setScalar(1+Math.sin(t*0.4)*0.003);

      // Exo↔Pol switching animation
      if (this.exoPolPlaying && this.exoPolCurve && this.exoPolSphere) {
        this.exoPolT += dt * 0.3;
        if (this.exoPolT > 2) this.exoPolT -= 2;
        const progress = this.exoPolT <= 1 ? this.exoPolT : 2 - this.exoPolT;
        this.exoPolSphere.position.copy(this.exoPolCurve.getPoint(progress));
      }

      // Label overlap culling
      this._updateLabelVisibility();

      // FPS counter
      this.fpsFrames++;
      if (t - this.fpsTime >= 1) {
        if (this.fpsCallback) this.fpsCallback(this.fpsFrames);
        this.fpsFrames = 0;
        this.fpsTime = t;
      }

      this.renderer.render(this.scene,this.camera);
      requestAnimationFrame(tick);
    };
    tick();
  }

  // ==================== PUBLIC API ====================

  setToggle(name,on) {
    const map = { dna:'dna', metals:'metals', dntp:'dntp', labels:'labels',
      mutations:'mutations', surface:'surface', hbonds:'hbonds', exoPolPath:'exoPolPath',
      sheetHbonds:'sheetHbonds', zincMotifs:'zincMotifs', fesCluster:'fesCluster',
      exoSite:'exoSite', accessory:'accessory', pcna:'pcna' };
    if (name==='rotate') { this.autoRotate=on; return; }
    if (name==='exoPolPath') { this.exoPolPlaying=on; }
    const g = this.groups[map[name]]; if (g) g.visible=on;
  }

  highlightDomain(key, active) {
    this._highlightedDomain = active ? key : null;

    for (const dk of DOMAIN_ORDER) {
      const d = this.domainMeshes[dk];
      if (!d?.ribbons) continue;

      const isSelected = (dk === key && active);
      const hasSelection = active; // some domain is selected

      for (const mesh of d.ribbons) {
        if (isSelected) {
          // Highlighted domain: bright glow
          mesh.material.emissiveIntensity = 0.55;
          mesh.material.opacity = 1.0;
          mesh.material.transparent = false;
        } else if (hasSelection) {
          // Other domains: dim and transparent
          mesh.material.emissiveIntensity = 0.02;
          mesh.material.opacity = 0.15;
          mesh.material.transparent = true;
          mesh.material.depthWrite = false;
        } else {
          // No selection: restore defaults
          mesh.material.emissiveIntensity = this.colorMode === 'default' ? 0.10 : 0.12;
          mesh.material.opacity = 1.0;
          mesh.material.transparent = false;
          mesh.material.depthWrite = true;
        }
        mesh.material.needsUpdate = true;
      }
    }
  }

  focusMutation(mutId) {
    const marker=this.mutMarkers[mutId]; if (!marker) return;
    const s={th:this.cameraTheta,ph:this.cameraPhi,r:this.cameraRadius};
    const e={th:-1.2,ph:0,r:35};
    const t0=performance.now();
    const animate=()=>{
      const p=Math.min(1,(performance.now()-t0)/800);
      const ease=p<0.5?2*p*p:1-Math.pow(-2*p+2,2)/2;
      this.cameraTheta=s.th+(e.th-s.th)*ease;
      this.cameraPhi=s.ph+(e.ph-s.ph)*ease;
      this.cameraRadius=s.r+(e.r-s.r)*ease;
      this._updateCamera();
      if (p<1) requestAnimationFrame(animate);
    };
    animate();
    marker.userData.pulseBoost=1;
  }

  focusDomain(key) {
    const dom = DOMAINS[key];
    if (!dom) return;
    const domCenter = new THREE.Vector3(...dom.center);
    const s = { th:this.cameraTheta, ph:this.cameraPhi, r:this.cameraRadius,
      tx:this.cameraTarget.x, ty:this.cameraTarget.y, tz:this.cameraTarget.z };
    const e = { th:s.th, ph:s.ph, r:35, tx:domCenter.x, ty:domCenter.y, tz:domCenter.z };
    const t0 = performance.now();
    const animate = () => {
      const p = Math.min(1, (performance.now()-t0)/800);
      const ease = p<0.5 ? 2*p*p : 1-Math.pow(-2*p+2,2)/2;
      this.cameraTheta = s.th+(e.th-s.th)*ease;
      this.cameraPhi = s.ph+(e.ph-s.ph)*ease;
      this.cameraRadius = s.r+(e.r-s.r)*ease;
      this.cameraTarget.set(
        s.tx+(e.tx-s.tx)*ease, s.ty+(e.ty-s.ty)*ease, s.tz+(e.tz-s.tz)*ease);
      this._updateCamera();
      if (p<1) requestAnimationFrame(animate);
    };
    animate();
  }

  resetView() {
    const s = { th:this.cameraTheta, ph:this.cameraPhi, r:this.cameraRadius,
      tx:this.cameraTarget.x, ty:this.cameraTarget.y, tz:this.cameraTarget.z };
    const e = { th:-0.3, ph:0.15, r:55, tx:0, ty:0, tz:0 };
    const t0 = performance.now();
    const animate = () => {
      const p = Math.min(1, (performance.now()-t0)/600);
      const ease = p<0.5 ? 2*p*p : 1-Math.pow(-2*p+2,2)/2;
      this.cameraTheta = s.th+(e.th-s.th)*ease;
      this.cameraPhi = s.ph+(e.ph-s.ph)*ease;
      this.cameraRadius = s.r+(e.r-s.r)*ease;
      this.cameraTarget.set(
        s.tx+(e.tx-s.tx)*ease, s.ty+(e.ty-s.ty)*ease, s.tz+(e.tz-s.tz)*ease);
      this._updateCamera();
      if (p<1) requestAnimationFrame(animate);
    };
    animate();
  }

  // ==================== COLOR MODES ====================

  setColorMode(mode) {
    this.colorMode = mode;
    for (const key of DOMAIN_ORDER) {
      const dm = this.domainMeshes[key];
      if (!dm?.ribbons) continue;

      for (const mesh of dm.ribbons) {
        const ssType = mesh.userData.ssType;
        let color;

        switch (mode) {
          case 'bfactor': {
            // Use per-residue B-factor from pipeline if available
            let bval;
            if (PER_RESIDUE_SCORES?.bfactor && mesh.userData.startIdx !== undefined) {
              const dom = DOMAINS[key];
              const bb = this.backboneCache[key];
              const totalRes = bb ? bb.length : 1;
              const midIdx = Math.floor((mesh.userData.startIdx + mesh.userData.endIdx) / 2);
              const residue = dom.range[0] + Math.round((midIdx / Math.max(1, totalRes - 1)) * (dom.range[1] - dom.range[0]));
              bval = PER_RESIDUE_SCORES.bfactor[residue - 1] ?? BFACTOR_MAP[ssType] ?? 0.5;
            } else {
              bval = BFACTOR_MAP[ssType] || 0.5;
            }
            color = new THREE.Color().setHSL(0.66 - bval * 0.66, 0.8, 0.5);
            break;
          }
          case 'conservation': {
            // Use per-residue conservation from pipeline if available
            let score;
            if (PER_RESIDUE_SCORES?.conservation && mesh.userData.startIdx !== undefined) {
              const dom = DOMAINS[key];
              const bb = this.backboneCache[key];
              const totalRes = bb ? bb.length : 1;
              const midIdx = Math.floor((mesh.userData.startIdx + mesh.userData.endIdx) / 2);
              const residue = dom.range[0] + Math.round((midIdx / Math.max(1, totalRes - 1)) * (dom.range[1] - dom.range[0]));
              score = PER_RESIDUE_SCORES.conservation[residue - 1] ?? CONSERVATION[key] ?? 0.5;
            } else {
              score = CONSERVATION[key] || 0.5;
            }
            color = new THREE.Color().lerpColors(
              new THREE.Color(0x2288cc), new THREE.Color(0x8822aa), score);
            break;
          }
          case 'heatmap': {
            // Use per-residue variant density from pipeline if available
            let density;
            if (PER_RESIDUE_SCORES?.variant_density && mesh.userData.startIdx !== undefined) {
              const dom = DOMAINS[key];
              const bb = this.backboneCache[key];
              const totalRes = bb ? bb.length : 1;
              const midIdx = Math.floor((mesh.userData.startIdx + mesh.userData.endIdx) / 2);
              const residue = dom.range[0] + Math.round((midIdx / Math.max(1, totalRes - 1)) * (dom.range[1] - dom.range[0]));
              density = PER_RESIDUE_SCORES.variant_density[residue - 1] ?? VARIANT_DENSITY[key] ?? 0.1;
            } else {
              density = VARIANT_DENSITY[key] || 0.1;
            }
            color = new THREE.Color().lerpColors(
              new THREE.Color(0x2244aa), new THREE.Color(0xff2200), density);
            break;
          }
          case 'electrostatic': {
            // Use per-residue charge from pipeline if available
            let charge;
            if (PER_RESIDUE_SCORES?.electrostatic_charge && mesh.userData.startIdx !== undefined) {
              const dom = DOMAINS[key];
              const bb = this.backboneCache[key];
              const totalRes = bb ? bb.length : 1;
              const midIdx = Math.floor((mesh.userData.startIdx + mesh.userData.endIdx) / 2);
              const residue = dom.range[0] + Math.round((midIdx / Math.max(1, totalRes - 1)) * (dom.range[1] - dom.range[0]));
              charge = PER_RESIDUE_SCORES.electrostatic_charge[residue - 1] ?? 0;
            } else {
              const charges = { ntd:0, exo:-0.3, palm:-0.8, pdomain:-0.2, fingers:0.5, thumb:0.6, inactpol:0.1, ctd:0.2 };
              charge = charges[key] || 0;
            }
            if (charge >= 0) color = new THREE.Color().lerpColors(
              new THREE.Color(0xffffff), new THREE.Color(0x2244ff), charge);
            else color = new THREE.Color().lerpColors(
              new THREE.Color(0xffffff), new THREE.Color(0xff2222), -charge);
            break;
          }
          case 'pathogenicity': {
            const bb = this.backboneCache[key];
            const totalRes = bb ? bb.length : 1;
            const midIdx = mesh.userData.startIdx !== undefined
              ? Math.floor((mesh.userData.startIdx + mesh.userData.endIdx) / 2) : 0;
            // Use per-residue pathogenicity from pipeline if available
            let score;
            if (PER_RESIDUE_SCORES?.pathogenicity) {
              const dom = DOMAINS[key];
              const residue = dom.range[0] + Math.round((midIdx / Math.max(1, totalRes - 1)) * (dom.range[1] - dom.range[0]));
              score = PER_RESIDUE_SCORES.pathogenicity[residue - 1] ?? getPathogenicityScore(key, midIdx, totalRes);
            } else {
              score = getPathogenicityScore(key, midIdx, totalRes);
            }
            if (score <= 0.5) {
              const t = score / 0.5;
              color = new THREE.Color().lerpColors(
                new THREE.Color(0x228888), new THREE.Color(0xeeeeee), t);
            } else {
              const t = (score - 0.5) / 0.5;
              color = new THREE.Color().lerpColors(
                new THREE.Color(0xeeeeee), new THREE.Color(0xcc2222), t);
            }
            break;
          }
          default:
            color = PAL[DOMAINS[key].color];
        }

        mesh.material.color.copy(color);
        mesh.material.emissive.copy(color);
        mesh.material.emissiveIntensity = mode === 'default' ? 0.10 : 0.12;
      }
    }

    // Update surface envelope for electrostatic mode
    if (this.groups.surface) {
      this.groups.surface.children.forEach(mesh => {
        const dk = mesh.userData.surfaceDomain;
        if (!dk) return;
        if (mode === 'electrostatic') {
          // Use per-domain charge from pipeline or fallback
          const chargesFallback = { ntd:0, exo:-0.3, palm:-0.8, pdomain:-0.2, fingers:0.5, thumb:0.6, inactpol:0.1, ctd:0.2 };
          const domCharges = (PER_RESIDUE_SCORES && CONSERVATION._charge) || chargesFallback;
          const charge = domCharges[dk] ?? chargesFallback[dk] ?? 0;
          let col;
          if (charge >= 0) col = new THREE.Color().lerpColors(
            new THREE.Color(0xffffff), new THREE.Color(0x2244ff), charge);
          else col = new THREE.Color().lerpColors(
            new THREE.Color(0xffffff), new THREE.Color(0xff2222), -charge);
          mesh.material.color.copy(col);
          mesh.material.opacity = 0.08;
        } else {
          mesh.material.color.copy(PAL[DOMAINS[dk].color]);
          mesh.material.opacity = 0.09;
        }
      });
    }
  }

  setColorScheme(scheme) {
    if (!COLOR_SCHEMES[scheme]) return;
    this.colorScheme = scheme;
    const s = COLOR_SCHEMES[scheme];
    for (const key of DOMAIN_ORDER) {
      PAL[DOMAINS[key].color].setHex(s[key]);
    }
    if (this.colorMode === 'default') this.setColorMode('default');
  }

  // ==================== SCREENSHOT ====================

  takeScreenshot() {
    this.renderer.render(this.scene, this.camera);
    const dataURL = this.renderer.domElement.toDataURL('image/png');
    const a = document.createElement('a');
    a.href = dataURL;
    a.download = `POLE-structure-${Date.now()}.png`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  }

  setFPSCallback(cb) { this.fpsCallback = cb; }
}
