// ==========================================================================
// POLE c.138del (p.Leu46Phefs*8) -Mutation-Specific 3D Visualization
// Cartoon Ribbon Representation
//
// Biological context:
//   c.138del: single nucleotide deletion at CDS position 138
//   Codon 46 (CTG → Leu) shifts to TGn → Phe (first altered residue)
//   Premature stop at position 54 → only 2.4% of protein translated
//   All catalytic domains absent → complete loss of function
//   Likely NMD substrate (PTC in exon 3, >50 nt upstream of final EJC)
//
// Visualization modes:
//   1. Comparison: full-length wild-type (ghost ribbons) vs truncated mutant
//   2. Domain loss: 3D bar chart of retained vs lost domains
//   3. Sequence: reading frame analysis showing frameshift
//
// References:
//   UniProt Q07864 (POLE_HUMAN)
//   PDB: 9F6D/9F6E/9F6F (Roske & Yeeles 2024), 9B8S (He et al. 2024), 4M8O (Hogg et al. 2014)
//   Mur et al., Genome Med 15:85 (2023) — POLE/POLD1 ACMG classification
// ==========================================================================

import * as THREE from 'three';
import { buildSSElement, segmentBySSType } from './ribbon.js';
import { PAL as BASE_PAL } from './pole-palette.js';
import { SCALE, makeRng, generateLinker } from './pole-geometry.js';

// ==================== PALETTE ====================

const PAL = {
  ...BASE_PAL,
  ghost:    new THREE.Color(0x2a3040),
  intact:   new THREE.Color(0x7abba8),
  frameshift: new THREE.Color(0xb8a04a),
  stop:     new THREE.Color(0xc45c5c),
};

// ==================== SEQUENCE DATA (defaults, overridden by pipeline JSON) ====================

let WT_CODONS = [
  { codon: 'GAG', aa: 'Glu', pos: 43 }, { codon: 'AGC', aa: 'Ser', pos: 44 },
  { codon: 'TCC', aa: 'Ser', pos: 45 }, { codon: 'CTG', aa: 'Leu', pos: 46 },
  { codon: 'GAG', aa: 'Glu', pos: 47 }, { codon: 'TTG', aa: 'Leu', pos: 48 },
  { codon: 'AGG', aa: 'Arg', pos: 49 }, { codon: 'GCA', aa: 'Ala', pos: 50 },
  { codon: 'GCT', aa: 'Ala', pos: 51 }, { codon: 'TTA', aa: 'Leu', pos: 52 },
  { codon: 'CCT', aa: 'Pro', pos: 53 }, { codon: 'TAG', aa: '*',   pos: 54 },
];

let MUT_CODONS = [
  { codon: 'GAG', aa: 'Glu', pos: 43, status: 'normal' },
  { codon: 'AGC', aa: 'Ser', pos: 44, status: 'normal' },
  { codon: 'TCC', aa: 'Ser', pos: 45, status: 'normal' },
  { codon: 'TGA', aa: 'Phe', pos: 46, status: 'frameshift' },
  { codon: 'GGA', aa: 'Gly', pos: 47, status: 'frameshift' },
  { codon: 'GTT', aa: 'Val', pos: 48, status: 'frameshift' },
  { codon: 'GAG', aa: 'Glu', pos: 49, status: 'frameshift' },
  { codon: 'GGC', aa: 'Gly', pos: 50, status: 'frameshift' },
  { codon: 'AGC', aa: 'Ser', pos: 51, status: 'frameshift' },
  { codon: 'TTT', aa: 'Phe', pos: 52, status: 'frameshift' },
  { codon: 'ACC', aa: 'Thr', pos: 53, status: 'frameshift' },
  { codon: 'TAG', aa: '*',   pos: 54, status: 'stop' },
];

let DOMAINS = [
  { key: 'ntd',     label: 'NTD',                         start: 1,    end: 280,  color: PAL.ntd },
  { key: 'exo',     label: "3\u2032\u21925\u2032 Exonuclease", start: 268,  end: 471,  color: PAL.exo },
  { key: 'palm',    label: 'Palm',                         start: 600,  end: 870,  color: PAL.palm },
  { key: 'fingers', label: 'Fingers',                      start: 870,  end: 1020, color: PAL.fingers },
  { key: 'thumb',   label: 'Thumb',                        start: 1020, end: 1190, color: PAL.thumb },
  { key: 'ctd',     label: 'CTD',                          start: 1900, end: 2286, color: PAL.ctd },
];

let TOTAL_RESIDUES = 2286;
let TRUNCATION_SITE = 54;

// ==================== PIPELINE DATA LOADER ====================

/**
 * Load frameshift mutation data from pipeline JSON.
 * Falls back to hardcoded defaults on failure.
 */
export async function loadMutationPipelineData(basePath = './data') {
  try {
    const resp = await fetch(`${basePath}/pole_mutations.json`);
    if (!resp.ok) throw new Error(resp.status);
    const data = await resp.json();

    if (data?.frameshift) {
      const fs = data.frameshift;
      if (fs.wt_codons) WT_CODONS = fs.wt_codons;
      if (fs.mut_codons) MUT_CODONS = fs.mut_codons;
      if (fs.truncation_site) TRUNCATION_SITE = fs.truncation_site;
      if (fs.total_residues) TOTAL_RESIDUES = fs.total_residues;
    }

    console.info('[POLEPipeline] Mutation viewer data loaded');
    return data;
  } catch (e) {
    console.warn('[POLEPipeline] Mutation data not available, using defaults:', e.message);
    return null;
  }
}

// ==================== BACKBONE GENERATION ====================

function generateHelixPoints(cur, dir, n) {
  const axis = dir.clone().normalize();
  let pA = new THREE.Vector3(1, 0, 0);
  if (Math.abs(axis.dot(pA)) > 0.85) pA.set(0, 1, 0);
  const pB = new THREE.Vector3().crossVectors(axis, pA).normalize();
  pA.crossVectors(pB, axis).normalize();
  // Handedness check: ensure right-handed winding
  const check = new THREE.Vector3().crossVectors(pA, pB);
  if (check.dot(axis) < 0) pB.negate();
  const pts = [];
  for (let i = 0; i < n; i++) {
    const ang = i * (2 * Math.PI / 3.6);
    const p = cur.clone()
      .add(axis.clone().multiplyScalar(i * 1.5 * SCALE))
      .add(pA.clone().multiplyScalar(Math.cos(ang) * 2.3 * SCALE))
      .add(pB.clone().multiplyScalar(Math.sin(ang) * 2.3 * SCALE));
    pts.push({ pos: p, ss: 'H' });
  }
  return pts;
}

function generateStrandPoints(cur, dir, n) {
  // Proper 120° pleated-sheet rotation using perpendicular frame
  const sDir = dir.clone().normalize();
  let sA = new THREE.Vector3(1, 0, 0);
  if (Math.abs(sDir.dot(sA)) > 0.85) sA.set(0, 1, 0);
  const sB = new THREE.Vector3().crossVectors(sDir, sA).normalize();
  sA.crossVectors(sB, sDir).normalize();
  const pts = [];
  for (let i = 0; i < n; i++) {
    const ang = i * (2 * Math.PI / 3);
    const p = cur.clone()
      .add(sDir.clone().multiplyScalar(i * 3.3 * SCALE))
      .add(sA.clone().multiplyScalar(Math.cos(ang) * 0.7 * SCALE))
      .add(sB.clone().multiplyScalar(Math.sin(ang) * 0.7 * SCALE));
    pts.push({ pos: p, ss: 'E' });
  }
  return pts;
}

function generateLoopPoints(cur, dir, center, n, rng) {
  // 3.8 Å step chain with dihedral-inspired angular deflection
  const target = center.clone().add(new THREE.Vector3(
    (rng() - 0.5) * 7, (rng() - 0.5) * 7, (rng() - 0.5) * 7));
  let loopCur = cur.clone();
  let loopDir = target.clone().sub(cur);
  if (loopDir.length() < 0.01) loopDir.set(rng() - 0.5, rng() - 0.5, rng() - 0.5);
  loopDir.normalize();
  const pts = [];
  for (let i = 0; i < n; i++) {
    const toTarget = target.clone().sub(loopCur);
    if (toTarget.length() > 0.01) loopDir.lerp(toTarget.normalize(), 0.3).normalize();
    const phi = (rng() - 0.5) * Math.PI * 0.6;
    const psi = (rng() - 0.5) * Math.PI * 0.4;
    let perpSeed = new THREE.Vector3(1, 0, 0);
    if (Math.abs(loopDir.dot(perpSeed)) > 0.85) perpSeed.set(0, 1, 0);
    const perp1 = new THREE.Vector3().crossVectors(loopDir, perpSeed).normalize();
    const perp2 = new THREE.Vector3().crossVectors(loopDir, perp1).normalize();
    const step = loopDir.clone()
      .add(perp1.clone().multiplyScalar(Math.sin(phi) * 0.5))
      .add(perp2.clone().multiplyScalar(Math.sin(psi) * 0.3))
      .normalize().multiplyScalar(3.8 * SCALE);
    loopCur = loopCur.clone().add(step);
    pts.push({ pos: loopCur.clone(), ss: 'L' });
    loopDir = step.clone().normalize();
  }
  return pts;
}

// Domain backbone geometry (secondary structure elements)
const DOMAIN_SS = {
  ntd: [
    { t:'H',n:14 },{ t:'L',n:6 },{ t:'E',n:7 },{ t:'L',n:4 },{ t:'E',n:6 },{ t:'L',n:5 },
    { t:'H',n:18 },{ t:'L',n:5 },{ t:'E',n:8 },{ t:'L',n:4 },{ t:'E',n:7 },{ t:'L',n:6 },
    { t:'H',n:12 },{ t:'L',n:8 },{ t:'H',n:10 },{ t:'L',n:5 },{ t:'E',n:6 },{ t:'L',n:4 },
    { t:'H',n:16 },{ t:'L',n:7 },{ t:'H',n:10 },{ t:'L',n:8 },{ t:'E',n:5 },{ t:'L',n:6 },
    { t:'H',n:8 },{ t:'L',n:9 },
  ],
  exo: [
    { t:'E',n:6 },{ t:'L',n:4 },{ t:'H',n:14 },{ t:'L',n:5 },{ t:'E',n:7 },{ t:'L',n:3 },
    { t:'E',n:6 },{ t:'H',n:12 },{ t:'L',n:6 },{ t:'H',n:16 },{ t:'L',n:4 },{ t:'E',n:8 },
    { t:'L',n:5 },{ t:'H',n:10 },{ t:'L',n:6 },{ t:'E',n:7 },{ t:'L',n:3 },{ t:'H',n:14 },
    { t:'L',n:5 },{ t:'H',n:8 },{ t:'L',n:6 },{ t:'E',n:6 },{ t:'L',n:4 },{ t:'H',n:10 },
    { t:'L',n:4 },
  ],
  palm: [
    { t:'E',n:8 },{ t:'L',n:3 },{ t:'E',n:7 },{ t:'L',n:5 },{ t:'H',n:12 },{ t:'L',n:4 },
    { t:'E',n:9 },{ t:'L',n:3 },{ t:'E',n:7 },{ t:'L',n:6 },{ t:'H',n:10 },{ t:'L',n:4 },
    { t:'E',n:8 },{ t:'L',n:5 },{ t:'E',n:6 },{ t:'L',n:3 },{ t:'H',n:14 },{ t:'L',n:5 },
    { t:'E',n:7 },{ t:'L',n:4 },{ t:'H',n:8 },{ t:'L',n:6 },{ t:'E',n:9 },{ t:'L',n:3 },
    { t:'E',n:7 },{ t:'L',n:5 },{ t:'H',n:12 },{ t:'L',n:4 },
  ],
  fingers: [
    { t:'H',n:20 },{ t:'L',n:4 },{ t:'H',n:16 },{ t:'L',n:6 },{ t:'H',n:22 },{ t:'L',n:5 },
    { t:'H',n:14 },{ t:'L',n:4 },{ t:'H',n:12 },{ t:'L',n:5 },{ t:'H',n:10 },{ t:'L',n:6 },
    { t:'H',n:8 },{ t:'L',n:4 },
  ],
  thumb: [
    { t:'H',n:18 },{ t:'L',n:4 },{ t:'H',n:22 },{ t:'L',n:5 },{ t:'H',n:16 },{ t:'L',n:6 },
    { t:'H',n:20 },{ t:'L',n:4 },{ t:'H',n:14 },{ t:'L',n:5 },{ t:'H',n:12 },{ t:'L',n:6 },
    { t:'H',n:10 },{ t:'L',n:4 },
  ],
  ctd: [
    { t:'H',n:10 },{ t:'L',n:6 },{ t:'E',n:5 },{ t:'L',n:4 },{ t:'E',n:6 },{ t:'L',n:3 },
    { t:'H',n:14 },{ t:'L',n:5 },{ t:'E',n:7 },{ t:'L',n:4 },{ t:'E',n:5 },{ t:'L',n:6 },
    { t:'H',n:12 },{ t:'L',n:4 },{ t:'E',n:6 },{ t:'L',n:5 },{ t:'H',n:10 },{ t:'L',n:8 },
    { t:'E',n:5 },{ t:'L',n:4 },{ t:'H',n:8 },{ t:'L',n:6 },{ t:'E',n:7 },{ t:'L',n:3 },
    { t:'H',n:12 },{ t:'L',n:5 },{ t:'E',n:6 },{ t:'L',n:4 },{ t:'H',n:10 },{ t:'L',n:7 },
    { t:'E',n:5 },{ t:'L',n:5 },{ t:'H',n:14 },{ t:'L',n:6 },
  ],
};

const DOMAIN_CENTERS = {
  ntd:     [-18, 8, -5],
  exo:     [-10, 2, 7],
  palm:    [0, -1, 0],
  fingers: [7, 5, 4],
  thumb:   [3, -9, 7],
  ctd:     [16, 1, -5],
};

const DOMAIN_ORDER = ['ntd', 'exo', 'palm', 'fingers', 'thumb', 'ctd'];

function generateDomainBackbone(key) {
  const center = new THREE.Vector3(...DOMAIN_CENTERS[key]);
  const ss = DOMAIN_SS[key];
  const rng = makeRng(key.charCodeAt(0) * 137 + key.charCodeAt(1) * 31);
  const points = [];
  let cur = center.clone().add(new THREE.Vector3(
    (rng() - 0.5) * 5, (rng() - 0.5) * 5, (rng() - 0.5) * 5));
  let dir = new THREE.Vector3(rng() - 0.5, rng() - 0.5, rng() - 0.5).normalize();

  for (const el of ss) {
    let pts;
    if (el.t === 'H') pts = generateHelixPoints(cur, dir, el.n);
    else if (el.t === 'E') pts = generateStrandPoints(cur, dir, el.n);
    else pts = generateLoopPoints(cur, dir, center, el.n, rng);
    points.push(...pts);
    cur = points[points.length - 1].pos.clone();
    const toC = center.clone().sub(cur).normalize();
    dir = dir.clone().lerp(toC, 0.35 + rng() * 0.3).normalize();
    dir.add(new THREE.Vector3(
      (rng() - 0.5) * 0.5, (rng() - 0.5) * 0.5, (rng() - 0.5) * 0.5)).normalize();
  }
  return points;
}

// ==================== MAIN CLASS ====================

export class MutationViewer {
  constructor(containerEl, tooltipEl) {
    this.container = containerEl;
    this.tooltipEl = tooltipEl;
    this.autoRotate = true;
    this.isDragging = false;
    this.prevMouse = { x: 0, y: 0 };
    this.cameraTheta = -0.4;
    this.cameraPhi = 0.2;
    this.cameraRadius = 95;
    this.mode = 'comparison';
    this.cameraTarget = new THREE.Vector3(0, 0, 0);
    this.groups = {};
    this.fpsFrames = 0;
    this.fpsTime = 0;
    this.fpsCallback = null;
    this._init();
  }

  _init() {
    this._createScene();
    this._createLighting();
    this._buildComparison();
    this._buildDomainLossChart();
    this._buildDivider();
    this._setupControls();
    this._setupTouch();
    this._startRenderLoop();
  }

  _createScene() {
    this.scene = new THREE.Scene();
    this.scene.fog = new THREE.FogExp2(PAL.bg, 0.003);
    this.camera = new THREE.PerspectiveCamera(40, 1, 0.1, 2000);
    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
    this.renderer.toneMappingExposure = 1.4;
    this.container.appendChild(this.renderer.domElement);
    const resize = () => {
      const w = this.container.clientWidth, h = this.container.clientHeight;
      this.renderer.setSize(w, h, false);
      this.camera.aspect = w / h;
      this.camera.updateProjectionMatrix();
    };
    new ResizeObserver(resize).observe(this.container);
    resize();
    this._updateCamera();
  }

  _createLighting() {
    this.scene.add(new THREE.AmbientLight(0x8090a0, 0.8));
    const key = new THREE.DirectionalLight(0xeee8dd, 1.3);
    key.position.set(18, 25, 20);
    this.scene.add(key);
    const fill = new THREE.DirectionalLight(0x6a7a8a, 0.4);
    fill.position.set(-10, -15, 10);
    this.scene.add(fill);
    const rim = new THREE.DirectionalLight(0x7abba8, 0.5);
    rim.position.set(-20, 8, -15);
    this.scene.add(rim);
    this.truncLight = new THREE.PointLight(0xc45c5c, 0.7, 30);
    this.scene.add(this.truncLight);
  }

  // ==================== COMPARISON VIEW ====================

  _buildComparison() {
    this.groups.wildtype = new THREE.Group();
    this.groups.mutant = new THREE.Group();
    this.scene.add(this.groups.wildtype);
    this.scene.add(this.groups.mutant);

    this._buildWTGhost(new THREE.Vector3(-14, 0, 0));
    this._buildTruncatedMutant(new THREE.Vector3(14, 0, 0));
  }

  _buildWTGhost(offset) {
    const backbones = {};

    for (const key of DOMAIN_ORDER) {
      const bb = generateDomainBackbone(key);
      backbones[key] = bb;

      const segments = segmentBySSType(bb);
      for (let s = 0; s < segments.length; s++) {
        const seg = segments[s];
        const lo = Math.max(0, seg.startIdx - (s > 0 ? 1 : 0));
        const hi = Math.min(bb.length - 1, seg.endIdx + (s < segments.length - 1 ? 1 : 0));
        const pts = [];
        for (let i = lo; i <= hi; i++) pts.push(bb[i].pos.clone().add(offset));
        if (pts.length >= 2) {
          const mesh = buildSSElement(pts, seg.type, PAL.ghost, 0.15);
          this.groups.wildtype.add(mesh);
        }
      }
    }

    for (let i = 0; i < DOMAIN_ORDER.length - 1; i++) {
      const kA = DOMAIN_ORDER[i], kB = DOMAIN_ORDER[i + 1];
      const bbA = backbones[kA], bbB = backbones[kB];
      const domA = DOMAINS.find(d => d.key === kA), domB = DOMAINS.find(d => d.key === kB);
      const linkerN = Math.max(8, Math.floor((domB.start - domA.end) / 5));
      const linkerBB = generateLinker(
        bbA[bbA.length - 1].pos.clone().add(offset),
        bbB[0].pos.clone().add(offset), linkerN);
      const linkerPts = linkerBB.map(p => p.pos);
      if (linkerPts.length >= 2) {
        const mesh = buildSSElement(linkerPts, 'L', PAL.ghost, 0.15);
        this.groups.wildtype.add(mesh);
      }
    }

    const label = this._label('WILD-TYPE (2,286 aa)', '#7a8898', '#4a5868');
    label.scale.set(12, 2.5, 1);
    label.position.set(offset.x, 14, offset.z);
    this.groups.wildtype.add(label);

    for (const key of DOMAIN_ORDER) {
      const c = DOMAIN_CENTERS[key];
      const dl = this._label(key.toUpperCase(), '#4a5868', '#3a4858');
      dl.scale.set(4, 1.3, 1);
      dl.position.set(c[0] + offset.x, c[1] + 7 + offset.y, c[2] + offset.z);
      this.groups.wildtype.add(dl);
    }
  }

  _buildTruncatedMutant(offset) {
    const rng = makeRng(99);

    const elements = [
      { type: 'L', n: 4 },
      { type: 'H', n: 14 },
      { type: 'L', n: 7 },
      { type: 'E', n: 7 },
      { type: 'L', n: 6 },
      { type: 'H', n: 16 },
    ];

    const tracePoints = [];
    let cur = new THREE.Vector3(-3 + offset.x, 2 + offset.y, offset.z);
    let dir = new THREE.Vector3(1, -0.3, 0.2).normalize();
    const traceCenter = cur.clone();

    for (const el of elements) {
      const n = Math.min(el.n, 54 - tracePoints.length);
      if (n <= 0) break;
      let pts;
      if (el.type === 'H') pts = generateHelixPoints(cur, dir, n);
      else if (el.type === 'E') pts = generateStrandPoints(cur, dir, n);
      else pts = generateLoopPoints(cur, dir, traceCenter, n, rng);

      for (const pt of pts) {
        if (tracePoints.length >= 54) break;
        pt.pos.add(offset).sub(traceCenter).add(traceCenter);
        tracePoints.push(pt);
      }
      if (tracePoints.length >= 54) break;
      cur = tracePoints[tracePoints.length - 1].pos.clone();
      dir.add(new THREE.Vector3(
        (rng() - 0.5) * 0.5, (rng() - 0.5) * 0.5, (rng() - 0.5) * 0.3)).normalize();
    }

    // Color zones: 0-44 intact (teal), 45-52 frameshift (gold), 53 stop (red)
    const colorBounds = [
      { end: 44, color: PAL.intact },
      { end: 52, color: PAL.frameshift },
      { end: 53, color: PAL.stop },
    ];

    const ssSegs = segmentBySSType(tracePoints);
    for (const seg of ssSegs) {
      let lo = seg.startIdx;
      for (const cb of colorBounds) {
        if (lo > seg.endIdx) break;
        const hi = Math.min(seg.endIdx, cb.end);
        if (hi < lo) continue;
        const ptsLo = Math.max(0, lo - 1);
        const ptsHi = Math.min(tracePoints.length - 1, hi + 1);
        const pts = [];
        for (let i = ptsLo; i <= ptsHi; i++) pts.push(tracePoints[i].pos);
        if (pts.length >= 2) {
          const mesh = buildSSElement(pts, seg.type, cb.color);
          this.groups.mutant.add(mesh);
        }
        lo = hi + 1;
      }
    }

    if (tracePoints.length > 0) {
      const terminus = tracePoints[tracePoints.length - 1].pos;

      const stopSphere = new THREE.Mesh(
        new THREE.SphereGeometry(0.4, 16, 16),
        new THREE.MeshPhysicalMaterial({
          color: 0xc45c5c, emissive: 0xc45c5c, emissiveIntensity: 0.7,
          roughness: 0.2, metalness: 0.5,
        })
      );
      stopSphere.position.copy(terminus);
      this.groups.mutant.add(stopSphere);

      const ring = new THREE.Mesh(
        new THREE.RingGeometry(0.6, 0.85, 24),
        new THREE.MeshBasicMaterial({
          color: PAL.stop, side: THREE.DoubleSide, transparent: true, opacity: 0.6 })
      );
      ring.position.copy(terminus);
      this.stopRing = ring;
      this.groups.mutant.add(ring);

      this.truncLight.position.copy(terminus);

      const tLabel = this._label('PTC @ res 54 / 2,286', '#c45c5c', '#c45c5c');
      tLabel.scale.set(12, 2.5, 1);
      tLabel.position.copy(terminus.clone().add(new THREE.Vector3(0, 3.5, 0)));
      this.groups.mutant.add(tLabel);

      if (tracePoints.length > 45) {
        const res46 = tracePoints[45].pos;
        const fsLabel = this._label('c.138del frameshift', '#b8a04a', '#b8a04a');
        fsLabel.scale.set(10, 2, 1);
        fsLabel.position.copy(res46.clone().add(new THREE.Vector3(-1, 2.5, 0)));
        this.groups.mutant.add(fsLabel);

        const dashGeo = new THREE.BufferGeometry().setFromPoints([
          res46.clone().add(new THREE.Vector3(-1, 2, 0)), res46]);
        const dashMat = new THREE.LineDashedMaterial({
          color: 0xb8a04a, dashSize: 0.2, gapSize: 0.15, transparent: true, opacity: 0.4 });
        const dash = new THREE.Line(dashGeo, dashMat);
        dash.computeLineDistances();
        this.groups.mutant.add(dash);
      }
    }

    this._buildInlineGhost(offset);

    const label = this._label('c.138del \u00b7 p.Leu46Phefs*8', '#c45c5c', '#c45c5c');
    label.scale.set(14, 2.8, 1);
    label.position.set(offset.x, 14, offset.z);
    this.groups.mutant.add(label);
  }

  _buildInlineGhost(offset) {
    for (const key of DOMAIN_ORDER) {
      const bb = generateDomainBackbone(key);
      const segments = segmentBySSType(bb);
      for (let s = 0; s < segments.length; s++) {
        const seg = segments[s];
        const lo = Math.max(0, seg.startIdx - (s > 0 ? 1 : 0));
        const hi = Math.min(bb.length - 1, seg.endIdx + (s < segments.length - 1 ? 1 : 0));
        const pts = [];
        for (let i = lo; i <= hi; i++) pts.push(bb[i].pos.clone().add(offset));
        if (pts.length >= 2) {
          const mesh = buildSSElement(pts, seg.type, PAL.ghost, 0.06);
          this.groups.mutant.add(mesh);
        }
      }
    }
  }

  // ==================== DIVIDER ====================

  _buildDivider() {
    this.groups.truncMarker = new THREE.Group();
    this.scene.add(this.groups.truncMarker);
    const geo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, -15, 0), new THREE.Vector3(0, 15, 0),
    ]);
    const mat = new THREE.LineDashedMaterial({
      color: 0x7abba8, dashSize: 0.4, gapSize: 0.3, transparent: true, opacity: 0.25,
    });
    const line = new THREE.Line(geo, mat);
    line.computeLineDistances();
    this.groups.truncMarker.add(line);
    const vsLabel = this._label('vs', '#7a8898', '#4a5868');
    vsLabel.scale.set(3, 1.5, 1);
    vsLabel.position.set(0, 14, 0);
    this.groups.truncMarker.add(vsLabel);
  }

  // ==================== DOMAIN LOSS CHART ====================

  _buildDomainLossChart() {
    this.groups.domainLoss = new THREE.Group();
    this.groups.domainLoss.visible = false;
    this.scene.add(this.groups.domainLoss);

    const barW = 2.5, spacing = 3.5, maxH = 12;
    for (let i = 0; i < DOMAINS.length; i++) {
      const d = DOMAINS[i];
      const x = (i - DOMAINS.length / 2 + 0.5) * spacing;
      const domSize = d.end - d.start;
      const h = (domSize / TOTAL_RESIDUES) * maxH * 8;

      const fullBar = new THREE.Mesh(
        new THREE.BoxGeometry(barW, h, barW * 0.6),
        new THREE.MeshPhysicalMaterial({
          color: d.color, roughness: 0.6, transparent: true, opacity: 0.12,
        })
      );
      fullBar.position.set(x, h / 2 - 8, 0);
      this.groups.domainLoss.add(fullBar);

      const retained = Math.max(0, Math.min(TRUNCATION_SITE, d.end) - d.start);
      const retFrac = retained / domSize;
      const retH = retFrac * h;
      if (retH > 0.01) {
        const retBar = new THREE.Mesh(
          new THREE.BoxGeometry(barW, retH, barW * 0.6),
          new THREE.MeshPhysicalMaterial({
            color: retFrac > 0.5 ? 0x7abba8 : 0xc45c5c, roughness: 0.4, metalness: 0.1,
          })
        );
        retBar.position.set(x, retH / 2 - 8, 0);
        this.groups.domainLoss.add(retBar);
      }

      const dl = this._label(d.label, '#d0d8e0', '#' + d.color.getHexString());
      dl.scale.set(5, 1.5, 1);
      dl.position.set(x, -9.5, 0);
      this.groups.domainLoss.add(dl);

      const pct = Math.round(retFrac * 100);
      const pl = this._label(`${pct}%`, pct > 0 ? '#7abba8' : '#c45c5c', pct > 0 ? '#7abba8' : '#c45c5c');
      pl.scale.set(3, 1.2, 1);
      pl.position.set(x, h - 7.5, 0);
      this.groups.domainLoss.add(pl);
    }

    const title = this._label('DOMAIN RETENTION ANALYSIS', '#d0d8e0', '#7abba8');
    title.scale.set(16, 3, 1);
    title.position.set(0, 8, 0);
    this.groups.domainLoss.add(title);

    const sub = this._label('54 / 2,286 residues translated (2.4%)', '#c45c5c', '#c45c5c');
    sub.scale.set(18, 2.5, 1);
    sub.position.set(0, 5.5, 0);
    this.groups.domainLoss.add(sub);
  }

  // ==================== LABEL HELPER ====================

  _label(text, textColor, accentColor) {
    const w = 400, h = 96;
    const dpr = 2; // HiDPI for crisp text
    const c = document.createElement('canvas');
    c.width = w * dpr; c.height = h * dpr;
    const ctx = c.getContext('2d');
    ctx.scale(dpr, dpr);
    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = 'rgba(8,11,17,0.9)';
    ctx.fillRect(2, 18, w - 4, 48);
    ctx.fillStyle = accentColor;
    ctx.fillRect(2, 18, 5, 48);
    ctx.font = '600 24px "JetBrains Mono", monospace';
    ctx.fillStyle = textColor;
    ctx.textBaseline = 'middle';
    ctx.fillText(text, 18, 42);
    const tex = new THREE.CanvasTexture(c);
    tex.colorSpace = THREE.SRGBColorSpace;
    const sprite = new THREE.Sprite(
      new THREE.SpriteMaterial({ map: tex, transparent: true, depthTest: false })
    );
    sprite.scale.set(9, 2.8, 1);
    return sprite;
  }

  // ==================== CONTROLS ====================

  _setupControls() {
    const el = this.renderer.domElement;
    this._shiftDown = false;
    el.addEventListener('pointerdown', (e) => {
      this.isDragging = true;
      this._shiftDown = e.shiftKey || e.button === 1;
      this.prevMouse = { x: e.clientX, y: e.clientY };
    });
    window.addEventListener('pointerup', () => { this.isDragging = false; this._shiftDown = false; });
    window.addEventListener('pointermove', (e) => {
      if (!this.isDragging) return;
      const dx = e.clientX - this.prevMouse.x;
      const dy = e.clientY - this.prevMouse.y;
      if (this._shiftDown || e.shiftKey) {
        const panScale = this.cameraRadius * 0.002;
        const rx = Math.cos(this.cameraTheta);
        const rz = -Math.sin(this.cameraTheta);
        this.cameraTarget.x -= dx * panScale * rx;
        this.cameraTarget.z -= dx * panScale * rz;
        this.cameraTarget.y += dy * panScale;
      } else {
        this.cameraTheta -= dx * 0.005;
        this.cameraPhi = Math.max(-1.2, Math.min(1.2, this.cameraPhi + dy * 0.005));
      }
      this.prevMouse = { x: e.clientX, y: e.clientY };
      this._updateCamera();
    });
    el.addEventListener('wheel', (e) => {
      e.preventDefault();
      this.cameraRadius = Math.max(20, Math.min(200, this.cameraRadius + e.deltaY * 0.04));
      this._updateCamera();
    }, { passive: false });
    el.style.cursor = 'grab';
  }

  _setupTouch() {
    const el = this.renderer.domElement;
    let lastPinchDist = 0;

    el.addEventListener('touchstart', e => {
      e.preventDefault();
      if (e.touches.length === 1) {
        this.isDragging = true;
        this.prevMouse = { x: e.touches[0].clientX, y: e.touches[0].clientY };
      } else if (e.touches.length === 2) {
        this.isDragging = false;
        lastPinchDist = Math.hypot(
          e.touches[0].clientX - e.touches[1].clientX,
          e.touches[0].clientY - e.touches[1].clientY);
      }
    }, { passive: false });

    let lastTwoFingerMid = null;
    el.addEventListener('touchmove', e => {
      e.preventDefault();
      if (e.touches.length === 1 && this.isDragging) {
        const dx = e.touches[0].clientX - this.prevMouse.x;
        const dy = e.touches[0].clientY - this.prevMouse.y;
        this.cameraTheta -= dx * 0.005;
        this.cameraPhi = Math.max(-1.2, Math.min(1.2, this.cameraPhi + dy * 0.005));
        this.prevMouse = { x: e.touches[0].clientX, y: e.touches[0].clientY };
        this._updateCamera();
      } else if (e.touches.length === 2) {
        const dist = Math.hypot(
          e.touches[0].clientX - e.touches[1].clientX,
          e.touches[0].clientY - e.touches[1].clientY);
        const delta = lastPinchDist - dist;
        this.cameraRadius = Math.max(20, Math.min(200, this.cameraRadius + delta * 0.1));
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
    }, { passive: false });

    el.addEventListener('touchend', () => { this.isDragging = false; lastTwoFingerMid = null; });
  }

  _updateCamera() {
    const { cameraTheta: th, cameraPhi: ph, cameraRadius: r, cameraTarget: tgt } = this;
    this.camera.position.set(
      tgt.x + Math.cos(ph) * Math.sin(th) * r,
      tgt.y + Math.sin(ph) * r,
      tgt.z + Math.cos(ph) * Math.cos(th) * r
    );
    this.camera.lookAt(tgt);
  }

  // ==================== RENDER LOOP ====================

  _startRenderLoop() {
    this.clock = new THREE.Clock();
    const tick = () => {
      const dt = this.clock.getDelta();
      const t = this.clock.getElapsedTime();
      if (this.autoRotate) {
        this.cameraTheta += dt * 0.04;
        this._updateCamera();
      }
      if (this.truncLight) this.truncLight.intensity = 0.5 + Math.sin(t * 2) * 0.25;
      if (this.stopRing) {
        this.stopRing.lookAt(this.camera.position);
        this.stopRing.material.opacity = 0.4 + Math.sin(t * 1.8) * 0.2;
      }

      // FPS counter
      this.fpsFrames++;
      if (t - this.fpsTime >= 1) {
        if (this.fpsCallback) this.fpsCallback(this.fpsFrames);
        this.fpsFrames = 0;
        this.fpsTime = t;
      }

      this.renderer.render(this.scene, this.camera);
      requestAnimationFrame(tick);
    };
    tick();
  }

  // ==================== PUBLIC API ====================

  setMode(mode) {
    this.mode = mode;
    this.groups.wildtype.visible = mode === 'comparison';
    this.groups.mutant.visible = mode === 'comparison';
    this.groups.truncMarker.visible = mode === 'comparison';
    this.groups.domainLoss.visible = mode === 'domainLoss';
  }

  setAutoRotate(on) { this.autoRotate = on; }

  resetView() {
    const s = { th: this.cameraTheta, ph: this.cameraPhi, r: this.cameraRadius,
      tx: this.cameraTarget.x, ty: this.cameraTarget.y, tz: this.cameraTarget.z };
    const e = { th: -0.4, ph: 0.2, r: 60, tx: 0, ty: 0, tz: 0 };
    const t0 = performance.now();
    const animate = () => {
      const p = Math.min(1, (performance.now() - t0) / 600);
      const ease = p < 0.5 ? 2*p*p : 1 - Math.pow(-2*p+2, 2) / 2;
      this.cameraTheta = s.th + (e.th - s.th) * ease;
      this.cameraPhi = s.ph + (e.ph - s.ph) * ease;
      this.cameraRadius = s.r + (e.r - s.r) * ease;
      this.cameraTarget.set(
        s.tx + (e.tx - s.tx) * ease, s.ty + (e.ty - s.ty) * ease, s.tz + (e.tz - s.tz) * ease);
      this._updateCamera();
      if (p < 1) requestAnimationFrame(animate);
    };
    animate();
  }

  takeScreenshot() {
    this.renderer.render(this.scene, this.camera);
    const dataURL = this.renderer.domElement.toDataURL('image/png');
    const a = document.createElement('a');
    a.href = dataURL;
    a.download = `POLE-mutation-${Date.now()}.png`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  }

  setFPSCallback(cb) { this.fpsCallback = cb; }

  static getSequenceData() { return { wildtype: WT_CODONS, mutant: MUT_CODONS }; }

  static getDomainLossStats() {
    return DOMAINS.map(d => {
      const retained = Math.max(0, Math.min(TRUNCATION_SITE, d.end) - d.start);
      const total = d.end - d.start;
      return { key: d.key, label: d.label, retained, total,
        pct: Math.round((retained / total) * 100), lost: retained === 0 };
    });
  }
}
