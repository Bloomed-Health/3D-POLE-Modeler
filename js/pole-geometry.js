// ==========================================================================
// POLE Viewer Geometry Utilities
//
// Pure functions for procedural backbone generation: PRNG, domain backbone
// geometry, and inter-domain linker generation. No `this` references —
// these are shared between pole-viewer.js and mutation-viewer.js.
// ==========================================================================

import * as THREE from 'three';

export const SCALE = 0.5;

/** Deterministic PRNG (Park-Miller LCG). */
export function makeRng(seed) {
  let s = seed;
  return () => { s = (s * 16807) % 2147483647; return (s & 0x7fffffff) / 2147483647; };
}

// ==================== INDIVIDUAL SS ELEMENT GENERATORS ====================

/**
 * Generate alpha-helix backbone points.
 * @param {THREE.Vector3} cur - Start position
 * @param {THREE.Vector3} dir - Direction vector
 * @param {number} n - Number of residues
 * @param {number} [sc] - Scale factor (default: SCALE)
 * @returns {Array<{pos: THREE.Vector3, ss: string}>}
 */
export function generateHelixPoints(cur, dir, n, sc) {
  if (sc === undefined) sc = SCALE;
  const axis = dir.clone().normalize();
  let pA = new THREE.Vector3(1,0,0);
  if (Math.abs(axis.dot(pA)) > 0.85) pA.set(0,1,0);
  const pB = new THREE.Vector3().crossVectors(axis,pA).normalize();
  pA.crossVectors(pB,axis).normalize();
  const check = new THREE.Vector3().crossVectors(pA, pB);
  if (check.dot(axis) < 0) pB.negate();
  const pts = [];
  for (let i = 0; i < n; i++) {
    const ang = i * (2*Math.PI/3.6);
    const p = cur.clone()
      .add(axis.clone().multiplyScalar(i * 1.5 * sc))
      .add(pA.clone().multiplyScalar(Math.cos(ang) * 2.3 * sc))
      .add(pB.clone().multiplyScalar(Math.sin(ang) * 2.3 * sc));
    pts.push({ pos:p, ss:'H' });
  }
  return pts;
}

/**
 * Generate beta-strand backbone points.
 * @param {THREE.Vector3} cur - Start position
 * @param {THREE.Vector3} dir - Direction vector
 * @param {number} n - Number of residues
 * @param {number} [sc] - Scale factor (default: SCALE)
 * @returns {Array<{pos: THREE.Vector3, ss: string}>}
 */
export function generateStrandPoints(cur, dir, n, sc) {
  if (sc === undefined) sc = SCALE;
  const sDir = dir.clone().normalize();
  let sA = new THREE.Vector3(1,0,0);
  if (Math.abs(sDir.dot(sA)) > 0.85) sA.set(0,1,0);
  const sB = new THREE.Vector3().crossVectors(sDir, sA).normalize();
  sA.crossVectors(sB, sDir).normalize();
  const pts = [];
  for (let i = 0; i < n; i++) {
    const ang = i * (2 * Math.PI / 3);
    const p = cur.clone()
      .add(sDir.clone().multiplyScalar(i * 3.3 * sc))
      .add(sA.clone().multiplyScalar(Math.cos(ang) * 0.7 * sc))
      .add(sB.clone().multiplyScalar(Math.sin(ang) * 0.7 * sc));
    pts.push({ pos:p, ss:'E' });
  }
  return pts;
}

/**
 * Generate loop/coil backbone points.
 * @param {THREE.Vector3} cur - Start position
 * @param {THREE.Vector3} dir - Direction vector
 * @param {THREE.Vector3} center - Target center for loop attraction
 * @param {number} n - Number of residues
 * @param {function} rng - PRNG function
 * @param {number} [sc] - Scale factor (default: SCALE)
 * @returns {Array<{pos: THREE.Vector3, ss: string}>}
 */
export function generateLoopPoints(cur, dir, center, n, rng, sc) {
  if (sc === undefined) sc = SCALE;
  const target = center.clone().add(new THREE.Vector3(
    (rng()-0.5)*7,(rng()-0.5)*7,(rng()-0.5)*7));
  let loopCur = cur.clone();
  let loopDir = target.clone().sub(cur);
  if (loopDir.length() < 0.01) loopDir.set(rng()-0.5, rng()-0.5, rng()-0.5);
  loopDir.normalize();
  const pts = [];
  for (let i = 0; i < n; i++) {
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
      .normalize().multiplyScalar(3.8 * sc);
    loopCur = loopCur.clone().add(step);
    pts.push({ pos: loopCur.clone(), ss:'L' });
    loopDir = step.clone().normalize();
  }
  return pts;
}

// ==================== DOMAIN BACKBONE GENERATOR ====================

/**
 * Generate backbone points for a domain's secondary structure elements.
 * @param {string} key - Domain key (used for seed)
 * @param {object} domains - Domain definitions object (must have .center and .ss)
 * @param {number} [scale] - Scale factor (default: SCALE)
 * @returns {Array<{pos: THREE.Vector3, ss: string}>}
 */
export function generateDomainBackbone(key, domains, scale) {
  const dom = domains[key];
  const sc = scale !== undefined ? scale : SCALE;
  const center = new THREE.Vector3(...dom.center);
  const rng = makeRng(key.charCodeAt(0) * 137 + key.charCodeAt(1) * 31);
  const points = [];
  let cur = center.clone().add(new THREE.Vector3((rng()-0.5)*5,(rng()-0.5)*5,(rng()-0.5)*5));
  let dir = new THREE.Vector3(rng()-0.5,rng()-0.5,rng()-0.5).normalize();

  for (const el of dom.ss) {
    let pts;
    if (el.t === 'H') pts = generateHelixPoints(cur, dir, el.n, sc);
    else if (el.t === 'E') pts = generateStrandPoints(cur, dir, el.n, sc);
    else pts = generateLoopPoints(cur, dir, center, el.n, rng, sc);
    points.push(...pts);
    cur = points[points.length-1].pos.clone();
    const toC = center.clone().sub(cur).normalize();
    dir = dir.clone().lerp(toC, 0.35+rng()*0.3).normalize();
    dir.add(new THREE.Vector3((rng()-0.5)*0.5,(rng()-0.5)*0.5,(rng()-0.5)*0.5)).normalize();
  }
  return points;
}

/**
 * Generate inter-domain linker points.
 * @param {THREE.Vector3} fromPos
 * @param {THREE.Vector3} toPos
 * @param {number} n - Number of linker points
 * @returns {Array<{pos: THREE.Vector3, ss: string}>}
 */
export function generateLinker(fromPos, toPos, n) {
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
