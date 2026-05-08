// ==========================================================================
// Shared Ribbon Geometry Utilities
// Cartoon ribbon rendering for protein secondary structure visualization.
// Used by both pole-viewer.js and mutation-viewer.js.
// ==========================================================================

import * as THREE from 'three';

/**
 * Sweep an elliptical cross-section along a CatmullRom curve.
 * Returns indexed BufferGeometry with positions and normals.
 */
export function sweepRibbon(curvePoints, halfW, halfH, profileN, taperFn) {
  if (curvePoints.length < 2) return new THREE.BufferGeometry();

  const curve = new THREE.CatmullRomCurve3(curvePoints);
  const segments = Math.max(curvePoints.length * 4, 16);
  const frames = curve.computeFrenetFrames(segments, false);

  const vertices = [];
  const normals = [];
  const indices = [];
  const ringSize = profileN + 1;

  for (let i = 0; i <= segments; i++) {
    const t = i / segments;
    const pos = curve.getPoint(t);
    const N = frames.normals[i];
    const B = frames.binormals[i];
    const taper = taperFn ? taperFn(t) : 1;
    const w = halfW * taper;
    const h = halfH * taper;

    for (let j = 0; j <= profileN; j++) {
      const angle = (j / profileN) * Math.PI * 2;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);

      vertices.push(
        pos.x + N.x * cosA * w + B.x * sinA * h,
        pos.y + N.y * cosA * w + B.y * sinA * h,
        pos.z + N.z * cosA * w + B.z * sinA * h
      );

      // Ellipse outward normal in Frenet frame
      const nx = N.x * halfH * cosA + B.x * halfW * sinA;
      const ny = N.y * halfH * cosA + B.y * halfW * sinA;
      const nz = N.z * halfH * cosA + B.z * halfW * sinA;
      const nl = Math.sqrt(nx * nx + ny * ny + nz * nz) || 1;
      normals.push(nx / nl, ny / nl, nz / nl);
    }
  }

  for (let i = 0; i < segments; i++) {
    for (let j = 0; j < profileN; j++) {
      const a = i * ringSize + j;
      const b = a + ringSize;
      const c = a + 1;
      const d = b + 1;
      indices.push(a, b, c, c, b, d);
    }
  }

  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
  geo.setAttribute('normal', new THREE.Float32BufferAttribute(normals, 3));
  geo.setIndex(indices);
  return geo;
}

/** Strand taper: widens from 1x to 1.8x at 60%, then tapers to 0 (arrowhead). */
export function strandTaperFn(t) {
  if (t < 0.6) return 1 + 0.8 * (t / 0.6);
  return 1.8 * (1 - (t - 0.6) / 0.4);
}

/**
 * Build a mesh for one secondary-structure element.
 * @param {THREE.Vector3[]} points - Backbone control points
 * @param {'H'|'E'|'L'} ssType - Secondary structure type
 * @param {THREE.Color} color - Base color
 * @param {number} [opacity] - If < 1, enables transparency
 */
export function buildSSElement(points, ssType, color, opacity) {
  let halfW, halfH, profileN, taperFn;
  if (ssType === 'H') { halfW = 0.55; halfH = 0.12; profileN = 8; taperFn = null; }
  else if (ssType === 'E') { halfW = 0.50; halfH = 0.12; profileN = 6; taperFn = strandTaperFn; }
  else { halfW = 0.10; halfH = 0.10; profileN = 6; taperFn = null; }

  const geo = sweepRibbon(points, halfW, halfH, profileN, taperFn);
  const matOpts = {
    color,
    roughness: 0.55,
    metalness: 0.05,
    clearcoat: 0.15,
    side: THREE.DoubleSide,
    emissive: color,
    emissiveIntensity: 0.10,
  };
  if (opacity != null && opacity < 1) {
    matOpts.transparent = true;
    matOpts.opacity = opacity;
    matOpts.depthWrite = false;
    matOpts.emissiveIntensity = 0.08;
  }
  return new THREE.Mesh(geo, new THREE.MeshPhysicalMaterial(matOpts));
}

/** Group backbone points into contiguous secondary-structure runs. */
export function segmentBySSType(backbone) {
  if (!backbone.length) return [];
  const segments = [];
  let startIdx = 0;
  for (let i = 1; i < backbone.length; i++) {
    if (backbone[i].ss !== backbone[startIdx].ss) {
      segments.push({ type: backbone[startIdx].ss, startIdx, endIdx: i - 1 });
      startIdx = i;
    }
  }
  segments.push({ type: backbone[startIdx].ss, startIdx, endIdx: backbone.length - 1 });
  return segments;
}
