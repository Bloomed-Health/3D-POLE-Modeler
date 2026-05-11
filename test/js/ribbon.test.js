import { describe, it, expect } from 'vitest';
import { strandTaperFn, segmentBySSType, sweepRibbon, buildSSElement } from '../../js/ribbon.js';
import * as THREE from 'three';

describe('strandTaperFn', () => {
  it('returns 1.0 at t=0', () => {
    expect(strandTaperFn(0)).toBeCloseTo(1.0);
  });

  it('returns 1.8 at t=0.6', () => {
    expect(strandTaperFn(0.6)).toBeCloseTo(1.8);
  });

  it('returns 0.0 at t=1.0', () => {
    expect(strandTaperFn(1.0)).toBeCloseTo(0.0);
  });

  it('is monotonically increasing from 0 to 0.6', () => {
    let prev = strandTaperFn(0);
    for (let t = 0.1; t <= 0.6; t += 0.1) {
      const val = strandTaperFn(t);
      expect(val).toBeGreaterThanOrEqual(prev);
      prev = val;
    }
  });

  it('is monotonically decreasing from 0.6 to 1.0', () => {
    let prev = strandTaperFn(0.6);
    for (let t = 0.7; t <= 1.0; t += 0.1) {
      const val = strandTaperFn(t);
      expect(val).toBeLessThanOrEqual(prev);
      prev = val;
    }
  });
});

describe('segmentBySSType', () => {
  it('returns empty array for empty input', () => {
    expect(segmentBySSType([])).toEqual([]);
  });

  it('returns single segment for uniform type', () => {
    const backbone = [
      { ss: 'H', x: 0, y: 0, z: 0 },
      { ss: 'H', x: 1, y: 0, z: 0 },
      { ss: 'H', x: 2, y: 0, z: 0 },
    ];
    const result = segmentBySSType(backbone);
    expect(result).toHaveLength(1);
    expect(result[0]).toEqual({ type: 'H', startIdx: 0, endIdx: 2 });
  });

  it('segments alternating types correctly', () => {
    const backbone = [
      { ss: 'H' }, { ss: 'H' },
      { ss: 'L' },
      { ss: 'E' }, { ss: 'E' }, { ss: 'E' },
    ];
    const result = segmentBySSType(backbone);
    expect(result).toHaveLength(3);
    expect(result[0]).toEqual({ type: 'H', startIdx: 0, endIdx: 1 });
    expect(result[1]).toEqual({ type: 'L', startIdx: 2, endIdx: 2 });
    expect(result[2]).toEqual({ type: 'E', startIdx: 3, endIdx: 5 });
  });

  it('handles single element input', () => {
    const backbone = [{ ss: 'L' }];
    const result = segmentBySSType(backbone);
    expect(result).toHaveLength(1);
    expect(result[0]).toEqual({ type: 'L', startIdx: 0, endIdx: 0 });
  });
});

describe('sweepRibbon', () => {
  it('returns empty geometry for fewer than 2 points', () => {
    const geo = sweepRibbon([new THREE.Vector3(0, 0, 0)], 0.5, 0.1, 6, null);
    expect(geo).toBeInstanceOf(THREE.BufferGeometry);
    expect(geo.getAttribute('position')).toBeUndefined();
  });

  it('returns valid BufferGeometry for 4+ points', () => {
    const points = [
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(1, 0, 0),
      new THREE.Vector3(2, 1, 0),
      new THREE.Vector3(3, 1, 1),
    ];
    const geo = sweepRibbon(points, 0.5, 0.12, 6, null);
    expect(geo).toBeInstanceOf(THREE.BufferGeometry);
    expect(geo.getAttribute('position')).not.toBeNull();
    expect(geo.getAttribute('normal')).not.toBeNull();
    expect(geo.getIndex()).not.toBeNull();
    expect(geo.getAttribute('position').count).toBeGreaterThan(0);
  });

  it('applies taper function when provided', () => {
    const points = [
      new THREE.Vector3(0, 0, 0),
      new THREE.Vector3(1, 0, 0),
      new THREE.Vector3(2, 0, 0),
      new THREE.Vector3(3, 0, 0),
    ];
    const geoNoTaper = sweepRibbon(points, 0.5, 0.12, 6, null);
    const geoTaper = sweepRibbon(points, 0.5, 0.12, 6, strandTaperFn);
    // Both should produce valid geometry but with different vertex positions
    expect(geoNoTaper.getAttribute('position').count).toBe(
      geoTaper.getAttribute('position').count
    );
    // Vertex data should differ due to taper
    const posNoTaper = geoNoTaper.getAttribute('position').array;
    const posTaper = geoTaper.getAttribute('position').array;
    let differs = false;
    for (let i = 0; i < posNoTaper.length; i++) {
      if (Math.abs(posNoTaper[i] - posTaper[i]) > 1e-6) {
        differs = true;
        break;
      }
    }
    expect(differs).toBe(true);
  });
});

describe('buildSSElement', () => {
  const points = [
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(1, 0.5, 0),
    new THREE.Vector3(2, 0, 0.5),
    new THREE.Vector3(3, -0.5, 0),
  ];
  const color = new THREE.Color(0x4488aa);

  it('returns a Mesh instance', () => {
    const mesh = buildSSElement(points, 'H', color, 1.0);
    expect(mesh).toBeInstanceOf(THREE.Mesh);
  });

  it('has correct material properties for helix', () => {
    const mesh = buildSSElement(points, 'H', color, 1.0);
    expect(mesh.material.roughness).toBeCloseTo(0.55);
    expect(mesh.material.metalness).toBeCloseTo(0.05);
    expect(mesh.material.side).toBe(THREE.DoubleSide);
  });

  it('enables transparency when opacity < 1', () => {
    const mesh = buildSSElement(points, 'E', color, 0.5);
    expect(mesh.material.transparent).toBe(true);
    expect(mesh.material.opacity).toBeCloseTo(0.5);
    expect(mesh.material.depthWrite).toBe(false);
  });

  it('does not enable transparency when opacity is 1', () => {
    const mesh = buildSSElement(points, 'L', color, 1.0);
    expect(mesh.material.transparent).toBeFalsy();
  });
});
