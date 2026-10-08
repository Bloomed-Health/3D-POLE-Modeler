/**
 * Shared camera animation utilities.
 * Drop-in replacement for duplicated easing logic in pole-viewer.js and mutation-viewer.js.
 */

export function easeInOutQuad(p) {
  return p < 0.5 ? 2 * p * p : 1 - Math.pow(-2 * p + 2, 2) / 2;
}

/**
 * Smoothly animate camera parameters on a viewer that exposes:
 *   cameraTheta, cameraPhi, cameraRadius, cameraTarget (THREE.Vector3), _updateCamera()
 *
 * @param {object} viewer
 * @param {object} endState - any subset of { th, ph, r, tx, ty, tz }
 * @param {number} duration - ms (default 800)
 */
export function animateCamera(viewer, endState, duration = 800) {
  const s = {
    th: viewer.cameraTheta, ph: viewer.cameraPhi, r: viewer.cameraRadius,
    tx: viewer.cameraTarget.x, ty: viewer.cameraTarget.y, tz: viewer.cameraTarget.z,
  };
  const e = { th: s.th, ph: s.ph, r: s.r, tx: s.tx, ty: s.ty, tz: s.tz, ...endState };
  const t0 = performance.now();
  const tick = () => {
    const p = Math.min(1, (performance.now() - t0) / duration);
    const ease = easeInOutQuad(p);
    viewer.cameraTheta = s.th + (e.th - s.th) * ease;
    viewer.cameraPhi = s.ph + (e.ph - s.ph) * ease;
    viewer.cameraRadius = s.r + (e.r - s.r) * ease;
    viewer.cameraTarget.set(
      s.tx + (e.tx - s.tx) * ease,
      s.ty + (e.ty - s.ty) * ease,
      s.tz + (e.tz - s.tz) * ease);
    viewer._updateCamera();
    if (p < 1) requestAnimationFrame(tick);
  };
  tick();
}
