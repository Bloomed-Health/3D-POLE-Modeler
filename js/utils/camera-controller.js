/**
 * Shared orbit/pan/zoom camera controls for Three.js viewers.
 *
 * The viewer must expose: cameraTheta, cameraPhi, cameraRadius,
 * cameraTarget (THREE.Vector3), _updateCamera(), isDragging, prevMouse.
 *
 * @param {object} viewer
 * @param {HTMLCanvasElement} canvas
 * @param {object} opts
 * @param {number}   opts.minRadius     - min zoom (default 20)
 * @param {number}   opts.maxRadius     - max zoom (default 200)
 * @param {Function} opts.onPointerDown - called on pointerdown; return false to suppress drag
 * @param {Function} opts.onLongPress   - called after 500ms single-touch hold
 * @param {Function} opts.onTouchEnd    - called on touchend
 * @returns {{ destroy: Function }}
 */
export function attachCameraControls(viewer, canvas, opts = {}) {
  const minR = opts.minRadius ?? 20;
  const maxR = opts.maxRadius ?? 200;

  let shiftDown = false;
  let lastPinchDist = 0;
  let lastTwoFingerMid = null;
  let longPressTimer = null;

  // --- Pointer events (mouse / stylus) ---

  function onPointerDown(e) {
    if (opts.onPointerDown && opts.onPointerDown(e) === false) return;
    viewer.isDragging = true;
    shiftDown = e.shiftKey || e.button === 1;
    viewer.prevMouse = { x: e.clientX, y: e.clientY };
  }

  function onPointerUp() {
    viewer.isDragging = false;
    shiftDown = false;
  }

  function onPointerMove(e) {
    if (!viewer.isDragging) return;
    const dx = e.clientX - viewer.prevMouse.x;
    const dy = e.clientY - viewer.prevMouse.y;
    if (shiftDown || e.shiftKey) {
      const panScale = viewer.cameraRadius * 0.002;
      const rx = Math.cos(viewer.cameraTheta);
      const rz = -Math.sin(viewer.cameraTheta);
      viewer.cameraTarget.x -= dx * panScale * rx;
      viewer.cameraTarget.z -= dx * panScale * rz;
      viewer.cameraTarget.y += dy * panScale;
    } else {
      viewer.cameraTheta -= dx * 0.005;
      viewer.cameraPhi = Math.max(-1.2, Math.min(1.2, viewer.cameraPhi + dy * 0.005));
    }
    viewer.prevMouse = { x: e.clientX, y: e.clientY };
    viewer._updateCamera();
  }

  function onWheel(e) {
    e.preventDefault();
    viewer.cameraRadius = Math.max(minR, Math.min(maxR, viewer.cameraRadius + e.deltaY * 0.04));
    viewer._updateCamera();
  }

  // --- Touch events ---

  function onTouchStart(e) {
    e.preventDefault();
    clearTimeout(longPressTimer);
    if (e.touches.length === 1) {
      viewer.isDragging = true;
      viewer.prevMouse = { x: e.touches[0].clientX, y: e.touches[0].clientY };
      if (opts.onLongPress) {
        longPressTimer = setTimeout(() => {
          opts.onLongPress({ clientX: e.touches[0].clientX, clientY: e.touches[0].clientY });
        }, 500);
      }
    } else if (e.touches.length === 2) {
      viewer.isDragging = false;
      lastPinchDist = Math.hypot(
        e.touches[0].clientX - e.touches[1].clientX,
        e.touches[0].clientY - e.touches[1].clientY);
    }
  }

  function onTouchMove(e) {
    e.preventDefault();
    clearTimeout(longPressTimer);
    if (e.touches.length === 1 && viewer.isDragging) {
      const dx = e.touches[0].clientX - viewer.prevMouse.x;
      const dy = e.touches[0].clientY - viewer.prevMouse.y;
      viewer.cameraTheta -= dx * 0.005;
      viewer.cameraPhi = Math.max(-1.2, Math.min(1.2, viewer.cameraPhi + dy * 0.005));
      viewer.prevMouse = { x: e.touches[0].clientX, y: e.touches[0].clientY };
      viewer._updateCamera();
    } else if (e.touches.length === 2) {
      const dist = Math.hypot(
        e.touches[0].clientX - e.touches[1].clientX,
        e.touches[0].clientY - e.touches[1].clientY);
      const delta = lastPinchDist - dist;
      viewer.cameraRadius = Math.max(minR, Math.min(maxR, viewer.cameraRadius + delta * 0.1));
      lastPinchDist = dist;

      // Two-finger pan
      const midX = (e.touches[0].clientX + e.touches[1].clientX) / 2;
      const midY = (e.touches[0].clientY + e.touches[1].clientY) / 2;
      if (lastTwoFingerMid) {
        const panDx = midX - lastTwoFingerMid.x;
        const panDy = midY - lastTwoFingerMid.y;
        const panScale = viewer.cameraRadius * 0.002;
        const rx = Math.cos(viewer.cameraTheta);
        const rz = -Math.sin(viewer.cameraTheta);
        viewer.cameraTarget.x -= panDx * panScale * rx;
        viewer.cameraTarget.z -= panDx * panScale * rz;
        viewer.cameraTarget.y += panDy * panScale;
      }
      lastTwoFingerMid = { x: midX, y: midY };
      viewer._updateCamera();
    }
  }

  function onTouchEnd() {
    clearTimeout(longPressTimer);
    viewer.isDragging = false;
    lastTwoFingerMid = null;
    if (opts.onTouchEnd) opts.onTouchEnd();
  }

  // Attach
  canvas.addEventListener('pointerdown', onPointerDown);
  window.addEventListener('pointerup', onPointerUp);
  window.addEventListener('pointermove', onPointerMove);
  canvas.addEventListener('wheel', onWheel, { passive: false });
  canvas.addEventListener('touchstart', onTouchStart, { passive: false });
  canvas.addEventListener('touchmove', onTouchMove, { passive: false });
  canvas.addEventListener('touchend', onTouchEnd);

  return {
    destroy() {
      canvas.removeEventListener('pointerdown', onPointerDown);
      window.removeEventListener('pointerup', onPointerUp);
      window.removeEventListener('pointermove', onPointerMove);
      canvas.removeEventListener('wheel', onWheel);
      canvas.removeEventListener('touchstart', onTouchStart);
      canvas.removeEventListener('touchmove', onTouchMove);
      canvas.removeEventListener('touchend', onTouchEnd);
    },
  };
}
