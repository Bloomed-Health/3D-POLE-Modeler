/**
 * Download a Three.js canvas as PNG screenshot.
 */

export function downloadScreenshot(renderer, scene, camera, prefix = 'POLE') {
  renderer.render(scene, camera);
  const dataURL = renderer.domElement.toDataURL('image/png');
  const a = document.createElement('a');
  a.href = dataURL;
  a.download = `${prefix}-${Date.now()}.png`;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
}
