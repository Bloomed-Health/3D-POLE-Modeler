/* global molstar */
// ==========================================================================
// Shared Mol* viewer initialization factory.
// Creates a Mol* Viewer instance with common options (dark background,
// hidden controls) and runs initial resize. Each page provides overrides.
// ==========================================================================

const BASE_OPTIONS = {
  layoutShowControls: false,
  layoutShowRemoteState: false,
  layoutShowSequence: false,
  layoutShowLog: false,
  layoutShowLeftPanel: false,
  viewportShowExpand: false,
  viewportShowSelectionMode: false,
  viewportShowAnimation: false,
  canvas3d: {
    attribs: { preserveDrawingBuffer: true },
  },
};

const DARK_BG = 0x0E141E;

/**
 * Create a Mol* Viewer with shared defaults.
 * @param {string} containerId  DOM id for the viewer container
 * @param {object} [options]    Mol* Viewer options merged over BASE_OPTIONS
 * @returns {Promise<object>}   The created molstar.Viewer instance
 */
export async function createMolstarViewer(containerId, options = {}) {
  const viewer = await molstar.Viewer.create(containerId, {
    ...BASE_OPTIONS,
    ...options,
  });

  // Dark background
  try {
    viewer.plugin.canvas3d?.setProps({
      renderer: { backgroundColor: DARK_BG },
    });
  } catch (e) { /* ignore if props differ between versions */ }

  // Force the plugin to detect its container size
  try { viewer.handleResize(); } catch (e) { /* ignore */ }

  return viewer;
}
