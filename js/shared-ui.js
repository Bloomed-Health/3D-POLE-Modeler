// ==========================================================================
// Shared UI utilities used by all POLE viewer pages.
// Sidebar toggle, loading progress, provenance footer.
// ==========================================================================

/**
 * Initialize the loading progress bar.
 * @returns {function(number, string=): void} setProgress(pct, text)
 */
export function initProgress() {
  const progressFill = document.getElementById('progress-fill');
  const loadingText = document.getElementById('loading-text');
  return function setProgress(pct, text) {
    progressFill.style.width = pct + '%';
    if (text) loadingText.textContent = text;
  };
}

/**
 * Wire up mobile sidebar toggle, close buttons, backdrop, and tab buttons.
 * Expects DOM elements: sidebar-toggle, sidebar-left, sidebar-right,
 * sidebar-backdrop, sidebar-close-left, sidebar-close-right,
 * and .mobile-sidebar-tabs button[data-sidebar].
 */
export function initSidebarToggle() {
  const sidebarToggle = document.getElementById('sidebar-toggle');
  const sidebarLeft = document.getElementById('sidebar-left');
  const sidebarRight = document.getElementById('sidebar-right');
  const sidebarBackdrop = document.getElementById('sidebar-backdrop');

  function openSidebar(which) {
    sidebarLeft.classList.toggle('open', which === 'left');
    sidebarRight.classList.toggle('open', which === 'right');
    sidebarBackdrop.classList.add('open');
    document.querySelectorAll('.mobile-sidebar-tabs button').forEach(btn => {
      btn.classList.toggle('active', btn.dataset.sidebar === which);
    });
  }
  function closeSidebars() {
    sidebarLeft.classList.remove('open');
    sidebarRight.classList.remove('open');
    sidebarBackdrop.classList.remove('open');
  }

  sidebarToggle.addEventListener('click', () => {
    if (sidebarLeft.classList.contains('open') || sidebarRight.classList.contains('open')) closeSidebars();
    else openSidebar('left');
  });
  document.getElementById('sidebar-close-left').addEventListener('click', closeSidebars);
  document.getElementById('sidebar-close-right').addEventListener('click', closeSidebars);
  sidebarBackdrop.addEventListener('click', closeSidebars);
  document.querySelectorAll('.mobile-sidebar-tabs button').forEach(btn => {
    btn.addEventListener('click', () => openSidebar(btn.dataset.sidebar));
  });
}

/**
 * Load canonical data and populate the [data-pole-provenance] footer element.
 * Degrades silently if the canonical JSON is unavailable.
 */
export async function initProvenanceFooter() {
  try {
    const { loadPoleData, provenanceFooter } = await import('./data/poleData.js');
    const data = await loadPoleData();
    const el = document.querySelector('[data-pole-provenance]');
    if (el) el.textContent = provenanceFooter(data);
  } catch (e) { /* degrade silently if canonical JSON unavailable */ }
}
