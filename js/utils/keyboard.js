/**
 * Shared keyboard shortcut dispatcher.
 * Filters out events from input/textarea fields.
 *
 * @param {Object<string, {action: Function, noModifier?: boolean}>} actions
 *   Keys are lowercase key names. Each value has an `action(e)` function
 *   and optional `noModifier` flag to suppress when Ctrl/Cmd is held.
 */
export function initKeyboardShortcuts(actions) {
  document.addEventListener('keydown', e => {
    if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA') return;
    const key = e.key.toLowerCase();
    const handler = actions[key];
    if (handler) {
      if (handler.noModifier && (e.ctrlKey || e.metaKey)) return;
      e.preventDefault();
      handler.action(e);
    }
  });
}
