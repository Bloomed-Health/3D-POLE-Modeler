import { describe, it, expect, vi, beforeEach } from 'vitest';
import { loadPipelineData } from '../../js/pole-viewer.js';

describe('loadPipelineData', () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it('returns null when fetch fails', async () => {
    vi.stubGlobal('fetch', vi.fn().mockRejectedValue(new Error('Network error')));
    const result = await loadPipelineData('./data');
    expect(result).toBeNull();
  });

  it('returns null when response is not ok', async () => {
    vi.stubGlobal('fetch', vi.fn().mockResolvedValue({
      ok: false,
      status: 404,
      json: () => Promise.resolve({}),
    }));
    const result = await loadPipelineData('./data');
    expect(result).toBeNull();
  });
});
