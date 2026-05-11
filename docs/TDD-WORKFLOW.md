# TDD Workflow Guide

## Red-Green-Refactor Cycle

This project follows Test-Driven Development (TDD) for all new JavaScript functionality.

### The Cycle

1. **Red** — Write a failing test that describes the desired behavior
2. **Green** — Write the minimum code to make the test pass
3. **Refactor** — Clean up the implementation while keeping tests green

### Running Tests

```bash
# Run all tests once
npm test

# Watch mode (re-runs on file changes)
npm run test:watch

# Generate coverage report
npm run test:coverage
```

## Project-Specific Examples

### Example 1: Adding a new geometry helper

**Red** — Write the test first:
```javascript
// test/js/ribbon.test.js
import { newHelper } from '../../js/ribbon.js';

describe('newHelper', () => {
  it('returns expected geometry for valid input', () => {
    const result = newHelper(inputPoints);
    expect(result.getAttribute('position').count).toBeGreaterThan(0);
  });
});
```

**Green** — Implement in `js/ribbon.js`:
```javascript
export function newHelper(points) {
  // Minimum implementation to pass the test
}
```

**Refactor** — Optimize, extract constants, improve naming.

### Example 2: Adding a data transformation

**Red** — Describe the contract:
```javascript
// test/js/pole-viewer.test.js
describe('transformScores', () => {
  it('normalizes scores to 0-1 range', () => {
    const raw = [0, 50, 100];
    const result = transformScores(raw);
    expect(result).toEqual([0, 0.5, 1.0]);
  });
});
```

## Guidelines for Testable Code

### Extract Pure Functions

Large modules like `pole-viewer.js` contain complex logic. To make code testable:

1. **Identify pure computations** — Functions that take input and return output without side effects
2. **Extract to separate modules** — Move geometry math, data transforms, and scoring logic into importable functions
3. **Keep DOM/Three.js setup in the main module** — Side-effectful code stays in the viewer class

### What to Test

| Category | Test? | Example |
|----------|-------|---------|
| Pure math/geometry | Yes | `strandTaperFn`, `sweepRibbon` |
| Data transformations | Yes | Score normalization, domain lookup |
| API contracts | Yes | `loadPipelineData` return shape |
| DOM manipulation | Sparingly | Use jsdom if needed |
| Three.js rendering | No | Visual tests are fragile |

### Mocking External Dependencies

Use Vitest's built-in mocking for `fetch`, DOM APIs, etc.:

```javascript
import { vi } from 'vitest';

vi.stubGlobal('fetch', vi.fn().mockResolvedValue({
  ok: true,
  json: () => Promise.resolve({ data: [] }),
}));
```

## Coverage

Coverage reports are generated in `coverage/`. Aim for:
- **100%** on pure utility functions (ribbon.js)
- **High** on data loading and transformation logic
- **Reasonable** on integration paths

Coverage reports are excluded from git via `.gitignore`.
