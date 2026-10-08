/**
 * Simple FPS counter for Three.js render loops.
 * Replaces duplicated fpsFrames/fpsTime/fpsCallback pattern.
 */

export class FPSCounter {
  constructor() {
    this.frames = 0;
    this.time = 0;
    this.callback = null;
  }

  tick(t) {
    this.frames++;
    if (t - this.time >= 1) {
      if (this.callback) this.callback(this.frames);
      this.frames = 0;
      this.time = t;
    }
  }

  setCallback(cb) { this.callback = cb; }
}
