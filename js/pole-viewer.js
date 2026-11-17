// ==========================================================================
// POLE p261 -3D Structural Viewer
// Cartoon Ribbon Representation
//
// The protein fold is rendered as cartoon ribbons -helices as wide flat
// spirals, strands as flat arrows with arrowheads, loops as thin tubes.
// DNA template/primer rendered as backbone tubes with base-pair rungs.
//
// References:
//   UniProt Q07864 (POLE_HUMAN)
//   PDB: 9F6D / 9F6E / 9F6F (Roske & Yeeles 2024), 9B8S (He et al. 2024),
//        4M8O (Hogg et al. 2014)
//   Roske & Yeeles, Nat Struct Mol Biol 31:1921 (2024) — human Pol ε structures
//   He, Wang, Yao, O'Donnell & Li, Nat Commun 15:7847 (2024) — human Pol ε–PCNA
//   Mur et al., Genome Med 15:85 (2023) — POLE/POLD1 ACMG classification
//   Yuan, Georgescu, Schauer, O'Donnell & Li, Nat Commun 11:3156 (2020) — yeast Pol ε
//   Robinson, Coorens et al., Nat Genet 53:1434 (2021) — normal-tissue mutation rates
//   Hogg et al., Nat Struct Mol Biol 21:49 (2014) — S. cerevisiae Pol2 NTD (4M8O)
//
// NOTE: This viewer renders a procedural cartoon-ribbon schematic, not atomic
// coordinates from PDB depositions. Color modes (B-factor, conservation, charge,
// pathogenicity) use approximate domain-level or synthetic per-residue values.
// For coordinate-derived analysis, load PDB 9F6D into Mol* or PyMOL.
// ==========================================================================

import * as THREE from 'three';
import { buildSSElement, segmentBySSType } from './ribbon.js';
import { DOMAINS as BASE_DOMAINS, DOMAIN_ORDER, MUTATIONS as BASE_MUTATIONS } from './data/domainConfig.js';
import { PAL, COLOR_SCHEMES } from './pole-palette.js';
import { SCALE, generateDomainBackbone, generateLinker } from './pole-geometry.js';
import { animateCamera } from './utils/animate.js';
import { FPSCounter } from './utils/fps-counter.js';
import { downloadScreenshot } from './utils/screenshot.js';
import { attachCameraControls } from './utils/camera-controller.js';


// ==================== DOMAIN ARCHITECTURE ====================

// Schematic-specific domain extensions (center, ss, header are viewer-only)
const SCHEMATIC_EXTENSIONS = {
  ntd: { center:[-18,8,-5], color:'ntd',
    header:'// NTD (1-222) // Structural scaffold',
    ss:[ {t:'H',n:14},{t:'L',n:6},{t:'E',n:7},{t:'L',n:4},{t:'E',n:6},{t:'L',n:5},{t:'H',n:18},{t:'L',n:5},{t:'E',n:8},{t:'L',n:4},{t:'E',n:7},{t:'L',n:6},{t:'H',n:12},{t:'L',n:8},{t:'H',n:10},{t:'L',n:5},{t:'E',n:6},{t:'L',n:4},{t:'H',n:16},{t:'L',n:7},{t:'H',n:10},{t:'L',n:8},{t:'E',n:5},{t:'L',n:6},{t:'H',n:8},{t:'L',n:9} ] },
  exo: { center:[-10,2,7], color:'exo',
    header:"// EXO (223-517) // 3'\u21925' proofreading / D275/E277(ExoI) D368(ExoII,catalytic) D462(ExoIII)",
    ss:[ {t:'E',n:6},{t:'L',n:4},{t:'H',n:14},{t:'L',n:5},{t:'E',n:7},{t:'L',n:3},{t:'E',n:6},{t:'H',n:12},{t:'L',n:6},{t:'H',n:16},{t:'L',n:4},{t:'E',n:8},{t:'L',n:5},{t:'H',n:10},{t:'L',n:6},{t:'E',n:7},{t:'L',n:3},{t:'H',n:14},{t:'L',n:5},{t:'H',n:8},{t:'L',n:6},{t:'E',n:6},{t:'L',n:4},{t:'H',n:10},{t:'L',n:4} ] },
  palm: { center:[0,-1,0], color:'palm',
    header:'// PALM (518-870) // Catalytic core / motifs A(D626) C(D860/D862)',
    ss:[ {t:'E',n:8},{t:'L',n:3},{t:'E',n:7},{t:'L',n:5},{t:'H',n:12},{t:'L',n:4},{t:'E',n:9},{t:'L',n:3},{t:'E',n:7},{t:'L',n:6},{t:'H',n:10},{t:'L',n:4},{t:'E',n:8},{t:'L',n:5},{t:'E',n:6},{t:'L',n:3},{t:'H',n:14},{t:'L',n:5},{t:'E',n:7},{t:'L',n:4},{t:'H',n:8},{t:'L',n:6},{t:'E',n:9},{t:'L',n:3},{t:'E',n:7},{t:'L',n:5},{t:'H',n:12},{t:'L',n:4} ] },
  pdomain: { center:[4,1,-3], color:'pdomain',
    header:'// P-DOM (640-760) // [4Fe-4S] cluster / PCNA contact',
    ss:[ {t:'H',n:12},{t:'L',n:4},{t:'H',n:10},{t:'L',n:5},{t:'E',n:6},{t:'L',n:3},{t:'H',n:8},{t:'L',n:4} ] },
  fingers: { center:[7,5,4], color:'fingers',
    header:'// FNG (870-1020) // O-helix / Watson-Crick fidelity gate',
    ss:[ {t:'H',n:20},{t:'L',n:4},{t:'H',n:16},{t:'L',n:6},{t:'H',n:22},{t:'L',n:5},{t:'H',n:14},{t:'L',n:4},{t:'H',n:12},{t:'L',n:5},{t:'H',n:10},{t:'L',n:6},{t:'H',n:8},{t:'L',n:4} ] },
  thumb: { center:[3,-9,7], color:'thumb',
    header:'// THM (1020-1189) // Helical bundle / grips duplex DNA',
    ss:[ {t:'H',n:18},{t:'L',n:4},{t:'H',n:22},{t:'L',n:5},{t:'H',n:16},{t:'L',n:6},{t:'H',n:20},{t:'L',n:4},{t:'H',n:14},{t:'L',n:5},{t:'H',n:12},{t:'L',n:6},{t:'H',n:10},{t:'L',n:4} ] },
  inactpol: { center:[10,-4,1], color:'inactpol',
    header:'// INACT (1190-1900) // Non-catalytic Pol fold / subunit scaffold',
    ss:[ {t:'H',n:14},{t:'L',n:6},{t:'E',n:7},{t:'L',n:4},{t:'H',n:10},{t:'L',n:5},{t:'E',n:6},{t:'L',n:4},{t:'H',n:12},{t:'L',n:6},{t:'E',n:8},{t:'L',n:3},{t:'H',n:8},{t:'L',n:5},{t:'E',n:7},{t:'L',n:4},{t:'H',n:10},{t:'L',n:6},{t:'E',n:6},{t:'L',n:4},{t:'H',n:14},{t:'L',n:5},{t:'E',n:5},{t:'L',n:4},{t:'H',n:8},{t:'L',n:6} ] },
  ctd: { center:[16,1,-5], color:'ctd',
    header:'// CTD (1900-2286) // CysA/CysB Zn-finger / POLE2/3/4 anchor',
    ss:[ {t:'H',n:10},{t:'L',n:6},{t:'E',n:5},{t:'L',n:4},{t:'E',n:6},{t:'L',n:3},{t:'H',n:14},{t:'L',n:5},{t:'E',n:7},{t:'L',n:4},{t:'E',n:5},{t:'L',n:6},{t:'H',n:12},{t:'L',n:4},{t:'E',n:6},{t:'L',n:5},{t:'H',n:10},{t:'L',n:8},{t:'E',n:5},{t:'L',n:4},{t:'H',n:8},{t:'L',n:6},{t:'E',n:7},{t:'L',n:3},{t:'H',n:12},{t:'L',n:5},{t:'E',n:6},{t:'L',n:4},{t:'H',n:10},{t:'L',n:7},{t:'E',n:5},{t:'L',n:5},{t:'H',n:14},{t:'L',n:6} ] },
};

// Merge shared domain definitions with schematic-specific extensions
const DOMAINS = {};
for (const key of DOMAIN_ORDER) {
  DOMAINS[key] = { ...BASE_DOMAINS[key], ...SCHEMATIC_EXTENSIONS[key] };
}

// Extend shared mutations with schematic-specific offsets
const MUTATION_OFFSETS = {
  P286R: [2.5,1.2,2.0], V411L: [-2.8,-0.5,1.8], S297F: [0.5,2.8,-1.2],
  L424V: [-1.5,2.0,2.5], D287E: [3.0,-1.0,1.0], P436R: [-2.0,1.5,-1.5],
  M444K: [1.0,-2.0,2.5], S459F: [-3.0,0.5,-0.5], F367S: [2.0,0.8,-2.0],
  A456P: [-1.0,-1.5,1.5],
};
let MUTATIONS = BASE_MUTATIONS.map(m => ({ ...m, offset: MUTATION_OFFSETS[m.id] || [0,0,0] }));


// ==================== VIEWER CLASS ====================

export class POLEViewer {
  constructor(containerEl, tooltipEl) {
    this.container = containerEl;
    this.tooltipEl = tooltipEl;
    this.autoRotate = !window.matchMedia('(prefers-reduced-motion: reduce)').matches;
    this.isDragging = false;
    this.prevMouse = { x:0, y:0 };
    this.cameraTheta = -0.3;
    this.cameraPhi = 0.15;
    this.cameraRadius = 155;
    this.hoverTarget = null;
    this.groups = {};
    this.domainMeshes = {};
    this.mutMarkers = {};
    this.activeSitePos = new THREE.Vector3();
    this.cameraTarget = new THREE.Vector3(0, 0, 0);

    // New state for features
    this.backboneCache = {};
    this.colorMode = 'default';
    this.colorScheme = 'default';
    this.measureMode = false;
    this.measurePoints = [];
    this.measureLine = null;
    this.measureLabel = null;
    this.exoPolT = 0;
    this.exoPolPlaying = false;
    this._fps = new FPSCounter();
    this.labelMinPriority = 2;
    this._labelSprites = [];

    this._init();
  }

  _init() {
    this._createScene();
    this._createLighting();
    this._buildProtein();
    this._buildDNA();
    this._buildActiveSite();
    this._buildZincMotifs();
    this._buildFeSCluster();
    this._buildExoSite();
    this._buildAccessorySubunits();
    this._buildPCNA();
    this._buildMutationMarkers();
    this._buildAnnotations();
    this._buildHydrogenBonds();
    this._buildBetaSheetHBonds();
    this._buildExoPolPath();
    this._buildParticles();
    this._setupControls();
    this._startRenderLoop();
  }

  _createScene() {
    this.scene = new THREE.Scene();
    this.scene.fog = new THREE.FogExp2(PAL.bg, 0.003);
    this.camera = new THREE.PerspectiveCamera(40,1,0.1,2000);
    this.renderer = new THREE.WebGLRenderer({ antialias:true, alpha:true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio,2));
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
    this.renderer.toneMappingExposure = 1.4;
    this.renderer.domElement.setAttribute('role', 'img');
    this.renderer.domElement.setAttribute('aria-label', '3D protein structure schematic');
    this.container.appendChild(this.renderer.domElement);
    const resize = () => {
      const w = this.container.clientWidth, h = this.container.clientHeight;
      this.renderer.setSize(w,h,false);
      this.camera.aspect = w/h;
      this.camera.updateProjectionMatrix();
    };
    new ResizeObserver(resize).observe(this.container);
    resize();
    this._updateCamera();
  }

  _createLighting() {
    this.scene.add(new THREE.AmbientLight(0x8090a0, 0.8));
    const key = new THREE.DirectionalLight(0xeee8dd, 1.3);
    key.position.set(18,25,20);
    this.scene.add(key);
    const fill = new THREE.DirectionalLight(0x6a7a8a, 0.4);
    fill.position.set(-10,-15,10);
    this.scene.add(fill);
    const rim = new THREE.DirectionalLight(0x7abba8, 0.5);
    rim.position.set(-20,8,-15);
    this.scene.add(rim);
    this.activeSiteLight = new THREE.PointLight(0xc45c8c, 0.6, 18);
    this.scene.add(this.activeSiteLight);
  }

  // ==================== PROTEIN: cartoon ribbon representation ====================

  _buildProtein() {
    this.groups.protein = new THREE.Group();
    this.groups.surface = new THREE.Group();
    this.groups.surface.visible = true;
    this.scene.add(this.groups.protein);
    this.groups.protein.add(this.groups.surface);

    const backbones = {};

    for (const key of DOMAIN_ORDER) {
      const dom = DOMAINS[key];
      const bb = generateDomainBackbone(key, DOMAINS);
      backbones[key] = bb;
      this.backboneCache[key] = bb;
      const domColor = PAL[dom.color];

      // Build ribbon meshes per SS segment
      const segments = segmentBySSType(bb);
      const ribbons = [];
      for (let s = 0; s < segments.length; s++) {
        const seg = segments[s];
        // 1-point overlap at boundaries for smooth joins
        const lo = Math.max(0, seg.startIdx - (s > 0 ? 1 : 0));
        const hi = Math.min(bb.length - 1, seg.endIdx + (s < segments.length - 1 ? 1 : 0));
        const pts = [];
        for (let i = lo; i <= hi; i++) pts.push(bb[i].pos);
        if (pts.length >= 2) {
          const mesh = buildSSElement(pts, seg.type, domColor);
          mesh.userData.ssType = seg.type;
          mesh.userData.domainKey = key;
          mesh.userData.startIdx = seg.startIdx;
          mesh.userData.endIdx = seg.endIdx;
          this.groups.protein.add(mesh);
          ribbons.push(mesh);
        }
      }

      // Raycasting proxy
      const bbox = new THREE.Box3();
      for (const r of bb) bbox.expandByPoint(r.pos);
      const center = new THREE.Vector3();
      bbox.getCenter(center);
      const size = new THREE.Vector3();
      bbox.getSize(size);
      const proxy = new THREE.Mesh(
        new THREE.SphereGeometry(Math.max(size.x,size.y,size.z)*0.5, 8, 8),
        new THREE.MeshBasicMaterial({ visible:false })
      );
      proxy.position.copy(center);
      proxy.userData = { domain:key, name:dom.name, detail:`res ${dom.range[0]}\u2013${dom.range[1]}` };
      this.groups.protein.add(proxy);
      this.domainMeshes[key] = { proxy, ribbons };
    }

    // Inter-domain linkers as thin tubes
    for (let i = 0; i < DOMAIN_ORDER.length - 1; i++) {
      const kA = DOMAIN_ORDER[i], kB = DOMAIN_ORDER[i+1];
      const bbA = backbones[kA], bbB = backbones[kB];
      const linkerN = Math.max(8, Math.floor((DOMAINS[kB].range[0] - DOMAINS[kA].range[1])/5));
      const linkerBB = generateLinker(bbA[bbA.length-1].pos, bbB[0].pos, linkerN);
      const linkerPts = linkerBB.map(p => p.pos);
      if (linkerPts.length >= 2) {
        const mesh = buildSSElement(linkerPts, 'L', PAL.linker);
        this.groups.protein.add(mesh);
      }
    }

    // Surface envelope (toggle)
    for (const key of DOMAIN_ORDER) {
      const c = DOMAINS[key].center;
      const col = PAL[DOMAINS[key].color];
      const mesh = new THREE.Mesh(
        new THREE.IcosahedronGeometry(6,3),
        new THREE.MeshPhysicalMaterial({
          color:col, transparent:true, opacity:0.09, roughness:0.9,
          side:THREE.DoubleSide, depthWrite:false })
      );
      mesh.position.set(c[0],c[1],c[2]);
      mesh.userData.surfaceDomain = key;
      this.groups.surface.add(mesh);
    }
  }

  // ==================== DNA: backbone tubes + base-pair rungs ====================

  _buildDNA() {
    this.groups.dna = new THREE.Group();
    this.scene.add(this.groups.dna);

    const axisPts = [];
    for (let i = 0; i < 80; i++) {
      const t = i/79;
      axisPts.push(new THREE.Vector3(
        -8+t*26, -9+Math.sin(t*Math.PI)*5+t*4,
        8-t*10+Math.sin(t*Math.PI*1.3)*2.5));
    }
    const axisCurve = new THREE.CatmullRomCurve3(axisPts);

    const helixR = 10.0*SCALE*0.5;
    const turnsPerUnit = 6;
    const segments = 400;
    const primerCutoff = 0.52;
    const grooveOffset = Math.PI*(154/180);
    const strand1=[], strand2=[];
    const basePairs = [];

    for (let i = 0; i <= segments; i++) {
      const t = i/segments;
      const center = axisCurve.getPoint(t);
      const tangent = axisCurve.getTangent(t).normalize();
      let up = Math.abs(tangent.y)<0.95 ? new THREE.Vector3(0,1,0) : new THREE.Vector3(1,0,0);
      const normal = new THREE.Vector3().crossVectors(tangent,up).normalize();
      const binormal = new THREE.Vector3().crossVectors(tangent,normal).normalize();
      const angle = t*turnsPerUnit*Math.PI*2;

      const o1 = normal.clone().multiplyScalar(Math.cos(angle)*helixR)
        .add(binormal.clone().multiplyScalar(Math.sin(angle)*helixR));
      const o2 = normal.clone().multiplyScalar(Math.cos(angle+grooveOffset)*helixR)
        .add(binormal.clone().multiplyScalar(Math.sin(angle+grooveOffset)*helixR));

      strand1.push(center.clone().add(o1));
      strand2.push(center.clone().add(o2));

      if (i % 14 === 0 && t < primerCutoff) {
        basePairs.push({ p1: center.clone().add(o1), p2: center.clone().add(o2) });
      }
    }

    // Template backbone tube
    if (strand1.length > 2) {
      const curve = new THREE.CatmullRomCurve3(strand1);
      const geo = new THREE.TubeGeometry(curve, 200, 0.12, 8, false);
      const mat = new THREE.MeshPhysicalMaterial({
        color: PAL.dnaTemplate, roughness:0.5, metalness:0.05, clearcoat:0.15
      });
      this.groups.dna.add(new THREE.Mesh(geo,mat));
    }
    // Primer backbone tube
    const primerPts = strand2.slice(0, Math.floor(strand2.length*primerCutoff));
    if (primerPts.length > 2) {
      const curve = new THREE.CatmullRomCurve3(primerPts);
      const geo = new THREE.TubeGeometry(curve, 120, 0.12, 8, false);
      const mat = new THREE.MeshPhysicalMaterial({
        color: PAL.dnaPrimer, roughness:0.5, metalness:0.05, clearcoat:0.15
      });
      this.groups.dna.add(new THREE.Mesh(geo,mat));
    }

    // Base-pair rungs (cylinders)
    for (const bp of basePairs) {
      const dir = bp.p2.clone().sub(bp.p1);
      const len = dir.length();
      const rung = new THREE.Mesh(
        new THREE.CylinderGeometry(0.04, 0.04, len, 6),
        new THREE.MeshPhysicalMaterial({ color: PAL.dnaBase, roughness:0.5, metalness:0.05 })
      );
      rung.position.copy(bp.p1.clone().add(bp.p2).multiplyScalar(0.5));
      rung.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0), dir.normalize());
      this.groups.dna.add(rung);
    }

    const pidx = Math.floor(segments*primerCutoff);
    this.activeSitePos.copy(strand2[pidx] || new THREE.Vector3());
  }

  // ==================== ACTIVE SITE: ball-and-stick + annotation ====================

  _buildActiveSite() {
    this.groups.metals = new THREE.Group();
    this.groups.dntp = new THREE.Group();
    this.scene.add(this.groups.metals);
    this.scene.add(this.groups.dntp);

    const asp = this.activeSitePos;
    const mgApos = asp.clone().add(new THREE.Vector3(-0.6, 0.30, 0.10));
    const mgBpos = asp.clone().add(new THREE.Vector3(1.30,-0.25, 0.50));

    for (const [pos, label] of [[mgApos,'A'],[mgBpos,'B']]) {
      const core = new THREE.Mesh(
        new THREE.SphereGeometry(0.35,16,16),
        new THREE.MeshPhysicalMaterial({
          color:0x60a830, roughness:0.2, metalness:0.7,
          emissive:0x60a830, emissiveIntensity:0.5 })
      );
      core.position.copy(pos);
      core.userData = { name:`Mg\u00b2\u207a (Metal ${label})`,
        detail: label==='A' ? 'Activates 3\u2032-OH \u00b7 D626' : 'Coordinates \u03b2,\u03b3-phosphates \u00b7 D860/D862' };
      this.groups.metals.add(core);
    }

    // Octahedral coordination: 6 ligand spheres per metal
    const octDirs = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]];
    const ligDist = 2.1 * SCALE;
    for (const pos of [mgApos,mgBpos]) {
      for (const d of octDirs) {
        const lPos = pos.clone().add(new THREE.Vector3(...d).multiplyScalar(ligDist));
        const lig = this._atom(0xc03030, 0.08);
        lig.position.copy(lPos);
        this.groups.metals.add(lig);
        this.groups.metals.add(this._dash(pos,lPos,0x60a830));
      }
    }
    this.groups.metals.add(this._dash(mgApos,mgBpos,0x60a830));

    // Mg-Mg distance label
    const mgMgMid = mgApos.clone().add(mgBpos).multiplyScalar(0.5);
    const mgMgDist = mgApos.distanceTo(mgBpos) / SCALE;
    const mgMgLabelPos = mgMgMid.clone().add(new THREE.Vector3(6, -8, 6));
    const mgMgLabel = this._spriteLabel(`Mg\u2013Mg: ${mgMgDist.toFixed(1)} \u00c5`, '#60a830', '#60a830', 220);
    mgMgLabel.scale.set(7, 2, 1);
    mgMgLabel.position.copy(mgMgLabelPos);
    mgMgLabel.userData.labelPriority = 1;
    this.groups.metals.add(mgMgLabel);
    this.groups.metals.add(this._leaderLine(mgMgLabelPos, mgMgMid, PAL.magnesium));

    const catRes = [
      { l:'D626', pos:asp.clone().add(new THREE.Vector3(-1.2,0.8,-0.5)), mg:mgApos },
      { l:'D860', pos:asp.clone().add(new THREE.Vector3(-0.8,-0.6,-0.8)), mg:mgBpos },
      { l:'D862', pos:asp.clone().add(new THREE.Vector3(1.2,0.6,-0.4)), mg:mgBpos },
    ];
    for (const res of catRes) {
      const dir = res.pos.clone().sub(asp).normalize();
      const end = res.pos.clone().add(dir.multiplyScalar(0.8));
      this.groups.metals.add(this._stick(res.pos,end,0x707070,0.05));
      const perp = new THREE.Vector3(-dir.y,dir.x,dir.z*0.5).normalize();
      for (const s of [1,-1]) {
        const o = end.clone().add(perp.clone().multiplyScalar(s*0.35));
        this.groups.metals.add(this._atom(0xc03030,0.12));
        this.groups.metals.children[this.groups.metals.children.length-1].position.copy(o);
        this.groups.metals.add(this._stick(end,o,0xc03030,0.03));
        this.groups.metals.add(this._dash(o,res.mg,0x888888));
      }
      // Asp-Mg distance label -pushed away from center
      const distVal = res.pos.distanceTo(res.mg) / SCALE;
      const distMid = res.pos.clone().add(res.mg).multiplyScalar(0.5);
      const outDir = res.pos.clone().sub(asp).normalize();
      const dLabelPos = distMid.clone().add(outDir.multiplyScalar(6)).add(new THREE.Vector3(0, 4, 0));
      const dLabel = this._spriteLabel(`${res.l}: ${distVal.toFixed(1)} \u00c5`, '#aaaaaa', '#888888', 200);
      dLabel.scale.set(6, 1.8, 1);
      dLabel.position.copy(dLabelPos);
      dLabel.userData.labelPriority = 0;
      this.groups.metals.add(dLabel);
      this.groups.metals.add(this._leaderLine(dLabelPos, distMid, new THREE.Color(0x888888)));
    }

    const dBase = asp.clone().add(new THREE.Vector3(1.8,0.8,-0.6));
    const pA=dBase.clone(), pB=dBase.clone().add(new THREE.Vector3(0.65,0.15,0)),
          pG=dBase.clone().add(new THREE.Vector3(1.3,0.3,0.05));
    for (const pp of [pA,pB,pG]) {
      this.groups.dntp.add(this._atom(0xd07020,0.18));
      this.groups.dntp.children[this.groups.dntp.children.length-1].position.copy(pp);
    }
    for (const [a,b] of [[pA,pB],[pB,pG]]) {
      const mid=a.clone().add(b).multiplyScalar(0.5);
      this.groups.dntp.add(this._atom(0xc03030,0.09));
      this.groups.dntp.children[this.groups.dntp.children.length-1].position.copy(mid);
      this.groups.dntp.add(this._stick(a,mid,0xd07020,0.03));
      this.groups.dntp.add(this._stick(mid,b,0xd07020,0.03));
    }
    const sugarC = dBase.clone().add(new THREE.Vector3(-0.7,0,0));
    this.groups.dntp.add(this._atom(0x707070,0.12));
    this.groups.dntp.children[this.groups.dntp.children.length-1].position.copy(sugarC);
    this.groups.dntp.add(this._stick(sugarC,pA,0xd07020,0.03));
    const baseC = sugarC.clone().add(new THREE.Vector3(0,0.6,0));
    this.groups.dntp.add(this._atom(0x3050c0,0.1));
    this.groups.dntp.children[this.groups.dntp.children.length-1].position.copy(baseC);
    this.groups.dntp.add(this._stick(sugarC,baseC,0x3050c0,0.025));

    this.groups.dntp.add(this._dash(pB,mgBpos,0x60a830));
    this.groups.dntp.add(this._dash(pG,mgBpos,0x60a830));

    // Template nucleotide opposite the dNTP (Watson-Crick base pairing)
    const templateBase = baseC.clone().add(new THREE.Vector3(0, 1.7, 0));
    const templateAtom = this._atom(0x3050c0, 0.1);
    templateAtom.position.copy(templateBase);
    this.groups.dntp.add(templateAtom);
    const templateRing = templateBase.clone().add(new THREE.Vector3(0.3, 0.2, 0));
    const templateRingAtom = this._atom(0x3050c0, 0.08);
    templateRingAtom.position.copy(templateRing);
    this.groups.dntp.add(templateRingAtom);
    this.groups.dntp.add(this._stick(templateBase, templateRing, 0x3050c0, 0.025));

    // 3 Watson-Crick H-bond dashed lines
    const hbondColor = 0x88aacc;
    const hbondOffsets = [
      [0, 0, 0], [0.18, 0.06, 0.12], [-0.18, 0.06, -0.12]
    ];
    for (const off of hbondOffsets) {
      const from = baseC.clone().add(new THREE.Vector3(off[0], 0.1 + off[1], off[2]));
      const to = templateBase.clone().add(new THREE.Vector3(off[0], -0.1 + off[1], off[2]));
      const geo = new THREE.BufferGeometry().setFromPoints([from, to]);
      const mat = new THREE.LineDashedMaterial({
        color: hbondColor, dashSize: 0.08, gapSize: 0.06,
        transparent: true, opacity: 0.5
      });
      const line = new THREE.Line(geo, mat);
      line.computeLineDistances();
      this.groups.dntp.add(line);
    }

    // W-C base pair label -pushed out from model center
    const wcAnchor = baseC.clone().add(new THREE.Vector3(0, 0.7, 0));
    const wcLabelPos = wcAnchor.clone().add(new THREE.Vector3(8, 6, -6));
    const wcLabel = this._spriteLabel('W\u2013C base pair', '#88aacc', '#88aacc', 220);
    wcLabel.scale.set(7, 2, 1);
    wcLabel.position.copy(wcLabelPos);
    wcLabel.userData.labelPriority = 1;
    this.groups.dntp.add(wcLabel);
    this.groups.dntp.add(this._leaderLine(wcLabelPos, wcAnchor, new THREE.Color(0x88aacc)));

    this.groups.dntp.userData = { name:'Incoming dNTP',
      detail:'Templated nucleotide \u00b7 \u03b2,\u03b3-phosphates \u2192 Metal B' };
    this.activeSiteLight.position.copy(asp);
  }

  _atom(color,r) {
    return new THREE.Mesh(new THREE.SphereGeometry(r,12,12),
      new THREE.MeshPhysicalMaterial({ color, roughness:0.25, metalness:0.5,
        emissive:color, emissiveIntensity:0.25 }));
  }
  _stick(a,b,color,r) {
    const d=b.clone().sub(a);
    const m=new THREE.Mesh(new THREE.CylinderGeometry(r,r,d.length(),6),
      new THREE.MeshPhysicalMaterial({ color, roughness:0.4, metalness:0.3 }));
    m.position.copy(a.clone().add(b).multiplyScalar(0.5));
    m.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0),d.normalize());
    return m;
  }
  _dash(a,b,color) {
    const g=new THREE.BufferGeometry().setFromPoints([a,b]);
    const m=new THREE.LineDashedMaterial({ color, dashSize:0.2, gapSize:0.15, transparent:true, opacity:0.4 });
    const l=new THREE.Line(g,m); l.computeLineDistances(); return l;
  }

  // ==================== ZINC-BINDING MOTIFS (CTD) ====================

  _buildZincMotifs() {
    this.groups.zincMotifs = new THREE.Group();
    this.scene.add(this.groups.zincMotifs);

    const ctdCenter = new THREE.Vector3(...DOMAINS.ctd.center);
    const sites = [
      { name: 'CysA Zn-finger', offset: [-2.5, 1, -1.5] },
      { name: 'CysB Zn-finger', offset: [2, -1.5, 1] },
    ];
    const znPositions = [];

    for (const site of sites) {
      const znPos = ctdCenter.clone().add(new THREE.Vector3(...site.offset));
      znPositions.push(znPos);

      // Zn²⁺ ion
      const znSphere = new THREE.Mesh(
        new THREE.SphereGeometry(0.3, 16, 16),
        new THREE.MeshPhysicalMaterial({
          color: 0x7799aa, roughness: 0.2, metalness: 0.7,
          emissive: 0x7799aa, emissiveIntensity: 0.4
        })
      );
      znSphere.position.copy(znPos);
      znSphere.userData = { name: `Zn\u00b2\u207a (${site.name})`, detail: 'Tetrahedral Cys\u2084 coordination' };
      this.groups.zincMotifs.add(znSphere);

      // 4 coordinating cysteine sulfurs in tetrahedral geometry
      const tetDirs = [
        [1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]
      ];
      const bondLen = 2.3 * SCALE;
      for (const td of tetDirs) {
        const sDir = new THREE.Vector3(...td).normalize();
        const sPos = znPos.clone().add(sDir.clone().multiplyScalar(bondLen));

        // Yellow sulfur atom
        const sulfur = this._atom(0xcccc00, 0.12);
        sulfur.position.copy(sPos);
        this.groups.zincMotifs.add(sulfur);

        // Zn-S stick
        this.groups.zincMotifs.add(this._stick(znPos, sPos, 0x99aa88, 0.04));

        // C-alpha stub
        const caPos = sPos.clone().add(sDir.clone().multiplyScalar(0.5));
        const ca = this._atom(0x707070, 0.08);
        ca.position.copy(caPos);
        this.groups.zincMotifs.add(ca);
        this.groups.zincMotifs.add(this._stick(sPos, caPos, 0x707070, 0.025));
      }

      // Label -pushed away from model
      const labelPos = znPos.clone().add(new THREE.Vector3(
        site.offset[0] > 0 ? 4 : -4, 6, site.offset[2] > 0 ? 3 : -3));
      const label = this._spriteLabel(site.name, '#d0d8e0', '#7799aa', 260);
      label.scale.set(9, 2.5, 1);
      label.position.copy(labelPos);
      label.userData.labelPriority = 2;
      this.groups.zincMotifs.add(label);
      this.groups.zincMotifs.add(this._leaderLine(labelPos, znPos, new THREE.Color(0x7799aa)));
    }

    // Bridging dash between the two Zn sites
    if (znPositions.length === 2) {
      this.groups.zincMotifs.add(this._dash(znPositions[0], znPositions[1], 0x7799aa));
    }
  }

  // ==================== [4Fe-4S] CLUSTER (P DOMAIN) ====================

  _buildFeSCluster() {
    this.groups.fesCluster = new THREE.Group();
    this.scene.add(this.groups.fesCluster);

    const pCenter = new THREE.Vector3(...DOMAINS.pdomain.center);
    const clusterPos = pCenter.clone().add(new THREE.Vector3(0.5, 0.8, -0.5));

    // Cubane [4Fe-4S] core: alternating Fe and S at cube vertices
    const cubeR = 1.8 * SCALE;
    const fePositions = [
      clusterPos.clone().add(new THREE.Vector3( cubeR,  cubeR,  cubeR)),
      clusterPos.clone().add(new THREE.Vector3(-cubeR, -cubeR,  cubeR)),
      clusterPos.clone().add(new THREE.Vector3(-cubeR,  cubeR, -cubeR)),
      clusterPos.clone().add(new THREE.Vector3( cubeR, -cubeR, -cubeR)),
    ];
    const sPositions = [
      clusterPos.clone().add(new THREE.Vector3(-cubeR,  cubeR,  cubeR)),
      clusterPos.clone().add(new THREE.Vector3( cubeR, -cubeR,  cubeR)),
      clusterPos.clone().add(new THREE.Vector3( cubeR,  cubeR, -cubeR)),
      clusterPos.clone().add(new THREE.Vector3(-cubeR, -cubeR, -cubeR)),
    ];

    // Fe atoms (rust-orange)
    for (let i = 0; i < 4; i++) {
      const fe = new THREE.Mesh(
        new THREE.SphereGeometry(0.2, 12, 12),
        new THREE.MeshPhysicalMaterial({
          color: 0xcc6633, roughness: 0.3, metalness: 0.8,
          emissive: 0xcc6633, emissiveIntensity: 0.4
        })
      );
      fe.position.copy(fePositions[i]);
      fe.userData = { name: `Fe${i+1} ([4Fe-4S])`, detail: 'Iron-sulfur cluster / CysX motif' };
      this.groups.fesCluster.add(fe);
    }

    // Inorganic S atoms (yellow)
    for (let i = 0; i < 4; i++) {
      const s = new THREE.Mesh(
        new THREE.SphereGeometry(0.18, 12, 12),
        new THREE.MeshPhysicalMaterial({
          color: 0xcccc33, roughness: 0.4, metalness: 0.5,
          emissive: 0xcccc33, emissiveIntensity: 0.3
        })
      );
      s.position.copy(sPositions[i]);
      this.groups.fesCluster.add(s);
    }

    // Fe-S bonds within cubane
    for (const feP of fePositions) {
      for (const sP of sPositions) {
        if (feP.distanceTo(sP) < cubeR * 2.2) {
          this.groups.fesCluster.add(this._stick(feP, sP, 0x888866, 0.03));
        }
      }
    }

    // 4 coordinating cysteine sulfurs (CysX motif)
    const cysDirs = [[1.2,0.5,0],[-0.3,1.2,0.8],[0.5,-0.3,1.2],[-1.0,-0.8,-0.3]];
    for (let i = 0; i < 4; i++) {
      const cysS = fePositions[i].clone().add(
        new THREE.Vector3(cysDirs[i][0], cysDirs[i][1], cysDirs[i][2]).normalize().multiplyScalar(2.3 * SCALE)
      );
      const atom = this._atom(0xcccc33, 0.12);
      atom.position.copy(cysS);
      this.groups.fesCluster.add(atom);
      this.groups.fesCluster.add(this._stick(fePositions[i], cysS, 0xcccc33, 0.025));
    }

    // Label
    const labelPos = clusterPos.clone().add(new THREE.Vector3(6, 14, -6));
    const label = this._spriteLabel('[4Fe-4S] cluster', '#d0d8e0', '#cc6633', 240);
    label.scale.set(10, 2.8, 1);
    label.position.copy(labelPos);
    label.userData.labelPriority = 2;
    this.groups.fesCluster.add(label);
    this.groups.fesCluster.add(this._leaderLine(labelPos, clusterPos, new THREE.Color(0xcc6633)));
  }

  // ==================== EXONUCLEASE ACTIVE SITE ====================

  _buildExoSite() {
    this.groups.exoSite = new THREE.Group();
    this.scene.add(this.groups.exoSite);

    const exoCenter = new THREE.Vector3(...DOMAINS.exo.center);
    const mgOffsets = [[0.3, -0.5, 1.0], [-0.8, -0.3, 1.5]];
    const mgPositions = mgOffsets.map(o => exoCenter.clone().add(new THREE.Vector3(...o)));
    this.exoSitePos = exoCenter.clone();

    // Two Mg²⁺ ions with octahedral coordination
    for (let m = 0; m < 2; m++) {
      const mgPos = mgPositions[m];
      const mgSphere = new THREE.Mesh(
        new THREE.SphereGeometry(0.35, 16, 16),
        new THREE.MeshPhysicalMaterial({
          color: 0x60a830, roughness: 0.2, metalness: 0.7,
          emissive: 0x60a830, emissiveIntensity: 0.5
        })
      );
      mgSphere.position.copy(mgPos);
      mgSphere.userData = { name: `Mg\u00b2\u207a (EXO Metal ${m === 0 ? 'A' : 'B'})`,
        detail: 'Exonuclease two-metal-ion center' };
      this.groups.exoSite.add(mgSphere);

      // 6 octahedral ligand spheres
      const octDirs = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]];
      const ligDist = 2.1 * SCALE;
      for (const d of octDirs) {
        const lPos = mgPos.clone().add(new THREE.Vector3(...d).multiplyScalar(ligDist));
        const lig = this._atom(0xc03030, 0.08);
        lig.position.copy(lPos);
        this.groups.exoSite.add(lig);
        this.groups.exoSite.add(this._dash(mgPos, lPos, 0x60a830));
      }
    }

    // Mg-Mg bridging dash
    this.groups.exoSite.add(this._dash(mgPositions[0], mgPositions[1], 0x60a830));

    // Catalytic residues: D275, E277, D368, D462 (estimated)
    const catRes = [
      { l: 'D275', offset: [-1.5, 0.8, 0.5] },
      { l: 'E277', offset: [-1.2, -0.6, 1.8] },
      { l: 'D368', offset: [0.8, 0.5, 2.0] },
      { l: 'D462', offset: [1.2, -0.8, 0.8] },
    ];
    for (const res of catRes) {
      const resPos = exoCenter.clone().add(new THREE.Vector3(...res.offset));
      // Ball-and-stick
      const ca = this._atom(0x707070, 0.12);
      ca.position.copy(resPos);
      this.groups.exoSite.add(ca);

      const dir = resPos.clone().sub(exoCenter).normalize();
      const end = resPos.clone().add(dir.clone().multiplyScalar(0.6));
      this.groups.exoSite.add(this._stick(resPos, end, 0x707070, 0.04));

      // Carboxylate oxygens
      const perp = new THREE.Vector3(-dir.y, dir.x, dir.z * 0.5).normalize();
      for (const s of [1, -1]) {
        const oPos = end.clone().add(perp.clone().multiplyScalar(s * 0.3));
        const oxy = this._atom(0xc03030, 0.1);
        oxy.position.copy(oPos);
        this.groups.exoSite.add(oxy);
        this.groups.exoSite.add(this._stick(end, oPos, 0xc03030, 0.025));
        // Coordination dash to nearest metal
        const nearMg = mgPositions[0].distanceTo(oPos) < mgPositions[1].distanceTo(oPos)
          ? mgPositions[0] : mgPositions[1];
        this.groups.exoSite.add(this._dash(oPos, nearMg, 0x888888));
      }
    }

    // Label -pushed away from model
    const exoLabelPos = exoCenter.clone().add(new THREE.Vector3(-8, 16, 8));
    const label = this._spriteLabel('EXO SITE \u00b7 2\u00d7 Mg\u00b2\u207a', '#d0d8e0', '#b85c5c', 300);
    label.scale.set(11, 3, 1);
    label.position.copy(exoLabelPos);
    label.userData.labelPriority = 2;
    this.groups.exoSite.add(label);
    this.groups.exoSite.add(this._leaderLine(exoLabelPos, exoCenter, PAL.exo));
  }

  // ==================== ACCESSORY SUBUNITS ====================

  _buildAccessorySubunits() {
    this.groups.accessory = new THREE.Group();
    this.groups.accessory.visible = true;
    this.scene.add(this.groups.accessory);

    const ctdCenter = new THREE.Vector3(...DOMAINS.ctd.center);
    const subunits = [
      { name: 'POLE2 (p59)', offset: [8, 2, -3], radius: 5.0, color: 0x6699aa },
      { name: 'POLE3 (p17)', offset: [10, -3, 2], radius: 2.8, color: 0x99aa66 },
      { name: 'POLE4 (p12)', offset: [10, -3, -3], radius: 2.2, color: 0xaa8866 },
    ];

    for (const sub of subunits) {
      const pos = ctdCenter.clone().add(new THREE.Vector3(...sub.offset));

      // Ghost surface
      const ghost = new THREE.Mesh(
        new THREE.IcosahedronGeometry(sub.radius, 3),
        new THREE.MeshPhysicalMaterial({
          color: sub.color, transparent: true, opacity: 0.06,
          roughness: 0.9, side: THREE.DoubleSide, depthWrite: false
        })
      );
      ghost.position.copy(pos);
      ghost.userData = { name: sub.name, detail: 'Accessory subunit' };
      this.groups.accessory.add(ghost);

      // Wireframe outline
      const wire = new THREE.Mesh(
        new THREE.IcosahedronGeometry(sub.radius, 2),
        new THREE.MeshBasicMaterial({
          color: sub.color, wireframe: true, transparent: true, opacity: 0.08
        })
      );
      wire.position.copy(pos);
      this.groups.accessory.add(wire);

      // Label -pushed away from model
      const subLabelPos = pos.clone().add(new THREE.Vector3(
        sub.offset[0] > 9 ? 4 : -3, sub.radius + 5, sub.offset[2] > 0 ? 3 : -3));
      const label = this._spriteLabel(sub.name, '#d0d8e0', '#' + new THREE.Color(sub.color).getHexString(), 240);
      label.scale.set(9, 2.5, 1);
      label.position.copy(subLabelPos);
      label.userData.labelPriority = 2;
      this.groups.accessory.add(label);
      this.groups.accessory.add(this._leaderLine(subLabelPos, pos, new THREE.Color(sub.color)));

      // Dashed line from CTD center to subunit
      this.groups.accessory.add(this._dash(ctdCenter, pos, sub.color));
    }
  }

  // ==================== PCNA SLIDING CLAMP ====================

  _buildPCNA() {
    this.groups.pcna = new THREE.Group();
    this.groups.pcna.visible = true;
    this.scene.add(this.groups.pcna);

    const thumbCenter = new THREE.Vector3(...DOMAINS.thumb.center);
    const pcnaPos = thumbCenter.clone().add(new THREE.Vector3(5, -6, 10));

    // Torus geometry for PCNA ring
    const torusGeo = new THREE.TorusGeometry(4.0, 1.2, 16, 48);
    const ghostMat = new THREE.MeshPhysicalMaterial({
      color: 0x7abba8, transparent: true, opacity: 0.12,
      emissive: 0x7abba8, emissiveIntensity: 0.15,
      roughness: 0.8, side: THREE.DoubleSide, depthWrite: false
    });
    const ghostTorus = new THREE.Mesh(torusGeo, ghostMat);
    ghostTorus.position.copy(pcnaPos);
    ghostTorus.rotation.y = Math.PI / 2;
    ghostTorus.rotation.x = 0.2;
    ghostTorus.userData = { name: 'PCNA sliding clamp', detail: 'Homotrimeric processivity factor · Pol \u03b5 contacts all three protomers via PIP-box (Q1180), Thumb insertion (res. 1102\u20131122), and P domain \u2014 tripartite interface impervious to FEN1/Pol \u03b4 displacement (Roske & Yeeles 2024)' };
    this.groups.pcna.add(ghostTorus);

    // Wireframe torus
    const wireTorus = new THREE.Mesh(
      new THREE.TorusGeometry(4.0, 1.2, 8, 24),
      new THREE.MeshBasicMaterial({
        color: 0x7abba8, wireframe: true, transparent: true, opacity: 0.1
      })
    );
    wireTorus.position.copy(pcnaPos);
    wireTorus.rotation.y = Math.PI / 2;
    wireTorus.rotation.x = 0.2;
    this.groups.pcna.add(wireTorus);

    // 3 symmetry markers (PCNA homotrimer)
    for (let i = 0; i < 3; i++) {
      const ang = (i / 3) * Math.PI * 2;
      const markerPos = pcnaPos.clone().add(new THREE.Vector3(
        Math.cos(ang) * 4.0 * Math.cos(0.2),
        Math.sin(0.2) * Math.cos(ang) * 4.0,
        Math.sin(ang) * 4.0
      ));
      const marker = new THREE.Mesh(
        new THREE.SphereGeometry(0.25, 8, 8),
        new THREE.MeshPhysicalMaterial({
          color: 0x7abba8, emissive: 0x7abba8, emissiveIntensity: 0.5,
          roughness: 0.3
        })
      );
      marker.position.copy(markerPos);
      this.groups.pcna.add(marker);
    }

    // Tripartite interface labels
    const triLabels = [
      { text:'PIP-box (Q1180)', offset:[-7, -5, 3] },
      { text:'Thumb insert (1102\u20131122)', offset:[5, -4, -5] },
      { text:'P-domain contact', offset:[-5, 4, -6] },
    ];
    for (const tl of triLabels) {
      const tlPos = pcnaPos.clone().add(new THREE.Vector3(...tl.offset));
      const tlLabel = this._spriteLabel(tl.text, '#7abba8', '#7abba8', 240);
      tlLabel.scale.set(9, 2.5, 1);
      tlLabel.position.copy(tlPos);
      tlLabel.userData.labelPriority = 1;
      this.groups.pcna.add(tlLabel);
      this.groups.pcna.add(this._leaderLine(tlPos, pcnaPos, new THREE.Color(0x7abba8)));
    }

    // Main label
    const mainLabelPos = pcnaPos.clone().add(new THREE.Vector3(0, 10, 0));
    const mainLabel = this._spriteLabel('PCNA sliding clamp', '#d0d8e0', '#7abba8', 280);
    mainLabel.scale.set(11, 3, 1);
    mainLabel.position.copy(mainLabelPos);
    mainLabel.userData.labelPriority = 2;
    this.groups.pcna.add(mainLabel);
    this.groups.pcna.add(this._leaderLine(mainLabelPos, pcnaPos, new THREE.Color(0x7abba8)));

    // Dashed line from thumb center to PCNA
    this.groups.pcna.add(this._dash(thumbCenter, pcnaPos, 0x7abba8));
  }

  // ==================== MUTATION MARKERS ====================

  _buildMutationMarkers() {
    this.groups.mutations = new THREE.Group();
    this.scene.add(this.groups.mutations);

    for (const mut of MUTATIONS) {
      const domCenter = new THREE.Vector3(...DOMAINS[mut.domain].center);
      const pos = domCenter.clone().add(new THREE.Vector3(...mut.offset));
      const g = new THREE.Group();

      const ring = new THREE.Mesh(new THREE.RingGeometry(0.55,0.72,24),
        new THREE.MeshBasicMaterial({ color:PAL.mutMarker, side:THREE.DoubleSide,
          transparent:true, opacity:0.8 }));
      g.add(ring);
      g.add(new THREE.Mesh(new THREE.SphereGeometry(0.25,12,12),
        new THREE.MeshPhysicalMaterial({ color:PAL.mutMarker,
          emissive:PAL.mutMarker, emissiveIntensity:0.8, roughness:0.25, metalness:0.4 })));
      const halo = new THREE.Mesh(new THREE.SphereGeometry(1.0,12,12),
        new THREE.MeshBasicMaterial({ color:PAL.mutMarker, transparent:true, opacity:0.1, depthWrite:false }));
      g.add(halo);
      g.position.copy(pos);
      // Enrich tooltip detail with pipeline data (signatures, classification, DDG)
      let enrichedDetail = mut.detail;
      if (mut.classification) {
        enrichedDetail += ` / ${mut.classification.acmg_class} (P=${mut.classification.posterior})`;
      }
      if (mut.ddg !== undefined) {
        enrichedDetail += ` / \u0394\u0394G=${mut.ddg} kcal/mol`;
      }
      if (mut.signature) {
        const topSig = Object.entries(mut.signature).sort((a,b) => b[1]-a[1])[0];
        if (topSig) enrichedDetail += ` / ${topSig[0]}:${(topSig[1]*100).toFixed(0)}%`;
      }
      g.userData = { halo, ring, name:mut.label, detail:enrichedDetail, mut:mut.id, pulseBoost:0 };
      this.groups.mutations.add(g);
      this.mutMarkers[mut.id] = g;
    }
  }

  // ==================== ANNOTATIONS (domain headers, labels) ====================

  _buildAnnotations() {
    this.groups.labels = new THREE.Group();
    this.scene.add(this.groups.labels);

    // Stagger labels radially outward and high above the model
    const labelOffsets = {
      ntd:     { dx: -8, dy: 16, dz: -6 },
      exo:     { dx: -8, dy: 16, dz:  8 },
      palm:    { dx:-12, dy: 22, dz:-10 },
      pdomain: { dx:  6, dy: 18, dz: -8 },
      fingers: { dx: 10, dy: 20, dz:  8 },
      thumb:   { dx:-10, dy:-16, dz: 14 },
      inactpol:{ dx: 14, dy: 14, dz:  4 },
      ctd:     { dx:  8, dy: 16, dz: -6 },
    };

    for (const key of DOMAIN_ORDER) {
      const dom = DOMAINS[key];
      const col = PAL[dom.color];
      const hex = '#'+col.getHexString();
      const off = labelOffsets[key];
      const anchor = new THREE.Vector3(...dom.center);
      const labelPos = anchor.clone().add(new THREE.Vector3(off.dx, off.dy, off.dz));

      const label = this._spriteLabel(dom.abbrev, '#e8edf2', hex, 0, dom.range[0] + '\u2013' + dom.range[1]);
      label.position.copy(labelPos);
      label.userData.labelPriority = 3;
      this.groups.labels.add(label);
      this.groups.labels.add(this._leaderLine(labelPos, anchor, col));
    }

    const asAnchor = this.activeSitePos.clone();
    const asPos = asAnchor.clone().add(new THREE.Vector3(-10, 20, -10));
    const asLabel = this._spriteLabel('ACTIVE SITE \u00b7 2\u00d7 Mg\u00b2\u207a', '#d0d8e0',
      '#'+PAL.magnesium.getHexString());
    asLabel.scale.set(12,3,1);
    asLabel.position.copy(asPos);
    asLabel.userData.labelPriority = 3;
    this.groups.labels.add(asLabel);
    this.groups.labels.add(this._leaderLine(asPos, asAnchor, PAL.magnesium));

    const proofAnchor = new THREE.Vector3(-10, 2, 7);
    const proofPos = proofAnchor.clone().add(new THREE.Vector3(-8, 18, 8));
    const proofLabel = this._spriteLabel('PROOFREADING \u00b7 ExoI/II/III', '#d0d8e0', '#b85c5c');
    proofLabel.scale.set(12,3,1);
    proofLabel.position.copy(proofPos);
    proofLabel.userData.labelPriority = 2;
    this.groups.labels.add(proofLabel);
    this.groups.labels.add(this._leaderLine(proofPos, proofAnchor, PAL.exo));

    const annotations = [
      { text:'D626 (Motif A)', pos:[0,-1,0], offset:[-12,-14,10], color:'#6b8e6b' },
      { text:'D860/D862 (Motif C)', pos:[0,-1,0], offset:[12,-14,-8], color:'#6b8e6b' },
      { text:'O-helix', pos:[7,5,4], offset:[8,16,6], color:'#5b7f99' },
    ];
    for (const a of annotations) {
      const anchor = new THREE.Vector3(...a.pos);
      const lPos = anchor.clone().add(new THREE.Vector3(...a.offset));
      const lbl = this._spriteLabel(a.text, a.color, a.color, 280);
      lbl.scale.set(9,2.5,1);
      lbl.position.copy(lPos);
      lbl.userData.labelPriority = 1;
      this.groups.labels.add(lbl);
      this.groups.labels.add(this._leaderLine(lPos, anchor, new THREE.Color(a.color)));
    }
  }

  _leaderLine(from, to, color) {
    const geo = new THREE.BufferGeometry().setFromPoints([from, to]);
    const mat = new THREE.LineDashedMaterial({
      color, dashSize: 0.3, gapSize: 0.2,
      transparent: true, opacity: 0.45
    });
    const line = new THREE.Line(geo, mat);
    line.computeLineDistances();
    return line;
  }

  _spriteLabel(text, textColor, accentColor, width, sub) {
    // Fixed screen-size label: text-fitted pill, constant pixel height regardless of zoom.
    const dpr = 3, fs = 12, padX = 9, dot = 6, gap = 7, h = 22;
    const font = '500 ' + fs + 'px "IBM Plex Mono", ui-monospace, monospace';
    const subFont = '400 ' + fs + 'px "IBM Plex Mono", ui-monospace, monospace';
    const m = document.createElement('canvas').getContext('2d');
    m.font = font; const tw = m.measureText(text).width;
    let sw = 0; if (sub) { m.font = subFont; sw = m.measureText(sub).width + 9; }
    const w = Math.ceil(padX + dot + gap + tw + sw + padX);
    const c = document.createElement('canvas');
    c.width = w * dpr; c.height = h * dpr;
    const ctx = c.getContext('2d');
    ctx.scale(dpr, dpr);
    const r = 5;
    ctx.beginPath(); ctx.roundRect(0.5, 0.5, w - 1, h - 1, r);
    ctx.fillStyle = 'rgba(14,18,24,0.9)'; ctx.fill();
    ctx.globalAlpha = 0.6; ctx.strokeStyle = accentColor; ctx.lineWidth = 1; ctx.stroke();
    ctx.globalAlpha = 1;
    ctx.beginPath(); ctx.arc(padX + dot / 2, h / 2, dot / 2, 0, Math.PI * 2);
    ctx.fillStyle = accentColor; ctx.fill();
    ctx.textBaseline = 'middle';
    ctx.font = font; ctx.fillStyle = '#eef2f5';
    ctx.fillText(text, padX + dot + gap, h / 2 + 0.5);
    if (sub) { ctx.font = subFont; ctx.fillStyle = '#93a0ad'; ctx.fillText(sub, padX + dot + gap + tw + 9, h / 2 + 0.5); }
    const tex = new THREE.CanvasTexture(c);
    tex.colorSpace = THREE.SRGBColorSpace;
    tex.anisotropy = 4;
    const sprite = new THREE.Sprite(new THREE.SpriteMaterial({
      map: tex, transparent: true, depthTest: false, depthWrite: false, sizeAttenuation: false, fog: false }));
    sprite.renderOrder = 20;
    sprite.userData.labelPriority = 0;
    sprite.userData.baseOpacity = 1;
    sprite.userData.pxW = w; sprite.userData.pxH = h;
    sprite.material.opacity = 0;
    this._labelSprites.push(sprite);
    return sprite;
  }

  setLabelDensity(level) {
    this.labelMinPriority = ({ essential: 3, standard: 2, all: 0 })[level] ?? 2;
  }

  // ==================== LABEL OVERLAP CULLING ====================

  _updateLabelVisibility() {
    const el = this.renderer.domElement;
    const W = el.clientWidth, H = el.clientHeight;
    if (!W || !H) return;
    const P11 = this.camera.projectionMatrix.elements[5];
    const labelsOn = !this.groups.labels || this.groups.labels.visible;
    const minP = this.labelMinPriority ?? 0;
    const v = new THREE.Vector3();
    const items = [];
    for (const s of this._labelSprites) {
      if (!s.parent) continue;
      const ls = this.labelScale || 1;
      const sy = 2 * s.userData.pxH * ls / (P11 * H);
      s.scale.set(sy * s.userData.pxW / s.userData.pxH, sy, 1);
      let vis = true;
      for (let p = s.parent; p; p = p.parent) if (!p.visible) { vis = false; break; }
      const meas = !!s.userData.isMeasure;
      const pr = s.userData.labelPriority || 0;
      if (!vis || (!meas && (!labelsOn || pr < minP))) { s.visible = false; s.material.opacity = 0; continue; }
      v.setFromMatrixPosition(s.matrixWorld).project(this.camera);
      if (v.z > 1 || v.z < -1) { s.visible = false; continue; }
      items.push({ s, x: (v.x + 1) * W / 2, y: (-v.y + 1) * H / 2, z: v.z,
        pr: meas ? 99 : pr, hw: s.userData.pxW * ls / 2 + 4, hh: s.userData.pxH * ls / 2 + 3 });
    }
    items.sort((a, b) => b.pr - a.pr || a.z - b.z);
    const placed = [];
    for (const p of items) {
      let ov = false;
      for (const q of placed) {
        if (Math.abs(p.x - q.x) < p.hw + q.hw && Math.abs(p.y - q.y) < p.hh + q.hh) { ov = true; break; }
      }
      const o = p.s.material.opacity + ((ov ? 0 : 1) - p.s.material.opacity) * 0.18;
      p.s.material.opacity = o;
      p.s.visible = o > 0.02;
      if (!ov) placed.push(p);
    }
  }

  // ==================== HYDROGEN BONDS (i→i+4 in helices) ====================

  _buildHydrogenBonds() {
    this.groups.hbonds = new THREE.Group();
    this.groups.hbonds.visible = true;
    this.scene.add(this.groups.hbonds);

    // Batch all helix H-bond lines per domain into a single LineSegments call
    for (const key of DOMAIN_ORDER) {
      const bb = this.backboneCache[key];
      if (!bb) continue;
      const domColor = PAL[DOMAINS[key].color];
      const points = [];

      for (let i = 0; i < bb.length - 4; i++) {
        if (bb[i].ss === 'H' && bb[i+4].ss === 'H') {
          points.push(bb[i].pos, bb[i+4].pos);
        }
      }
      if (points.length === 0) continue;

      const geo = new THREE.BufferGeometry().setFromPoints(points);
      const mat = new THREE.LineBasicMaterial({
        color: domColor, transparent: true, opacity: 0.35
      });
      const lines = new THREE.LineSegments(geo, mat);
      this.groups.hbonds.add(lines);
    }
  }

  // ==================== β-SHEET HYDROGEN BONDS ====================

  _buildBetaSheetHBonds() {
    this.groups.sheetHbonds = new THREE.Group();
    this.groups.sheetHbonds.visible = true;
    this.scene.add(this.groups.sheetHbonds);

    // Batch all beta-sheet H-bond lines per domain into a single LineSegments call
    for (const key of DOMAIN_ORDER) {
      const bb = this.backboneCache[key];
      if (!bb) continue;
      const segments = segmentBySSType(bb);
      const strands = segments.filter(s => s.type === 'E');
      const sheetPts = [];

      // Find spatially proximal strand pairs
      for (let i = 0; i < strands.length; i++) {
        for (let j = i + 1; j < strands.length; j++) {
          const sA = strands[i], sB = strands[j];
          const midA = new THREE.Vector3();
          for (let k = sA.startIdx; k <= sA.endIdx; k++) midA.add(bb[k].pos);
          midA.divideScalar(sA.endIdx - sA.startIdx + 1);
          const midB = new THREE.Vector3();
          for (let k = sB.startIdx; k <= sB.endIdx; k++) midB.add(bb[k].pos);
          midB.divideScalar(sB.endIdx - sB.startIdx + 1);

          if (midA.distanceTo(midB) > 8.0) continue;

          const lenA = sA.endIdx - sA.startIdx + 1;
          const lenB = sB.endIdx - sB.startIdx + 1;
          const nBonds = Math.min(lenA, lenB);
          for (let k = 0; k < nBonds; k += 2) {
            const idxA = sA.startIdx + k;
            const idxB = sB.startIdx + Math.min(k, lenB - 1);
            sheetPts.push(bb[idxA].pos, bb[idxB].pos);
          }
        }
      }

      if (sheetPts.length === 0) continue;
      const geo = new THREE.BufferGeometry().setFromPoints(sheetPts);
      const mat = new THREE.LineBasicMaterial({
        color: PAL[DOMAINS[key].color], transparent: true, opacity: 0.3
      });
      this.groups.sheetHbonds.add(new THREE.LineSegments(geo, mat));
    }
  }

  // ==================== EXO↔POL SWITCHING PATH ====================

  _buildExoPolPath() {
    this.groups.exoPolPath = new THREE.Group();
    this.groups.exoPolPath.visible = true;
    this.scene.add(this.groups.exoPolPath);
    // Schematic waypoints for the four mismatch-proofreading intermediates
    // (Roske & Yeeles 2024: 9F6I post-insertion, 9F6J arrest, 9F6K frayed substrate, 9F6L mismatch excision).
    // Handover runs Thumb -> P domain (frayed) -> back to Thumb (excision); positions are illustrative only.
    const pol = this.activeSitePos.clone();
    const exo = new THREE.Vector3(...DOMAINS.exo.center);
    const pdom = new THREE.Vector3(...DOMAINS.pdomain.center);
    const thumb = new THREE.Vector3(...DOMAINS.thumb.center);
    const pts = [pol.clone(), pol.clone().lerp(thumb, 0.18), pdom.clone().add(new THREE.Vector3(0, 2.5, 0)), exo.clone()];
    this.proofPts = pts;
    this.proofCurve = new THREE.CatmullRomCurve3(pts, false, 'catmullrom', 0.5);
    this.exoPolCurve = null; // legacy straight tween disabled
    this._proofU = 0;
    this.groups.exoPolPath.add(new THREE.Mesh(new THREE.TubeGeometry(this.proofCurve, 90, 0.07, 6, false),
      new THREE.MeshPhysicalMaterial({ color:0x7abba8, transparent:true, opacity:0.35, roughness:0.5 })));
    const names = ['1 \u00b7 Post-insertion', '2 \u00b7 Arrest', '3 \u00b7 Frayed substrate', '4 \u00b7 Mismatch excision'];
    const pdb = ['9F6I', '9F6J', '9F6K', '9F6L'];
    const offs = [[-5, -6, 4], [6, -7, 3], [6, 6, -4], [-6, 7, 5]];
    this.proofNodes = [];
    pts.forEach((p, i) => {
      const node = new THREE.Mesh(new THREE.SphereGeometry(0.22, 12, 12),
        new THREE.MeshBasicMaterial({ color:0x7abba8, transparent:true, opacity:0.7 }));
      node.position.copy(p);
      this.groups.exoPolPath.add(node);
      this.proofNodes.push(node);
      const lp = p.clone().add(new THREE.Vector3(...offs[i]));
      const lbl = this._spriteLabel(names[i], '#e8edf2', '#7abba8', 0, pdb[i]);
      lbl.position.copy(lp);
      lbl.userData.labelPriority = 2;
      this.groups.exoPolPath.add(lbl);
      this.groups.exoPolPath.add(this._leaderLine(lp, p, new THREE.Color(0x7abba8)));
    });
    this.exoPolSphere = new THREE.Mesh(new THREE.SphereGeometry(0.42, 16, 16),
      new THREE.MeshPhysicalMaterial({ color:0xa6dccb, emissive:0x7abba8, emissiveIntensity:0.9, roughness:0.2, metalness:0.3 }));
    this.exoPolSphere.position.copy(pts[0]);
    this.groups.exoPolPath.add(this.exoPolSphere);
    const mid = pol.clone().add(exo).multiplyScalar(0.5).add(new THREE.Vector3(0, -6, 0));
    const distLabel = this._spriteLabel('Pol\u2194Exo ~40 \u00c5', '#e8edf2', '#7abba8', 0, 'cryo-EM, not model');
    distLabel.position.copy(mid);
    distLabel.userData.labelPriority = 1;
    this.groups.exoPolPath.add(distLabel);
  }

  setProofStep(i) {
    if (!this.proofCurve) return;
    const to = Math.max(0, Math.min(3, i)) / 3, from = this._proofU;
    const t0 = performance.now();
    this.proofNodes.forEach((n, k) => { n.scale.setScalar(k === i ? 1.8 : 1); n.material.opacity = k === i ? 1 : 0.55; });
    const step = () => {
      const p = Math.min(1, (performance.now() - t0) / 1100);
      const e = p < 0.5 ? 2 * p * p : 1 - Math.pow(-2 * p + 2, 2) / 2;
      this._proofU = from + (to - from) * e;
      this.exoPolSphere.position.copy(this.proofCurve.getPoint(this._proofU));
      if (p < 1) requestAnimationFrame(step);
    };
    step();
  }

  focusResidue(res) {
    const order = DOMAIN_ORDER.filter(k => res >= DOMAINS[k].range[0] && res <= DOMAINS[k].range[1])
      .sort((a, b) => (DOMAINS[a].range[1] - DOMAINS[a].range[0]) - (DOMAINS[b].range[1] - DOMAINS[b].range[0]));
    const key = order[0];
    if (!key) return null;
    const bb = this.backboneCache[key], [r0, r1] = DOMAINS[key].range;
    const idx = Math.round((res - r0) / Math.max(1, r1 - r0) * (bb.length - 1));
    const pos = bb[Math.max(0, Math.min(bb.length - 1, idx))].pos.clone();
    if (!this._resMarker) {
      const g = new THREE.Group();
      const ring = new THREE.Mesh(new THREE.RingGeometry(0.7, 0.9, 32),
        new THREE.MeshBasicMaterial({ color:0xf1d595, side:THREE.DoubleSide, transparent:true, opacity:0.95, depthTest:false }));
      ring.renderOrder = 15;
      g.add(ring);
      g.add(new THREE.Mesh(new THREE.SphereGeometry(0.28, 12, 12), new THREE.MeshBasicMaterial({ color:0xf1d595, depthTest:false })));
      g.userData.ring = ring;
      this.scene.add(g);
      this._resMarker = g;
    }
    this._resMarker.visible = true;
    this._resMarker.position.copy(pos);
    animateCamera(this, { r: 30, tx: pos.x, ty: pos.y, tz: pos.z });
    return { key, domains: order };
  }
  clearResidueMarker() { if (this._resMarker) this._resMarker.visible = false; }

  getView() { const t = this.cameraTarget; return [this.cameraTheta, this.cameraPhi, this.cameraRadius, t.x, t.y, t.z]; }
  setView(a) {
    if (!a || a.length < 6 || a.some(x => !isFinite(x))) return;
    [this.cameraTheta, this.cameraPhi, this.cameraRadius] = a;
    this.cameraTarget.set(a[3], a[4], a[5]);
    this._updateCamera();
  }
  setFigureMode(on) { if (this.scene.fog) this.scene.fog.color.set(on ? 0xf3f4f1 : PAL.bg); }
  pxPerAngstrom() {
    const H = this.renderer.domElement.clientHeight || 1;
    const worldPerPx = 2 * Math.tan(THREE.MathUtils.degToRad(this.camera.fov / 2)) * this.cameraRadius / H;
    return SCALE / worldPerPx;
  }

  _handleMutationClick(e) {
    if (!this.groups.mutations || !this.groups.mutations.visible) return false;
    const rect = this.renderer.domElement.getBoundingClientRect();
    const ndc = new THREE.Vector2(((e.clientX - rect.left) / rect.width) * 2 - 1, -((e.clientY - rect.top) / rect.height) * 2 + 1);
    this.raycaster.setFromCamera(ndc, this.camera);
    const hits = this.raycaster.intersectObjects(this.groups.mutations.children, true);
    if (!hits.length) return false;
    let o = hits[0].object;
    while (o && !(o.userData && o.userData.mut)) o = o.parent;
    if (!o) return false;
    o.userData.pulseBoost = 1;
    if (this.onMutationSelect) this.onMutationSelect(o.userData.mut);
    return true;
  }

  // ==================== PARTICLES ====================

  _buildParticles() {
    const n = 150;
    const pos = new Float32Array(n*3);
    for (let i=0;i<n;i++) { pos[i*3]=(Math.random()-0.5)*100; pos[i*3+1]=(Math.random()-0.5)*70; pos[i*3+2]=(Math.random()-0.5)*100; }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(pos,3));
    this.scene.add(new THREE.Points(geo,
      new THREE.PointsMaterial({ color:PAL.accent, size:0.04,
        transparent:true, opacity:0.15, depthWrite:false })));
  }

  // ==================== CONTROLS ====================

  _setupControls() {
    const el = this.renderer.domElement;
    this._pointerDownPos = { x: 0, y: 0 };

    attachCameraControls(this, el, {
      minRadius: 25,
      maxRadius: 200,
      onPointerDown: e => {
        if (this.measureMode) { this._handleMeasureClick(e); return false; }
        this._pointerDownPos = { x: e.clientX, y: e.clientY };
      },
      onLongPress: pos => this._hover(pos),
      onTouchEnd: () => { this.tooltipEl.style.opacity = '0'; },
    });

    el.addEventListener('click', e => {
      if (this.measureMode) return;
      const dx = e.clientX - this._pointerDownPos.x;
      const dy = e.clientY - this._pointerDownPos.y;
      if (dx * dx + dy * dy > 25) return;
      if (this._handleMutationClick(e)) return;
      this._handleDomainClick(e);
    });

    this.raycaster = new THREE.Raycaster();
    this.mouseNDC = new THREE.Vector2();
    el.addEventListener('mousemove', e => this._hover(e));
    el.style.cursor = 'grab';
  }

  _updateCamera() {
    const {cameraTheta:th,cameraPhi:ph,cameraRadius:r,cameraTarget:tgt} = this;
    this.camera.position.set(
      tgt.x + Math.cos(ph)*Math.sin(th)*r,
      tgt.y + Math.sin(ph)*r,
      tgt.z + Math.cos(ph)*Math.cos(th)*r);
    this.camera.lookAt(tgt);
  }

  _hover(e) {
    const rect = this.renderer.domElement.getBoundingClientRect();
    this.mouseNDC.set(((e.clientX-rect.left)/rect.width)*2-1,
      -((e.clientY-rect.top)/rect.height)*2+1);
    this.raycaster.setFromCamera(this.mouseNDC,this.camera);
    const targets = [];
    this.groups.protein.traverse(o => { if (o.isMesh&&o.userData.domain) targets.push(o); });
    this.groups.metals.traverse(o => { if (o.isMesh&&o.userData.name) targets.push(o); });
    this.groups.mutations.traverse(o => { if (o.isGroup&&o.userData.name) targets.push(o); });
    for (const gk of ['exoSite','zincMotifs','fesCluster','accessory','pcna']) {
      if (this.groups[gk] && this.groups[gk].visible) {
        this.groups[gk].traverse(o => { if (o.isMesh && o.userData.name) targets.push(o); });
      }
    }
    const hits = this.raycaster.intersectObjects(targets,true);
    if (hits.length) {
      let t=hits[0].object; while (t&&!t.userData.name) t=t.parent;
      if (t?.userData.name) {
        this.hoverTarget=t;
        let name = t.userData.name;
        let detail = t.userData.detail||'';

        // Per-residue enrichment for ribbon meshes
        const domainKey = t.userData.domainKey || t.userData.domain;
        if (domainKey && DOMAINS[domainKey] && hits[0].point) {
          const dom = DOMAINS[domainKey];
          const bb = this.backboneCache[domainKey];
          if (bb && bb.length > 0) {
            // Find closest backbone point to hit position
            const hitPt = hits[0].point;
            let closestIdx = 0;
            let closestDist = Infinity;
            for (let i = 0; i < bb.length; i++) {
              const d = hitPt.distanceTo(bb[i].pos);
              if (d < closestDist) { closestDist = d; closestIdx = i; }
            }
            const approxResidue = dom.range[0] + Math.round((closestIdx / Math.max(1, bb.length - 1)) * (dom.range[1] - dom.range[0]));
            const ssLabel = bb[closestIdx].ss === 'H' ? 'helix' : bb[closestIdx].ss === 'E' ? 'strand' : 'loop';
            name = `${dom.abbrev} · res ${approxResidue} · ${ssLabel}`;

            detail = dom.name;
          }
        }

        const hm = t.userData.mut || null;
        if (hm !== this._hoverMut) { this._hoverMut = hm; if (this.onMutationHover) this.onMutationHover(hm); }
        this.tooltipEl.querySelector('.t-name').textContent=name;
        this.tooltipEl.querySelector('.t-detail').textContent=detail;
        this.tooltipEl.style.opacity='1';
        this.tooltipEl.style.left=(e.clientX-rect.left+12)+'px';
        this.tooltipEl.style.top=(e.clientY-rect.top+12)+'px';
        this.renderer.domElement.style.cursor='pointer';
        return;
      }
    }
    this.hoverTarget=null;
    if (this._hoverMut) { this._hoverMut = null; if (this.onMutationHover) this.onMutationHover(null); }
    this.tooltipEl.style.opacity='0';
    this.renderer.domElement.style.cursor = this.measureMode ? 'crosshair' : (this.isDragging?'grabbing':'grab');
  }

  // ==================== CLICK-TO-SELECT DOMAIN ====================

  _handleDomainClick(e) {
    const rect = this.renderer.domElement.getBoundingClientRect();
    const ndc = new THREE.Vector2(
      ((e.clientX - rect.left) / rect.width) * 2 - 1,
      -((e.clientY - rect.top) / rect.height) * 2 + 1);
    this.raycaster.setFromCamera(ndc, this.camera);

    // Raycast against all ribbon meshes (they have domainKey in userData)
    const targets = [];
    this.groups.protein.traverse(o => {
      if (o.isMesh && (o.userData.domain || o.userData.domainKey)) targets.push(o);
    });
    const hits = this.raycaster.intersectObjects(targets);

    if (hits.length) {
      let obj = hits[0].object;
      const domainKey = obj.userData.domainKey || obj.userData.domain;
      if (domainKey) {
        const alreadySelected = (this._highlightedDomain === domainKey);
        if (alreadySelected) {
          // Deselect and zoom back out
          this.highlightDomain(domainKey, false);
          this.resetView();
        } else {
          // Select this domain and zoom in
          this.highlightDomain(domainKey, true);
          this.focusDomain(domainKey);
        }
        // Notify external listener (sidebar sync)
        if (this.onDomainSelect) {
          this.onDomainSelect(alreadySelected ? null : domainKey);
        }
        return;
      }
    }

    // Clicked empty space -deselect and zoom back out
    if (this._highlightedDomain) {
      this.highlightDomain(this._highlightedDomain, false);
      this.resetView();
      if (this.onDomainSelect) this.onDomainSelect(null);
    }
  }

  // ==================== DISTANCE MEASUREMENT ====================

  _handleMeasureClick(e) {
    const rect = this.renderer.domElement.getBoundingClientRect();
    const ndc = new THREE.Vector2(
      ((e.clientX-rect.left)/rect.width)*2-1,
      -((e.clientY-rect.top)/rect.height)*2+1);
    this.raycaster.setFromCamera(ndc, this.camera);

    const targets = [];
    this.groups.protein.traverse(o => { if (o.isMesh && o.visible) targets.push(o); });
    this.groups.dna.traverse(o => { if (o.isMesh) targets.push(o); });
    this.groups.metals.traverse(o => { if (o.isMesh) targets.push(o); });
    const hits = this.raycaster.intersectObjects(targets);

    if (!hits.length) return;

    this.measurePoints.push(hits[0].point.clone());

    // Show point marker
    const marker = new THREE.Mesh(
      new THREE.SphereGeometry(0.2, 8, 8),
      new THREE.MeshBasicMaterial({ color:0xffaa00 })
    );
    marker.position.copy(hits[0].point);
    marker.userData.isMeasure = true;
    this.scene.add(marker);

    if (this.measurePoints.length === 2) {
      const [p1, p2] = this.measurePoints;
      const dist = p1.distanceTo(p2) / SCALE;

      // Draw dashed line
      const geo = new THREE.BufferGeometry().setFromPoints([p1, p2]);
      const mat = new THREE.LineDashedMaterial({
        color:0xffaa00, dashSize:0.2, gapSize:0.1, linewidth:2
      });
      this.measureLine = new THREE.Line(geo, mat);
      this.measureLine.computeLineDistances();
      this.measureLine.userData.isMeasure = true;
      this.scene.add(this.measureLine);

      // Distance label
      const mid = p1.clone().add(p2).multiplyScalar(0.5);
      this.measureLabel = this._spriteLabel(
        `${dist.toFixed(1)} \u00c5`, '#ffaa00', '#ffaa00', 140);
      this.measureLabel.scale.set(5, 1.5, 1);
      this.measureLabel.position.copy(mid.clone().add(new THREE.Vector3(0, 1, 0)));
      this.measureLabel.userData.isMeasure = true;
      this.scene.add(this.measureLabel);

      this.measurePoints = [];
    }
  }

  toggleMeasure() {
    this.measureMode = !this.measureMode;
    this.measurePoints = [];
    this.renderer.domElement.style.cursor = this.measureMode ? 'crosshair' : 'grab';
    return this.measureMode;
  }

  clearMeasurements() {
    const toRemove = [];
    this.scene.traverse(o => { if (o.userData.isMeasure) toRemove.push(o); });
    for (const o of toRemove) this.scene.remove(o);
    this.measureLine = null;
    this.measureLabel = null;
    this.measurePoints = [];
  }

  // ==================== RENDER LOOP ====================

  _startRenderLoop() {
    this.clock = new THREE.Clock();
    const tick = () => {
      const dt=this.clock.getDelta(), t=this.clock.getElapsedTime();
      if (this.autoRotate) { this.cameraTheta+=dt*0.06*(this.rotateSpeed??1); this._updateCamera(); }

      // Mutation marker animation
      this.groups.mutations.children.forEach((g,i) => {
        if (!g.userData?.halo) return;
        const phase=t*1.2+i*0.6;
        g.userData.halo.scale.setScalar(1+Math.sin(phase)*0.15+(g.userData.pulseBoost||0)*0.4);
        g.userData.halo.material.opacity=0.1+Math.abs(Math.sin(phase))*0.08+(g.userData.pulseBoost||0)*0.25;
        g.userData.ring.lookAt(this.camera.position);
        g.userData.ring.rotateZ(t*0.4);
        if (g.userData.pulseBoost>0) g.userData.pulseBoost=Math.max(0,g.userData.pulseBoost-dt*0.8);
        const tgtS = g.userData.mut === this._emphMut ? 1.7 : 1;
        const cs = g.scale.x + (tgtS - g.scale.x) * Math.min(1, dt * 10);
        g.scale.setScalar(cs);
      });

      if (this._resMarker && this._resMarker.visible) {
        this._resMarker.userData.ring.lookAt(this.camera.position);
        this._resMarker.userData.ring.scale.setScalar(1 + Math.sin(t * 4) * 0.18);
      }
      this.activeSiteLight.intensity=0.5+Math.sin(t*2.5)*0.15;
      this.groups.protein.scale.setScalar(1+Math.sin(t*0.4)*0.003);

      // Exo↔Pol switching animation
      if (this.exoPolPlaying && this.exoPolCurve && this.exoPolSphere) {
        this.exoPolT += dt * 0.3;
        if (this.exoPolT > 2) this.exoPolT -= 2;
        const progress = this.exoPolT <= 1 ? this.exoPolT : 2 - this.exoPolT;
        this.exoPolSphere.position.copy(this.exoPolCurve.getPoint(progress));
      }

      // Label overlap culling
      this._updateLabelVisibility();

      // FPS counter
      this._fps.tick(t);

      this.renderer.render(this.scene,this.camera);
      requestAnimationFrame(tick);
    };
    tick();
  }

  // ==================== PUBLIC API ====================

  setToggle(name,on) {
    const map = { dna:'dna', metals:'metals', dntp:'dntp', labels:'labels',
      mutations:'mutations', surface:'surface', hbonds:'hbonds', exoPolPath:'exoPolPath',
      sheetHbonds:'sheetHbonds', zincMotifs:'zincMotifs', fesCluster:'fesCluster',
      exoSite:'exoSite', accessory:'accessory', pcna:'pcna' };
    if (name==='rotate') { this.autoRotate=on; return; }
    if (name==='exoPolPath') { this.exoPolPlaying=on; }
    const g = this.groups[map[name]]; if (g) g.visible=on;
  }

  highlightDomain(key, active) {
    this._highlightedDomain = active ? key : null;

    for (const dk of DOMAIN_ORDER) {
      const d = this.domainMeshes[dk];
      if (!d?.ribbons) continue;

      const isSelected = (dk === key && active);
      const hasSelection = active; // some domain is selected

      for (const mesh of d.ribbons) {
        if (isSelected) {
          // Highlighted domain: bright glow
          mesh.material.emissiveIntensity = 0.55;
          mesh.material.opacity = 1.0;
          mesh.material.transparent = false;
        } else if (hasSelection) {
          // Other domains: dim and transparent
          mesh.material.emissiveIntensity = 0.02;
          mesh.material.opacity = 0.15;
          mesh.material.transparent = true;
          mesh.material.depthWrite = false;
        } else {
          // No selection: restore defaults
          mesh.material.emissiveIntensity = this.colorMode === 'default' ? 0.10 : 0.12;
          mesh.material.opacity = 1.0;
          mesh.material.transparent = false;
          mesh.material.depthWrite = true;
        }
        mesh.material.needsUpdate = true;
      }
    }
  }

  focusMutation(mutId) {
    const marker=this.mutMarkers[mutId]; if (!marker) return;
    animateCamera(this, { th: -1.2, ph: 0, r: 35 });
    marker.userData.pulseBoost=1;
  }

  focusDomain(key) {
    const dom = DOMAINS[key];
    if (!dom) return;
    const domCenter = new THREE.Vector3(...dom.center);
    animateCamera(this, { r: 35, tx: domCenter.x, ty: domCenter.y, tz: domCenter.z });
  }

  resetView() {
    animateCamera(this, { th: -0.3, ph: 0.15, r: 55, tx: 0, ty: 0, tz: 0 }, 600);
  }

  // ==================== COLOR MODES ====================

  setColorMode(mode) {
    this.colorMode = mode;
    for (const key of DOMAIN_ORDER) {
      const dm = this.domainMeshes[key];
      if (!dm?.ribbons) continue;

      for (const mesh of dm.ribbons) {
        const color = PAL[DOMAINS[key].color];
        mesh.material.color.copy(color);
        mesh.material.emissive.copy(color);
        mesh.material.emissiveIntensity = 0.10;
      }
    }
  }

  setColorScheme(scheme) {
    if (!COLOR_SCHEMES[scheme]) return;
    this.colorScheme = scheme;
    const s = COLOR_SCHEMES[scheme];
    for (const key of DOMAIN_ORDER) {
      PAL[DOMAINS[key].color].setHex(s[key]);
    }
    if (this.colorMode === 'default') this.setColorMode('default');
  }

  // ==================== SCREENSHOT ====================

  takeScreenshot() {
    downloadScreenshot(this.renderer, this.scene, this.camera, 'POLE-structure');
  }

  setFPSCallback(cb) { this._fps.setCallback(cb); }
  setEmphasisMutation(id) { this._emphMut = id || null; }
}
