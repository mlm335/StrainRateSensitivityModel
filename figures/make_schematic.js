/* make_schematic.js ---------------------------------------------------------
   Editable PowerPoint schematic of the obstacles a dislocation meets, from the
   lattice up to the microstructure.  Everything is a native shape, so every
   line, circle and label can be moved in PowerPoint.

   node make_schematic.js
--------------------------------------------------------------------------- */
const pptxgen = require("pptxgenjs");
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";                      // 13.33 x 7.5 in

const NAVY = "21295C", DEEP = "065A82", TEAL = "1C7293";
const PART = "B85042", SOL = "C98A16", GREY = "5A6472", LIGHT = "EEF2F6", LINE = "C9D3DC";
const FAINT = "9AA7B4";
const BODY = "Calibri", HEAD = "Cambria";

/* ------------------------------------------------------------- helpers */
function txt(s, t, x, y, w, h, o = {}) {
  s.addText(t, Object.assign({
    x, y, w, h, isTextBox: true, margin: 0, fontFace: BODY, fontSize: 11,
    color: GREY, valign: "middle",
  }, o));
}
function line(s, x, y, w, h, o = {}) {
  s.addShape(pres.ShapeType.line, { x, y, w, h, line: Object.assign({ color: TEAL, width: 3 }, o) });
}
function dot(s, cx, cy, d, color) {
  s.addShape(pres.ShapeType.ellipse, { x: cx - d / 2, y: cy - d / 2, w: d, h: d,
    fill: { color }, line: { color, width: 1 } });
}
function panel(s, x, y, w, h, letter, title) {
  s.addShape(pres.ShapeType.roundRect, { x, y, w, h, rectRadius: 0.06,
    fill: { color: "FFFFFF" }, line: { color: LINE, width: 1 },
    shadow: { type: "outer", blur: 6, offset: 1, angle: 90, color: FAINT, opacity: 0.22 } });
  txt(s, [{ text: letter + "   ", options: { bold: true, color: NAVY } },
          { text: title, options: { color: NAVY } }],
      x + 0.24, y + 0.10, w - 0.4, 0.32, { fontFace: HEAD, fontSize: 14 });
}
function caption(s, t, x, y, w) {
  txt(s, t, x, y, w, 0.34, { fontSize: 10.5, italic: true, color: GREY, valign: "top" });
}

/* ==================================================== slide 1: four panels */
const s1 = pres.addSlide();
s1.background = { color: "FFFFFF" };
txt(s1, "What a dislocation runs into", 0.50, 0.20, 9.0, 0.50,
    { fontFace: HEAD, fontSize: 30, bold: true, color: NAVY });
txt(s1, "from the lattice itself, to individual obstacles, to the structure that holds them",
    0.52, 0.66, 11.5, 0.30, { fontSize: 13 });

const PW = 6.05, PH = 2.92, PAX = 0.50, PBX = 6.78, PAY = 1.08, PCY = 4.22;

/* ---------------------------------------------------------------- (a) */
panel(s1, PAX, PAY, PW, PH, "(a)", "Lattice friction and kink pairs");

// corrugated Peierls surface: eight valleys drawn in perspective
const VX = 0.98, VY = 3.46, DX = 0.26, DY = 0.185, VW = 3.15, NV = 8;
for (let i = 0; i < NV; i++) {
  line(s1, VX + i * DX, VY - i * DY, VW, 0, { color: LIGHT, width: 7 });
  line(s1, VX + i * DX, VY - i * DY, VW, 0,
    { color: i === 2 || i === 3 ? FAINT : LINE, width: i === 2 || i === 3 ? 1.6 : 1.2 });
}
// close the surface so the valleys read as one corrugated plane
line(s1, VX, VY, (NV - 1) * DX, -(NV - 1) * DY, { color: LINE, width: 1.2 });
line(s1, VX + VW, VY, (NV - 1) * DX, -(NV - 1) * DY, { color: LINE, width: 1.2 });
txt(s1, "Peierls valleys", 4.45, 1.52, 1.9, 0.26, { fontSize: 10.5, italic: true, align: "right" });
s1.addShape(pres.ShapeType.line, { x: 5.10, y: 1.78, w: 0.35, h: 0.18,
  line: { color: FAINT, width: 1, endArrowType: "triangle" } });

// the dislocation line: in one valley, with a kink pair thrown into the next
const yv4 = VY - 4 * DY, yv5 = VY - 5 * DY;
line(s1, 2.35, yv4, 0.98, 0);
line(s1, 3.33, yv4, DX, -DY);
line(s1, 3.33 + DX, yv5, 0.86, 0);
line(s1, 4.45, yv5, DX, DY);
line(s1, 4.45 + DX, yv4, 0.39, 0);
txt(s1, "kink pair", 2.62, 1.76, 1.10, 0.26,
    { fontSize: 10.5, bold: true, color: TEAL, align: "center" });
s1.addShape(pres.ShapeType.line, { x: 3.30, y: 2.02, w: 0.36, h: 0.52,
  line: { color: TEAL, width: 1 } });
txt(s1, "dislocation line", 1.00, 2.60, 1.30, 0.24,
    { fontSize: 10, color: TEAL, align: "left", fill: { color: "FFFFFF" } });
s1.addShape(pres.ShapeType.line, { x: 2.24, y: 2.72, w: 0.14, h: 0,
  line: { color: TEAL, width: 1 } });
caption(s1, "the line sits in a valley and moves forward by throwing a kink pair into the next one",
    PAX + 0.24, PAY + 2.46, PW - 0.5);

/* ---------------------------------------------------------------- (b) */
panel(s1, PBX, PAY, PW, PH, "(b)", "Waiting at discrete obstacles");

// applied stress
s1.addShape(pres.ShapeType.line, { x: 7.26, y: 3.22, w: 0, h: -0.72,
  line: { color: NAVY, width: 3, endArrowType: "triangle" } });
txt(s1, "applied\nstress", 6.82, 2.02, 0.90, 0.44, { fontSize: 10, color: NAVY, align: "center", valign: "top" });

// two obstacles with the segment pinned on them and bowing between
const yb = 2.86;
dot(s1, 8.90, yb, 0.52, PART);
dot(s1, 11.20, yb, 0.52, PART);
line(s1, 7.75, yb, 1.15, 0);
line(s1, 8.90, yb, 0.42, -0.30);
line(s1, 9.32, yb - 0.30, 0.38, -0.12);
line(s1, 9.70, yb - 0.42, 0.72, 0);
line(s1, 10.42, yb - 0.42, 0.38, 0.12);
line(s1, 10.80, yb - 0.30, 0.40, 0.30);
line(s1, 11.20, yb, 1.30, 0, { endArrowType: "triangle" });
txt(s1, "obstacles", 8.36, 3.16, 1.1, 0.24, { fontSize: 10, color: PART, align: "center" });
txt(s1, "bowed segment", 9.34, 2.10, 1.45, 0.24,
    { fontSize: 10, color: TEAL, align: "center", fill: { color: "FFFFFF" } });
s1.addShape(pres.ShapeType.line, { x: 10.06, y: 2.34, w: 0, h: 0.10, line: { color: TEAL, width: 1 } });

// the two escapes
const cl = (x, y, w, h) => s1.addShape(pres.ShapeType.line, { x, y, w, h,
  line: { color: DEEP, width: 2, dashType: "dash" } });
cl(11.02, yb, 0.14, -0.50); cl(11.16, yb - 0.50, 0.22, 0); cl(11.38, yb - 0.50, 0.14, 0.50);
txt(s1, "break away,\nor climb over", 11.60, 2.06, 1.30, 0.46,
    { fontSize: 10, color: DEEP, align: "center", valign: "top" });
s1.addShape(pres.ShapeType.line, { x: 11.90, y: 2.54, w: -0.52, h: 0.22,
  line: { color: DEEP, width: 1 } });
caption(s1, "the segment stops, bows out, and waits until it gets past the obstacle",
    PBX + 0.24, PAY + 2.46, PW - 0.5);

/* ---------------------------------------------------------------- (c) */
panel(s1, PAX, PCY, PW, PH, "(c)", "The structure that stores dislocations");

// a grain with sub-boundaries inside it
s1.addShape(pres.ShapeType.hexagon, { x: 0.82, y: 4.86, w: 2.05, h: 1.72,
  fill: { color: LIGHT }, line: { color: NAVY, width: 1.6 } });
[[1.34, 4.98, 0.30, 1.48], [1.92, 4.92, 0.34, 1.60], [2.42, 5.06, 0.22, 1.32]]
  .forEach(([x, y, w, h]) => s1.addShape(pres.ShapeType.line, { x, y, w, h,
    line: { color: DEEP, width: 1.1 } }));
[[1.06, 5.40, 1.70, 0.10], [1.02, 5.90, 1.78, 0.08]]
  .forEach(([x, y, w, h]) => s1.addShape(pres.ShapeType.line, { x, y, w, h,
    line: { color: DEEP, width: 1.1, dashType: "sysDot" } }));
txt(s1, "grain", 0.86, 6.58, 1.0, 0.24, { fontSize: 10, color: NAVY, align: "center" });
txt(s1, "sub-boundaries\n(sub-grain, lath)", 1.86, 6.56, 1.7, 0.42,
    { fontSize: 10, color: DEEP, align: "center", valign: "top" });

// expansion into a small volume
[[2.88, 5.00, 0.62, -0.28], [2.88, 6.34, 0.62, 0.10]].forEach(([x, y, w, h]) =>
  s1.addShape(pres.ShapeType.line, { x, y, w, h, line: { color: FAINT, width: 1, dashType: "dash" } }));
s1.addShape(pres.ShapeType.rect, { x: 3.54, y: 4.72, w: 2.62, h: 1.72,
  fill: { color: "FFFFFF" }, line: { color: GREY, width: 1.2 } });
s1.addShape(pres.ShapeType.line, { x: 3.54, y: 4.72, w: 0.26, h: -0.24, line: { color: GREY, width: 1 } });
s1.addShape(pres.ShapeType.line, { x: 6.16, y: 4.72, w: 0.26, h: -0.24, line: { color: GREY, width: 1 } });
s1.addShape(pres.ShapeType.line, { x: 3.80, y: 4.48, w: 2.62, h: 0, line: { color: GREY, width: 1 } });
s1.addShape(pres.ShapeType.line, { x: 6.42, y: 4.48, w: 0, h: 1.72, line: { color: GREY, width: 1 } });
s1.addShape(pres.ShapeType.line, { x: 6.16, y: 6.44, w: 0.26, h: -0.24, line: { color: GREY, width: 1 } });

// dislocation lines and two particle populations inside the volume
[[3.72, 5.02, 0.74, 0.30], [4.46, 5.32, 0.60, -0.22], [3.80, 5.86, 0.66, -0.24],
 [4.62, 6.02, 0.72, 0.22], [5.24, 5.10, 0.56, 0.34]]
  .forEach(([x, y, w, h]) => line(s1, x, y, w, h, { width: 1.8 }));
[[4.30, 5.72], [5.36, 5.62], [3.94, 6.20], [5.68, 6.08]].forEach(([x, y]) => dot(s1, x, y, 0.20, PART));
[[4.82, 5.06], [5.02, 5.84], [4.16, 5.34], [5.60, 5.32], [4.58, 6.26]].forEach(([x, y]) => dot(s1, x, y, 0.10, SOL));
txt(s1, "dislocations", 3.60, 6.50, 1.2, 0.24, { fontSize: 9.5, color: TEAL });
txt(s1, "coarse and fine particles", 4.62, 6.50, 1.70, 0.24, { fontSize: 9.5, color: PART, align: "right" });

/* ---------------------------------------------------------------- (d) */
panel(s1, PBX, PCY, PW, PH, "(d)", "Obstacle spacing sets the strength");

const yd = 5.66;
line(s1, 7.10, yd, 5.42, 0, { width: 2.6 });
const obs = [8.10, 9.86, 11.86];
obs.forEach(x => dot(s1, x, yd, 0.42, PART));
dot(s1, 10.70, yd - 0.72, 0.24, PART);
dot(s1, 9.00, yd + 0.74, 0.24, PART);

// spacing dimensions
const dim = (x1, x2, y, label) => {
  s1.addShape(pres.ShapeType.line, { x: x1, y, w: x2 - x1, h: 0,
    line: { color: GREY, width: 1, beginArrowType: "triangle", endArrowType: "triangle" } });
  s1.addShape(pres.ShapeType.line, { x: x1, y: y - 0.10, w: 0, h: 0.20, line: { color: GREY, width: 0.8 } });
  s1.addShape(pres.ShapeType.line, { x: x2, y: y - 0.10, w: 0, h: 0.20, line: { color: GREY, width: 0.8 } });
  txt(s1, label, (x1 + x2) / 2 - 0.70, y - 0.40, 1.40, 0.26,
      { fontSize: 10.5, align: "center", fill: { color: "FFFFFF" } });
};
dim(8.10, 9.86, 5.06, "spacing  L");
dim(9.86, 11.86, 6.34, "spacing  L′");

// obstacle radius
s1.addShape(pres.ShapeType.line, { x: 11.86, y: yd, w: 0.21, h: 0,
  line: { color: NAVY, width: 1.2, endArrowType: "triangle" } });
txt(s1, "radius  r", 11.52, yd - 0.44, 1.10, 0.24,
    { fontSize: 10.5, color: NAVY, align: "center", fill: { color: "FFFFFF" } });
txt(s1, "dislocation line", 7.06, yd - 0.34, 1.3, 0.24, { fontSize: 10, color: TEAL });
caption(s1, "closer and larger obstacles mean a shorter free length and a harder material",
    PBX + 0.24, PCY + 2.46, PW - 0.5);

s1.addNotes("Four-panel schematic: lattice friction and kink pairs, waiting at obstacles, the "
  + "structure that stores dislocations, and obstacle spacing. Every element is a native, editable shape.");

/* ============================================ slide 2: travel and wait */
function card(s, x, y, w, h, title, text, accent) {
  s.addShape(pres.ShapeType.roundRect, { x, y, w, h, rectRadius: 0.08, fill: { color: "FFFFFF" },
    line: { color: LINE, width: 1 },
    shadow: { type: "outer", blur: 6, offset: 1, angle: 90, color: FAINT, opacity: 0.25 } });
  s.addShape(pres.ShapeType.ellipse, { x: x + 0.22, y: y + 0.22, w: 0.34, h: 0.34,
    fill: { color: accent }, line: { color: accent } });
  txt(s, title, x + 0.66, y + 0.16, w - 0.85, 0.45, { fontFace: HEAD, fontSize: 15, bold: true, color: NAVY });
  txt(s, text, x + 0.24, y + 0.68, w - 0.48, h - 0.85,
      { fontSize: 12, valign: "top", lineSpacingMultiple: 1.1 });
}

const s2 = pres.addSlide();
s2.background = { color: "FFFFFF" };
txt(s2, "A dislocation travels, then waits", 0.55, 0.30, 9.5, 0.60,
    { fontFace: HEAD, fontSize: 34, bold: true, color: NAVY });
txt(s2, "every term in the model is a statement about how long a segment spends travelling and how long it spends held up",
    0.57, 0.92, 11.5, 0.40, { fontSize: 13.5 });

s2.addShape(pres.ShapeType.rect, { x: 0.55, y: 1.62, w: 12.2, h: 2.72,
  fill: { color: LIGHT }, line: { color: LINE, width: 1 } });
txt(s2, "slip plane", 0.72, 3.98, 1.4, 0.28, { fontSize: 11, italic: true });
s2.addShape(pres.ShapeType.line, { x: 0.90, y: 1.92, w: 11.5, h: 0,
  line: { color: GREY, width: 1, dashType: "dash", endArrowType: "triangle" } });
txt(s2, "glide direction", 10.70, 1.58, 1.8, 0.28, { fontSize: 11, align: "right" });

const yL = 3.12;
const seg = (x, y, w, h) => line(s2, x, y, w, h, { width: 3.5 });
seg(0.95, yL, 1.25, 0); seg(2.20, yL, 0, -0.42); seg(2.20, yL - 0.42, 0.80, 0);
seg(3.00, yL - 0.42, 0, 0.42); seg(3.00, yL, 0.95, 0);
txt(s2, "kink pair", 2.02, 2.26, 1.2, 0.26, { fontSize: 11, bold: true, color: TEAL, align: "center" });
s2.addShape(pres.ShapeType.line, { x: 2.60, y: 2.54, w: 0, h: 0.14, line: { color: TEAL, width: 1 } });

dot(s2, 4.88, 3.14, 0.55, PART);
seg(3.95, yL, 0.42, 0); seg(4.37, yL, 0.35, -0.40); seg(4.72, yL - 0.40, 0.32, 0);
seg(5.04, yL - 0.40, 0.35, 0.40); seg(5.39, yL, 0.50, 0);
txt(s2, "held at a particle, the line bows out", 3.85, 3.40, 2.2, 0.50,
    { fontSize: 11, color: NAVY, align: "center", valign: "top" });

s2.addShape(pres.ShapeType.line, { x: 7.05, y: 2.70, w: 0.46, h: 0.78, line: { color: NAVY, width: 2.5 } });
s2.addShape(pres.ShapeType.line, { x: 7.05, y: 3.48, w: 0.46, h: -0.78, line: { color: NAVY, width: 2.5 } });
seg(6.05, yL, 0.98, 0); seg(7.55, yL, 0.95, 0);
txt(s2, "forest junction", 6.53, 2.32, 1.5, 0.26, { fontSize: 11, color: NAVY, align: "center" });

[0, 1, 2, 3, 4].forEach(i => dot(s2, 8.76 + i * 0.22, yL - 0.16, 0.12, SOL));
seg(8.50, yL, 1.30, 0);
txt(s2, "solutes collect while it waits", 8.15, 2.34, 2.2, 0.26, { fontSize: 11, color: SOL, align: "center" });

dot(s2, 10.90, 3.14, 0.55, PART);
seg(9.80, yL, 0.66, 0);
const climb = (x, y, w, h) => s2.addShape(pres.ShapeType.line, { x, y, w, h,
  line: { color: DEEP, width: 2.5, dashType: "dash" } });
climb(10.46, yL, 0.28, -0.58); climb(10.74, yL - 0.58, 0.44, 0); climb(11.18, yL - 0.58, 0.28, 0.58);
line(s2, 11.46, yL, 1.12, 0, { width: 3.5, endArrowType: "triangle" });
txt(s2, "climbs over or breaks away,\nand travels again", 9.95, 3.40, 2.8, 0.50,
    { fontSize: 11, color: NAVY, align: "center", valign: "top" });

card(s2, 0.55, 4.62, 3.90, 2.35, "Travelling",
  "Between obstacles the segment glides by nucleating kink pairs on screw dislocations and sweeping them sideways. "
  + "This is slow and hard when it is cold, and almost free once the barrier is thermally overcome.", TEAL);
card(s2, 4.72, 4.62, 3.90, 2.35, "Waiting",
  "At forest junctions, particles and boundaries the line stops and bows out. "
  + "The longer it waits, the more solute atoms reach it and the stronger the pinning becomes.", PART);
card(s2, 8.88, 4.62, 3.87, 2.35, "Getting free",
  "Two escapes compete at every obstacle: thermal activation past it, or climb over it. "
  + "Whichever happens first releases the segment, and the travelling starts again.", DEEP);
s2.addNotes("The travel-and-wait picture behind the flow rule. All shapes are editable.");

/* ============================================ slide 3: mechanism cards */
const s3 = pres.addSlide();
s3.background = { color: "FFFFFF" };
txt(s3, "The mechanisms in play", 0.55, 0.30, 9.0, 0.60,
    { fontFace: HEAD, fontSize: 34, bold: true, color: NAVY });
txt(s3, "what sets the glide speed, what holds the line up, and what changes the structure while it deforms",
    0.57, 0.92, 11.5, 0.40, { fontSize: 13.5 });

const CW = 3.02, CH = 2.45, CX0 = 0.55, CY0 = 1.58, GX = 3.17, GY = 2.66;
const items = [
  ["Kink-pair glide", "How fast the segment travels through the lattice between obstacles. Sets the low-temperature strength and its strong rate and temperature dependence.", TEAL],
  ["Forest junctions", "Other dislocations threading the slip plane. They act as a back stress the moving line has to push against everywhere.", NAVY],
  ["Boundaries", "Grain, sub-grain and lath boundaries. They store dislocations and shorten the distance a segment can travel before it is stopped.", DEEP],
  ["Solute atmospheres", "Interstitial solutes diffusing to a waiting segment. Slower deformation means more time to age, which is dynamic strain aging.", SOL],
  ["Particle detachment", "Thermal activation lets a pinned segment break away from a second-phase particle once it has bowed out far enough.", PART],
  ["Local climb", "The other escape: the line climbs over the particle. It needs diffusion, so it takes over when it is hot and the loading is slow.", "8E4162"],
  ["Recovery", "Stored dislocations annihilate: by cross slip as the crystal deforms, and by climb, which keeps going whether or not it deforms.", "2C6E63"],
  ["Diffusional flow", "Matter moving through the grains and along the boundaries rather than dislocations moving. A separate carrier that only matters at very low stress.", GREY],
];
items.forEach((it, i) => card(s3, CX0 + (i % 4) * GX, CY0 + Math.floor(i / 4) * GY,
  CW, CH, it[0], it[1], it[2]));
s3.addNotes("One card per mechanism. Colours match the schematics on the earlier slides.");

pres.writeFile({ fileName: "dislocation_schematic.pptx" }).then(f => console.log("wrote", f));
