declare let angular: any;
declare let FS: any;
declare let Mousetrap: any;
declare let tinycolor: any;
declare let tinygradient: any;
function assert(condition: boolean, message = "") {
  if (!condition) {
    throw new Error(message);
  }
}
function forEachValue(o: any, fn: (v: any) => void) {
  for (let n in o) {
    fn(o[n]);
  }
}
function toPercent(v: number) {
  return (v * 100).toFixed(2) + "%";
}
function withCommas(v: number) {
  return v.toLocaleString();
}
class Y4MFile {
  constructor(public size: Size, public buffer: Uint8Array, public frames: Y4MFrame []) {
    // ...
  }
}
class Y4MFrame {
  constructor(public y: number, public cb: number, public cr: number) {
    // ...
  }
}
interface Window {
  Module: any;
}

interface CanvasRenderingContext2D {
  mozImageSmoothingEnabled: boolean;
  imageSmoothingEnabled;
}

interface Math {
  imul: (a: number, b: number) => number;
}

let colors = [
  "#E85EBE", "#009BFF", "#00FF00", "#0000FF", "#FF0000", "#01FFFE", "#FFA6FE",
  "#FFDB66", "#006401", "#010067", "#95003A", "#007DB5", "#FF00F6", "#FFEEE8",
  "#774D00", "#90FB92", "#0076FF", "#D5FF00", "#FF937E", "#6A826C", "#FF029D",
  "#FE8900", "#7A4782", "#7E2DD2", "#85A900", "#FF0056", "#A42400", "#00AE7E",
  "#683D3B", "#BDC6FF", "#263400", "#BDD393", "#00B917", "#9E008E", "#001544",
  "#C28C9F", "#FF74A3", "#01D0FF", "#004754", "#E56FFE", "#788231", "#0E4CA1",
  "#91D0CB", "#BE9970", "#968AE8", "#BB8800", "#43002C", "#DEFF74", "#00FFC6",
  "#FFE502", "#620E00", "#008F9C", "#98FF52", "#7544B1", "#B500FF", "#00FF78",
  "#FF6E41", "#005F39", "#6B6882", "#5FAD4E", "#A75740", "#A5FFD2", "#FFB167"
];

function hexToRGB(hex: string, alpha: number = 0) {
  let r = parseInt(hex.slice(1,3), 16),
      g = parseInt(hex.slice(3,5), 16),
      b = parseInt(hex.slice(5,7), 16),
      a = "";
  if (alpha) {
    a = ", 1";
  }
  return "rgb(" + r + ", " + g + ", " + b + a + ")";
}

function getLineOffset(lineWidth: number) {
  return lineWidth % 2 == 0 ? 0 : 0.5;
}

const mi_block_size_log2 = 3;

function alignPowerOfTwo(value: number, n: number) {
  return ((value) + ((1 << n) - 1)) & ~((1 << n) - 1);
}

function tileOffset(i: number, rowsOrCols: number, tileRowsOrColsLog2: number) {
  let sbRowsOrCols = alignPowerOfTwo(rowsOrCols, mi_block_size_log2) >> mi_block_size_log2;
  let offset = ((i * sbRowsOrCols) >> tileRowsOrColsLog2) << mi_block_size_log2;
  return Math.min(offset, rowsOrCols);
}

enum AOMAnalyzerPredictionMode {
  DC_PRED 		= 0,   // Average of above and left pixels
  V_PRED 			= 1,   // Vertical
  H_PRED 			= 2,   // Horizontal
  D45_PRED 		= 3,   // Directional 45  deg = round(arctan(1/1) * 180/pi)
  D135_PRED 	= 4,   // Directional 135 deg = 180 - 45
  D117_PRED 	= 5,   // Directional 117 deg = 180 - 63
  D153_PRED 	= 6,   // Directional 153 deg = 180 - 27
  D207_PRED 	= 7,   // Directional 207 deg = 180 + 27
  D63_PRED 		= 8,   // Directional 63  deg = round(arctan(2/1) * 180/pi)
  TM_PRED 		= 9,   // True-motion
  NEARESTMV 	= 10,
  NEARMV 			= 11,
  ZEROMV 			= 12,
  NEWMV 			= 13
}

enum AOMAnalyzerBlockSize {
  BLOCK_4X4   = 0,
  BLOCK_4X8   = 1,
  BLOCK_8X4   = 2,
  BLOCK_8X8   = 3,
  BLOCK_8X16  = 4,
  BLOCK_16X8  = 5,
  BLOCK_16X16 = 6,
  BLOCK_16X32 = 7,
  BLOCK_32X16 = 8,
  BLOCK_32X32 = 9,
  BLOCK_32X64 = 10,
  BLOCK_64X32 = 11,
  BLOCK_64X64 = 12
}

enum Property {
  GET_CODEC_BUILD_CONFIG,
  GET_CLPF_STRENGTH,
  GET_DERING_LEVEL
}

enum MIProperty {
  GET_MI_MV,
  GET_MI_MV_REFERENCE_FRAME,
  GET_MI_MODE,
  GET_MI_SKIP,
  GET_MI_BLOCK_SIZE,
  GET_MI_TRANSFORM_TYPE,
  GET_MI_TRANSFORM_SIZE,
  GET_MI_DERING_GAIN,
  GET_MI_BITS,
  GET_MI_AC_Y_DEQUANT,
  GET_MI_DC_Y_DEQUANT,
  GET_MI_AC_UV_DEQUANT,
  GET_MI_DC_UV_DEQUANT
}

enum AccountingProperty {
  GET_ACCCOUNTING_SYMBOL_COUNT,
  GET_ACCCOUNTING_SYMBOL_NAME,
  GET_ACCCOUNTING_SYMBOL_BITS,
  GET_ACCCOUNTING_SYMBOL_CONTEXT_X,
  GET_ACCCOUNTING_SYMBOL_CONTEXT_Y
}

enum AOMAnalyzerTransformMode {
  ONLY_4X4       = 0,
  ALLOW_8X8      = 1,
  ALLOW_16X16    = 2,
  ALLOW_32X32    = 3,
  TX_MODE_SELECT = 4
}

enum AOMAnalyzerTransformType {
  DCT_DCT        = 0,
  ADST_DCT       = 1,
  DCT_ADST       = 2,
  ADST_ADST      = 3
}

enum AOMAnalyzerTransformSize {
  TX_4X4         = 0,
  TX_8X8         = 1,
  TX_16X16       = 2,
  TX_32X32       = 3
}

/**
 * Maps AOMAnalyzerTransformSize enum to [w, h] log2 pairs.
 */
const TRANSFORM_SIZES = [
  [2, 2],
  [3, 3],
  [4, 4],
  [5, 5]
];

/**
 * Maps AOMAnalyzerBlockSize enum to [w, h] log2 pairs.
 */
const BLOCK_SIZES = [
  [2, 2],
  [2, 3],
  [3, 2],
  [3, 3],
  [3, 4],
  [4, 3],
  [4, 4],
  [4, 5],
  [5, 4],
  [5, 5],
  [5, 6],
  [6, 5],
  [6, 6]
];

class AccountingSymbol {
  constructor(public name: string, public bits: number, public samples: number, public x: number, public y: number) {
    // ...
  }
}

type AccountingSymbolMap = { [name: string]: AccountingSymbol };

class Accounting {
  symbols: AccountingSymbol [] = null;
  frameSymbols: AccountingSymbolMap = Object.create(null);
  constructor(symbols: AccountingSymbol [] = []) {
    this.symbols = symbols;
  }
  createFrameSymbols() {
    this.frameSymbols = Accounting.flatten(this.symbols);
  }
  createBlockSymbols(mi: Vector) {
    return Accounting.flatten(this.symbols.filter(symbol => {
      return symbol.x === mi.x && symbol.y === mi.y;
    }));
  }

  static flatten(sybmols: AccountingSymbol []): AccountingSymbolMap {
    let map = Object.create(null);
    sybmols.forEach(symbol => {
      let s = map[symbol.name];
      if (!s) {
        s = map[symbol.name] = new AccountingSymbol(symbol.name, 0, 0, symbol.x, symbol.y);
      }
      s.bits += symbol.bits;
      s.samples += symbol.samples;
    });
    let ret = Object.create(null);
    let names = [];
    for (let name in map) names.push(name);
    // Sort by bits.
    names.sort((a, b) => map[b].bits - map[a].bits);
    names.forEach(name => {
      ret[name] = map[name];
    });
    return ret;
  }
}

interface AOMInternal {
  _read_frame (): number;
  _get_plane (pli: number): number;
  _get_plane_stride (pli: number): number;
  _get_plane_width (pli: number): number;
  _get_plane_height (pli: number): number;
  _get_mi_cols_and_rows(): number;
  _get_tile_cols_and_rows_log2(): number;
  _get_frame_count(): number;
  _get_frame_width(): number;
  _get_frame_height(): number;
  _open_file(): number;

  _get_property(p: Property): number;
  _get_accounting_property(p: AccountingProperty, i: number): number;
  _get_mi_property(p: MIProperty, mi_col: number, mi_row: number, i: number): number;

  _get_predicted_plane_buffer(pli: number): number;
  _get_predicted_plane_stride(pli: number): number;

  HEAPU8: Uint8Array;
  UTF8ToString(p: number): string;
}

class AOM {
  native: AOMInternal = (<any>window).Module;
  HEAPU8: Uint8Array = this.native.HEAPU8;
  constructor () {

  }
  openFile() {
    return this.native._open_file();
  }
  readFrame() {
    return this.native._read_frame();
  }
  getPlane(pli: number): number {
    return this.native._get_plane(pli);
  }
  getPlaneStride(pli: number): number {
    return this.native._get_plane_stride(pli);
  }
  getPlaneWidth(pli: number): number {
    return this.native._get_plane_width(pli);
  }
  getPlaneHeight(pli: number): number {
    return this.native._get_plane_height(pli);
  }
  getProperty(p: Property) {
    return this.native._get_property(p);
  }
  getAccountingProperty(p: AccountingProperty, i: number = 0) {
    return this.native._get_accounting_property(p, i);
  }
  getString(i: number): string {
    return this.native.UTF8ToString(i);
  }
  getMIProperty(p: MIProperty, mi_col: number, mi_row: number, i: number = 0) {
    return this.native._get_mi_property(p, mi_col, mi_row, i);
  }
  getPredictedPlaneBuffer(pli: number): number {
    return this.native._get_predicted_plane_buffer(pli);
  }
  getPredictedPlaneStride(pli: number): number {
    return this.native._get_predicted_plane_stride(pli);
  }
  getFrameCount(): number {
    return this.native._get_frame_count();
  }
  getFrameSize(): Size {
    return new Size(this.native._get_frame_width(), this.native._get_frame_height());
  }
  getMIGridSize(): GridSize {
    let v = this.native._get_mi_cols_and_rows();
    let cols = v >> 16;
    let rows = this.native._get_mi_cols_and_rows() & 0xFF;
    return new GridSize(cols, rows);
  }
  getTileGridSizeLog2(): GridSize {
    let v = this.native._get_tile_cols_and_rows_log2();
    let cols = v >> 16;
    let rows = v & 0xFF;
    return new GridSize(cols, rows);
  }
}

class ErrorMetrics {
  constructor(public tss: number, public mse: number) {
    // ...
  }
  toString() {
    return "tss: " + this.tss + ", mse: " + this.mseToString();
  }
  mseToString() {
    return this.mse.toFixed(4);
  }
}

class Size {
  constructor(public w: number, public h: number) {
    // ...
  }
  clone() {
    return new Size(this.w, this.h);
  }
  equals(other: Size) {
    return this.w == other.w || this.h == other.h;
  }
  area(): number {
    return this.w * this.h;
  }
  multiplyScalar(scalar: number) {
		if (isFinite(scalar)) {
			this.w *= scalar;
			this.h *= scalar;
		} else {
			this.w = 0;
			this.h = 0;
		}
		return this;
	}
}

class Rectangle {
  constructor (public x: number, public y: number, public w: number, public h: number) {
    // ...
  }
  static createRectangleCenteredAtPoint(v: Vector, w: number, h: number) {
    return new Rectangle(v.x - w / 2, v.y - h / 2, w, h);
  }
  static createRectangleFromSize(size: Size) {
    return new Rectangle(0, 0, size.w, size.h);
  }
  set(x: number, y: number, w: number, h: number) {
    this.x = x;
    this.y = y;
    this.w = w;
    this.h = h;
    return this;
  }
  getCenter(): Vector {
    return new Vector(this.x + this.w / 2, this.y + this.h / 2);
  }
  clone(): Rectangle {
    return new Rectangle(this.x, this.y, this.w, this.h);
  }
  multiplyScalar(scalar: number) {
		if (isFinite(scalar)) {
      this.x *= scalar;
			this.y *= scalar;
			this.w *= scalar;
			this.h *= scalar;
		} else {
      this.x = 0;
			this.y = 0;
			this.w = 0;
			this.h = 0;
		}
		return this;
	}
}

class GridSize {
  constructor (public cols: number, public rows: number) {
    // ...
  }
}

/**
 * MI block coordinates in 8x8 units, or sub coordinates in 4x4 units in case of 4x4, 4x8, 8x8 blocks.
 */
class MICoordinates {
  c: number;
  r: number;
  constructor (c: number, r: number) {
    this.c = c;
    this.r = r;
  }
  set(c: number, r: number) {
    this.c = c;
    this.r = r;
    return this;
  }
}

class Vector {
  x: number;
  y: number;
  constructor (x: number, y: number) {
    this.x = x;
    this.y = y;
  }
  set(x: number, y: number) {
    this.x = x;
    this.y = y;
    return this;
  }
  lerp(v: Vector, alpha: number) {
    this.x += (v.x - this.x) * alpha;
		this.y += (v.y - this.y) * alpha;
    return this;
  }
  clone(): Vector {
    return new Vector(this.x, this.y);
  }
  lengthSq() {
		return this.x * this.x + this.y * this.y;
	}
	length() {
		return Math.sqrt(this.x * this.x + this.y * this.y);
	}
  normalize () {
		return this.divideScalar(this.length());
	}
  multiplyScalar(scalar) {
		if (isFinite(scalar)) {
			this.x *= scalar;
			this.y *= scalar;
		} else {
			this.x = 0;
			this.y = 0;
		}
		return this;
	}
	divide(v) {
		this.x /= v.x;
		this.y /= v.y;
		return this;
	}
	divideScalar(scalar) {
		return this.multiplyScalar(1 / scalar);
	}
  snap() {
    // TODO: Snap to nearest pixel
		this.x = this.x | 0;
    this.y = this.y | 0;
    return this;
	}
  sub(v: Vector): Vector {
		this.x -= v.x;
		this.y -= v.y;
		return this;
	}
  add(v: Vector): Vector {
		this.x += v.x;
		this.y += v.y;
		return this;
	}
  clampLength(min: number, max: number) {
		let length = this.length();
		this.multiplyScalar(Math.max(min, Math.min(max, length)) / length);
		return this;
	}
  toString(): string {
    return this.x + "," + this.y;
  }
}

function drawSplit(ctx, x, y, dx, dy) {
  ctx.beginPath();
  ctx.save();
  ctx.moveTo(x, y);
  ctx.lineTo(x + dx, y);
  ctx.moveTo(x, y);
  ctx.lineTo(x, y + dy);
  ctx.restore();
  ctx.closePath();
  ctx.stroke();
}

function drawVector(ctx: CanvasRenderingContext2D, a: Vector, b: Vector) {
  ctx.beginPath();
  ctx.moveTo(a.x, a.y);
  ctx.lineTo(b.x, b.y);
  ctx.closePath();
  ctx.stroke();
  return;
}

function drawLine(ctx: CanvasRenderingContext2D, x, y, dx, dy) {
  ctx.beginPath();
  ctx.moveTo(x, y);
  ctx.lineTo(x + dx, y + dy);
  ctx.closePath();
  ctx.stroke();
}

interface BlockVisitor {
  (size: AOMAnalyzerBlockSize | AOMAnalyzerTransformSize, coordinates: MICoordinates, subCoordinates: MICoordinates, bounds: Rectangle): void;
}

enum BlockVisitorMode {
  /** Visit all block paritions. */
  BLOCK,
  /** Visit all block transfrom paritions. */
  TRANSFORM,
  /** Visit just the super blocks. */
  SUPER_BLOCK
}

interface Decoder {
  description: string;
  path: string;
}

const BLOCK_SIZE = 8;

class AppCtrl {
  aom: AOM = null;

  decoders = {
    default: {
      description: "Default",
      path: "bin/decoder.js"
    },
    dering: {
      description: "Deringing",
      path: "bin/dering-decoder.js"
    }
  };
  selectedDecoder: string;
  frameSize: Size = new Size(128, 128);
  tileGridSize: GridSize = new GridSize(0, 0);
  fileSize: number = 0;
  fileBytes: Uint8Array;
  ratio: number = 1;
  scale: number = 1;

  options = {
    showY: {
      key: "y",
      description: "Y",
      detail: "Display Y image plane.",
      updatesImage: true,
      default: true,
      value: undefined
    },
    showU: {
      key: "u",
      description: "U",
      detail: "Display U image plane.",
      updatesImage: true,
      default: true,
      value: undefined
    },
    showV: {
      key: "v",
      description: "V",
      detail: "Display V image plane.",
      updatesImage: true,
      default: true,
      value: undefined
    },
    showOriginalImage: {
      key: "w",
      description: "Original Image",
      detail: "Display loaded .y4m file.",
      updatesImage: true,
      default: false,
      disabled: true,
      value: undefined
    },
    showDecodedImage: {
      key: "i",
      description: "Decoded Image",
      detail: "Display decoded image.",
      updatesImage: true,
      default: true,
      value: undefined
    },
    showPredictedImage: {
      key: "p",
      description: "Predicted Image",
      detail: "Display the predicted image, or the residual if the decoded image is displayed.",
      updatesImage: true,
      default: false,
      value: undefined
    },
    showSuperBlockGrid: {
      key: "g",
      description: "SB Grid",
      detail: "Display the 64x64 super block grid.",
      default: false,
      value: undefined
    },
    showTileGrid: {
      key: "l",
      description: "Tile Grid",
      detail: "Display tile grid.",
      default: false,
      value: undefined
    },
    showTransformSplit: {
      key: "t",
      description: "Transform Grid",
      detail: "Display transform blocks.",
      default: false,
      value: undefined
    },
    showBlockSplit: {
      key: "s",
      description: "Split Grid",
      detail: "Display block partitions.",
      default: false,
      value: undefined
    },
    showDering: {
      key: "d",
      description: "Dering",
      detail: "Display blocks where the deringing filter is applied.",
      default: false,
      value: undefined
    },
    showMotionVectors: {
      key: "m",
      description: "Motion Vectors",
      detail: "Display motion vectors, darker colors represent longer vectors.",
      default: false,
      value: undefined
    },
    showMode: {
      key: "o",
      description: "Mode",
      detail: "Display prediction modes.",
      default: false,
      value: undefined
    },
    showBits: {
      key: "b",
      description: "Bits",
      detail: "Display bits.",
      default: false,
      value: undefined
    },
    showSkip: {
      key: "k",
      description: "Skip",
      detail: "Display skip flags.",
      default: false,
      value: undefined
    },
    showInspector: {
      key: "2",
      description: "Inspector",
      detail: "Display block and frame details.",
      default: window.innerWidth > 1024,
      value: undefined
    },
    zoomLock: {
      key: "z",
      description: "Cursor Zoom Lock",
      detail: "Locks zoom at the current mouse position.",
      default: false,
      value: undefined,
      doNotShare: true
    }
  };

  frameNumber: number = -1;
  y4mFile: Y4MFile;

  progressValue = 0;
  progressMode = "determinate";

  get isPlaying() {
    return !!this.playInterval;
  }

  // If not zero, we are playing frames.
  playInterval: number = 0;

  container: HTMLDivElement;
  displayCanvas: HTMLCanvasElement;
  overlayCanvas: HTMLCanvasElement;
  zoomCanvas: HTMLCanvasElement;
  chartCanvas: HTMLCanvasElement;

  displayContext: CanvasRenderingContext2D = null;
  overlayContext: CanvasRenderingContext2D = null;
  zoomContext: CanvasRenderingContext2D = null;
  zoomWidth = 512;
  zoomLevel = 16;

  chartContext: CanvasRenderingContext2D = null;
  chartSeriesHeight = 16;
  mousePosition: Vector = new Vector(0, 0);
  imageData: ImageData = null;

  frameCanvas: HTMLCanvasElement;
  frameContext: CanvasRenderingContext2D = null;

  compositionCanvas: HTMLCanvasElement;
  compositionContext: CanvasRenderingContext2D = null;

  lastDecodeFrameTime: number = 0;

  colorOptions = {
    modeColor: {
      description: "Mode Color",
      detail: undefined,
      updatesImage: true,
      default: "hsl(79, 20%, 50%)",
      value: undefined
    },
    bitsColor: {
      description: "Bits Color",
      detail: undefined,
      updatesImage: true,
      default: "#9400D3",
      value: undefined
    },
    mv0Color: {
      description: "MV 0 Color",
      detail: undefined,
      updatesImage: true,
      default: "blue",
      value: undefined
    },
    mv1Color: {
      description: "MV 1 Color",
      detail: undefined,
      updatesImage: true,
      default: "red",
      value: undefined
    },
    skipColor: {
      description: "Skip Color",
      detail: undefined,
      updatesImage: true,
      default: "hsla(326, 53%, 42%, 0.2)",
      value: undefined
    },
    gridColor: {
      description: "Grid Color",
      detail: undefined,
      updatesImage: true,
      default: "rgba(55,55,55,1)",
      value: undefined
    },
    tileGridColor: {
      description: "Tile Color",
      detail: undefined,
      updatesImage: true,
      default: "hsl(0, 100%, 69%)",
      value: undefined
    },
    splitColor: {
      description: "Split Color",
      detail: undefined,
      updatesImage: true,
      default: "rgba(33,33,33,1)",
      value: undefined
    },
    transformColor: {
      description: "Transform Color",
      detail: undefined,
      updatesImage: true,
      default: "rgba(255,0,0,1)",
      value: undefined
    },
    crosshairColor: {
      description: "Crosshair Color",
      detail: undefined,
      updatesImage: true,
      default: "#ffffff",
      value: undefined
    }
  };
  installOptions() {
    for (let k in this.options) {
      this.options[k].value = this.options[k].default;
    }
    for (let k in this.colorOptions) {
      this.colorOptions[k].value = this.colorOptions[k].default;
    }
  }

  loadOptions() {
    let parameters = getUrlParameters();
    for (let name in this.options) {
      if (name in parameters) {
        // TODO: Check for boolean here..
        this.options[name].value = parameters[name] == "true";
      }
    }
    this.loadColorOptions();
  }

  saveColorsOptions(force: boolean = false) {
    for (let name in this.colorOptions) {
      let option = this.colorOptions[name];
      if (force || option.value != option.default) {
        localStorage[name] = option.value;
      }
    }
  }

  loadColorOptions() {
    for (let name in this.colorOptions) {
      if (name in localStorage) {
        this.colorOptions[name].value = localStorage[name];
      }
    }
  }

  crosshairLineWidth = 2;
  gridLineWidth = 3;
  tileGridLineWidth = 5;
  splitLineWidth = 1;
  transformLineWidth = 1;
  modeLineWidth = 2;

  sharingLink = "";

  $scope: any;
  $interval: any;
  $mdSidenav: any;

  uiFrameProperties: any;
  createUIFrameProperties() {
    let self = this;
    this.uiFrameProperties = {
      frameNumber: {
        description: "Frame Number",
        get value() {
          return self.frameNumber;
        }
      },
      scale: {
        description: "Scale",
        get value() {
          return self.scale;
        }
      },
      size: {
        description: "Frame Size",
        get value() {
          let s = self.frameSize;
          return s.w + " x " + s.h;
        }
      },
      tiles: {
        description: "Tiles",
        get value() {
          let s = self.tileGridSize;
          return s.cols + " x " + s.rows;
        }
      },
      fileSize: {
        description: "File Size",
        get value() {
          return self.fileSize;
        }
      },
      decodeTime: {
        description: "Decode Time",
        get value() {
          return self.lastDecodeFrameTime.toFixed(2) + "ms";
        }
      },
      frameError: {
        description: "Frame Error",
        get value() {
          let error = self.getFrameError();
          return error ? error.mseToString() : "N/A";
        }
      },
      clpfStrength: {
        description: "CLPF Strength",
        get value() {
          return self.aom.getProperty(Property.GET_CLPF_STRENGTH);
        }
      },
      deringLevel: {
        description: "Dering Level",
        get value() {
          return self.aom.getProperty(Property.GET_DERING_LEVEL);
        }
      }
    };
  }

  uiBlockProperties: any;
  createUIBlockProperties() {
    let self = this;
    function withMIUnderMouse(fn) {
      if (!self.aom) return;
      let mi = self.getMIUnderMouse();
      return fn(mi);
    }
    this.uiBlockProperties = {
      blockSize: {
        description: "Block Size",
        get value() {
          return withMIUnderMouse(mi => {
            return AOMAnalyzerBlockSize[self.aom.getMIProperty(MIProperty.GET_MI_BLOCK_SIZE, mi.x, mi.y)];
          });
        }
      },
      blockSkip: {
        description: "Block Skip",
        get value() {
          return withMIUnderMouse(mi => {
            self.aom.getMIProperty(MIProperty.GET_MI_SKIP, mi.x, mi.y);
          });
        }
      },
      predictionMode: {
        description: "Prediction Mode",
        get value() {
          return withMIUnderMouse(mi => {
            return AOMAnalyzerPredictionMode[self.aom.getMIProperty(MIProperty.GET_MI_MODE, mi.x, mi.y)];
          });
        }
      },
      deringGain: {
        description: "Dering Gain",
        get value() {
          return withMIUnderMouse(mi => {
            return String(self.aom.getMIProperty(MIProperty.GET_MI_DERING_GAIN, mi.x, mi.y));
          });
        }
      },
      motionVector0: {
        description: "Motion Vector 0",
        get value() {
          return withMIUnderMouse(mi => {
            var v = self.getMotionVector(mi.x, mi.y, 0);
            return v.toString() + " " + v.length().toFixed(2);
          });
        },
        get color() {
          return self.colorOptions.mv0Color.value;
        }
      },
      motionVector1: {
        description: "Motion Vector 1",
        get value() {
          return withMIUnderMouse(mi => {
            var v = self.getMotionVector(mi.x, mi.y, 1);
            return v.toString() + " " + v.length().toFixed(2);
          });
        },
        get color() {
          return self.colorOptions.mv1Color.value;
        }
      },
      referenceFrames: {
        description: "Reference Frames",
        get value() {
          return withMIUnderMouse(mi => {
            return self.aom.getMIProperty(MIProperty.GET_MI_MV_REFERENCE_FRAME, mi.x, mi.y, 0) + ", " +
                   self.aom.getMIProperty(MIProperty.GET_MI_MV_REFERENCE_FRAME, mi.x, mi.y, 1);
          });
        }
      },
      transformType: {
        description: "Transform Type",
        get value() {
          return withMIUnderMouse(mi => {
            return AOMAnalyzerTransformType[self.aom.getMIProperty(MIProperty.GET_MI_TRANSFORM_TYPE, mi.x, mi.y)];
          });
        }
      },
      transformSize: {
        description: "Transform Size",
        get value() {
          return withMIUnderMouse(mi => {
            return AOMAnalyzerTransformSize[self.aom.getMIProperty(MIProperty.GET_MI_TRANSFORM_SIZE, mi.x, mi.y)];
          });
        }
      },
      bits: {
        description: "Bits",
        get value() {
          return withMIUnderMouse(mi => {
            return self.aom.getMIProperty(MIProperty.GET_MI_BITS, mi.x, mi.y);
          });
        },
        get color() {
          return self.colorOptions.bitsColor.value;
        }
      },
      blockError: {
        description: "Block Error",
        get value() {
          return withMIUnderMouse(mi => {
            let error = self.getMIError(mi);
            return error ? error.toString() : "N/A";
          });
        }
      },
      yDequant: {
        description: "Y Dequant",
        get value() {
          return withMIUnderMouse(mi => {
            return self.aom.getMIProperty(MIProperty.GET_MI_DC_Y_DEQUANT, mi.x, mi.y) + "/" +
                   self.aom.getMIProperty(MIProperty.GET_MI_AC_Y_DEQUANT, mi.x, mi.y);
          });
        }
      },
      uvDequant: {
        description: "UV Dequant",
        get value() {
          return withMIUnderMouse(mi => {
            return self.aom.getMIProperty(MIProperty.GET_MI_DC_UV_DEQUANT, mi.x, mi.y) + "/" +
                   self.aom.getMIProperty(MIProperty.GET_MI_AC_UV_DEQUANT, mi.x, mi.y);
          });
        }
      }
    };
  }

  uiAccountingFrameProperties: any;
  uiAccountingBlockProperties: any;

  createUIAccountingProperties() {
    let self = this;
    this.uiAccountingFrameProperties = { };
    this.uiAccountingBlockProperties = { };
  }

  constructor($scope, $interval, $mdSidenav) {
    let self = this;
    this.$scope = $scope;
    this.$mdSidenav = $mdSidenav;
    // File input types don't have angular bindings, so we need set the
    // event handler on the scope object.
    $scope.fileInputNameChanged = function() {
      let input = <any>event.target;
      let reader = new FileReader();
      reader.onload = function() {
        let buffer = reader.result;
        self.openFileBytes(new Uint8Array(buffer));
        self.playFrameAsync(1, () => {
          self.drawFrame();
        });
      };
      reader.readAsArrayBuffer(input.files[0]);
    };
    $scope.y4mFileInputNameChanged = function() {
      let input = <any>event.target;
      let reader = new FileReader();
      reader.onload = function() {
        let buffer = reader.result;
        let y4mFile = self.loadY4MBytes(new Uint8Array(buffer));
        if (!y4mFile.size.equals(self.frameSize)) {
          alert("Y4M file frame size doesn't match current frame size.")
          return;
        }
        self.y4mFile = y4mFile;
        self.drawImages();
        // We can now show the image.
        self.options.showOriginalImage.disabled = false;
      };
      reader.readAsArrayBuffer(input.files[0]);
    };

    this.$interval = $interval;
    this.ratio = window.devicePixelRatio || 1;
    this.scale = this.ratio;
    this.scale = 1;

    this.container = <HTMLDivElement>document.getElementById("container");

    this.displayCanvas = <HTMLCanvasElement>document.getElementById("display");
    this.displayContext = this.displayCanvas.getContext("2d");

    this.overlayCanvas = <HTMLCanvasElement>document.getElementById("overlay");
    this.overlayContext = this.overlayCanvas.getContext("2d");

    this.zoomCanvas = <HTMLCanvasElement>document.getElementById("zoom");
    this.zoomContext = this.zoomCanvas.getContext("2d");
    this.zoomContext.mozImageSmoothingEnabled = false;
    this.zoomContext.imageSmoothingEnabled = false;

    this.chartCanvas = <HTMLCanvasElement>document.getElementById("chart");
    this.chartContext = this.chartCanvas.getContext("2d");
    this.chartContext.mozImageSmoothingEnabled = false;
    this.chartContext.imageSmoothingEnabled = false;

    window.addEventListener("resize", this.onResize.bind(this));

    this.overlayCanvas.addEventListener("mousemove", this.onMouseMove.bind(this));
    this.overlayCanvas.addEventListener("mousedown", this.onMouseDown.bind(this));

    this.frameCanvas = document.createElement("canvas");
    this.frameContext = this.frameCanvas.getContext("2d");

    this.compositionCanvas = document.createElement("canvas");
    this.compositionContext = this.compositionCanvas.getContext("2d");

    this.installOptions();

    this.loadOptions();

    this.installKeyboardShortcuts();

    this.createUIFrameProperties();
    this.createUIBlockProperties();
    this.createUIAccountingProperties();
    this.loadOptions();

    let parameters = getUrlParameters();
    let frames = parseInt(parameters.frameNumber) || 1;

    this.loadDecoder("default", () => {
      this.aom = new AOM();
      let file = "media/default.ivf";
      this.openFile(file, () => {
        this.playFrameAsync(frames, () => {
          this.drawFrame();
        })
      });
    });
  }

  installKeyboardShortcuts() {
    Mousetrap.bind(['.'], () => {
      this.playFrame();
      this.drawFrame();
      this.uiApply();
    });

    Mousetrap.bind(['space'], (e) => {
      this.uiPlayPause();
      e.preventDefault();
    });

    Mousetrap.bind([']'], () => {
      this.uiZoom(2);
      this.uiApply();
    });

    Mousetrap.bind(['['], () => {
      this.uiZoom(1 / 2);
      this.uiApply();
    });

    Mousetrap.bind(['x'], (e) => {
      this.uiResetLayers();
      e.preventDefault();
    });

    Mousetrap.bind(['r'], (e) => {
      this.uiReload();
      e.preventDefault();
    });

    Mousetrap.bind(['tab'], (e) => {
      this.uiToggleMenu();
      e.preventDefault();
    });

    let self = this;
    function toggle(name, event) {
      let option = this.options[name];
      option.value = !option.value;
      if (option.updatesImage) {
        self.drawImages();
      }
      self.drawMain();
      self.options.showInspector.value && self.drawInfo();
      self.uiApply();
      event.preventDefault();
    }

    let installedKeys = {};
    for (let name in this.options) {
      let option = this.options[name];
      if (option.key) {
        if (installedKeys[option.key]) {
          console.error("Key: " + option.key + " for " + option.description  + ", is already mapped to " + installedKeys[option.key].description);
        }
        installedKeys[option.key] = option;
        Mousetrap.bind([option.key], toggle.bind(this, name));
      }
    }
  }

  onMouseDown(event: MouseEvent) {
    this.handleMouseEvent(event);
    // TODO: Think about this some more.
    // this.options.zoomLock.value = !this.options.zoomLock.value;
    this.uiApply();
  }

  onMouseMove(event: MouseEvent) {
    this.handleMouseEvent(event);
  }

  onResize() {
    // Is browser zooming?
    let ratio = window.devicePixelRatio || 1;
    if (this.ratio !== ratio) {
      this.ratio = ratio;
      this.resetCanvases();
      this.drawFrame();
    }
  }

  handleMouseEvent(event: MouseEvent) {
    if (this.options.zoomLock.value) {
      return;
    }
    function getMousePosition(canvas: HTMLCanvasElement, event: MouseEvent) {
      let rect = canvas.getBoundingClientRect();
      return new Vector(
        event.clientX - rect.left,
        event.clientY - rect.top
      );
    }
    this.mousePosition = getMousePosition(this.overlayCanvas, event);
    this.updateBlockAccounting();
    this.options.showInspector.value && this.drawInfo();
    this.uiApply();
  }

  getMIBlockSize(c: number, r: number, miMinBlockSize: AOMAnalyzerBlockSize = AOMAnalyzerBlockSize.BLOCK_4X4): Size {
    let miBlockSize = this.aom.getMIProperty(MIProperty.GET_MI_BLOCK_SIZE, c, r);
    if (miBlockSize >= BLOCK_SIZES.length) {
      // TODO: This should not happen, figure out what is going on.
      return new Size(0, 0);
    }
    if (miBlockSize < miMinBlockSize) {
      miBlockSize = miMinBlockSize;
    }
    let w = 1 << BLOCK_SIZES[miBlockSize][0];
    let h = 1 << BLOCK_SIZES[miBlockSize][1];
    return new Size(w, h);
  }

  /**
   * Gets the coordinates of the MI block under the mousedown.
   */
  getMIUnderMouse(): Vector {
    let v = this.mousePosition;
    let c = (v.x / this.scale) >> 3;
    let r = (v.y / this.scale) >> 3;
    return this.getMI(c, r);
  }

  /**
   * Get the coordinates of the parent MI block if this block is
   * not an 8x8 block.
   */
  getMI(c: number, r: number): Vector {
    let blockSize = this.getMIBlockSize(c, r);
    c = c & ~((blockSize.w - 1) >> 3);
    r = r & ~((blockSize.h - 1) >> 3);
    return new Vector(c, r);
  }

  getMIBits(c: number, r: number): number {
    let mi = this.getMI(c, r);
    return this.aom.getMIProperty(MIProperty.GET_MI_BITS, mi.x, mi.y);
  }

  loadDecoder(decoder: string, next: () => any) {
    this.selectedDecoder = decoder;
    let s = document.createElement('script');
    s.onload = next;
    s.setAttribute('src', this.decoders[decoder].path);
    document.body.appendChild(s);
  }

  uiChangeDecoder() {
    alert("NYI");
    this.loadDecoder(this.selectedDecoder, () => {
      // ...
    });
  }

  loadY4MBytes(buffer: Uint8Array): Y4MFile {
    return this.parseY4MBytes(buffer);
  }

  parseY4MBytes(buffer: Uint8Array): Y4MFile {
    let header;
    let eol = "\n".charCodeAt(0);
    let offset = 0;
    while (offset < buffer.length) {
      if (buffer[offset++] == eol) {
        header = String.fromCharCode.apply(null, (buffer.subarray(0, offset - 1)));
        break;
      }
    }
    let parameters = header.split(" ");
    let size = new Size(0, 0);
    let colorSpace = "420";
    for (let i = 0; i < parameters.length; i++) {
      let parameter = parameters[i];
      if (parameter[0] == "W") {
        size.w = parseInt(parameter.substring(1));
      } else if (parameter[0] == "H") {
        size.h = parseInt(parameter.substring(1));
      } else if (parameter[0] == "C") {
        colorSpace = parameter.substring(1);
      }
    }
    if (colorSpace != "420" && colorSpace != "420jpeg") {
      console.error("Unsupported color space: " + colorSpace);
      return;
    }
    let y = size.w * size.h;
    let cb = ((size.w + 1) >> 1) * ((size.h + 1) >> 1);
    let cr = cb;
    let frameLength = y + cb + cr;
    let frames: Y4MFrame [] = [];
    while (offset < buffer.length) {
      let start = offset;
      while (offset < buffer.length) {
        if (buffer[offset++] == eol) {
          break;
        }
      }
      let frameHeader = String.fromCharCode.apply(null, (buffer.subarray(start, offset - 1)));
      if (frameHeader != "FRAME") {
        console.error("Cannot parse frame: ");
        return;
      }
      frames.push(new Y4MFrame(offset, offset + y, offset + y + cb));
      offset += frameLength;
    }
    return new Y4MFile(size, buffer, frames);
  }

  openFileBytes(buffer: Uint8Array) {
    this.fileBytes = buffer;
    this.fileSize = buffer.length;
    FS.writeFile("/tmp/input.ivf", buffer, { encoding: "binary" });
    this.aom.openFile();
    this.frameNumber = -1;
    this.frameSize = this.aom.getFrameSize();
    this.resetCanvases();
  }

  openFile(path: string, next: () => any = null) {
    let fileName = path.replace(/^.*[\\\/]/, '')
    document.title = fileName;
    this.downloadFile(path, (buffer: Uint8Array) => {
      this.openFileBytes(buffer);
      next();
    });
  }

  resetCanvases() {
    this.frameCanvas.width = this.compositionCanvas.width = this.frameSize.w;
		this.frameCanvas.height = this.compositionCanvas.height = this.frameSize.h;

    this.imageData = this.frameContext.createImageData(this.frameSize.w, this.frameSize.h);

    this.container.style.width = (this.frameSize.w * this.scale) + "px";
		this.container.style.height = (this.frameSize.h * this.scale) + "px";

    this.displayCanvas.style.width = (this.frameSize.w * this.scale) + "px";
		this.displayCanvas.style.height = (this.frameSize.h * this.scale) + "px";
    this.displayCanvas.width = this.frameSize.w * this.scale * this.ratio;
		this.displayCanvas.height = this.frameSize.h * this.scale * this.ratio;

    this.overlayCanvas.style.width = (this.frameSize.w * this.scale) + "px";
		this.overlayCanvas.style.height = (this.frameSize.h * this.scale) + "px";
    this.overlayCanvas.width = this.frameSize.w * this.scale * this.ratio;
		this.overlayCanvas.height = this.frameSize.h * this.scale * this.ratio;

    this.zoomCanvas.style.width = this.zoomWidth + "px";
		this.zoomCanvas.style.height = this.zoomWidth + "px";
    this.zoomCanvas.width = this.zoomWidth * this.ratio;
		this.zoomCanvas.height = this.zoomWidth * this.ratio;
  }

  downloadFile(path: string, next: (buffer: Uint8Array) => void) {
    let xhr = new XMLHttpRequest();
    let self = this;
    self.progressMode = "determinate";
    xhr.open("GET", path, true);
    xhr.responseType = "arraybuffer";
    xhr.send();
    xhr.addEventListener("progress", (e) => {
      let progress = (e.loaded / e.total) * 100;
      this.progressValue = progress;
      this.$scope.$apply();
    });
    xhr.addEventListener("load", function () {
      if (xhr.status != 200) {
        return;
      }
      next(new Uint8Array(this.response));
    });
  }

  showFileInputDialog() {
    angular.element(document.querySelector('#fileInput'))[0].click();
  }

  showY4MFileInputDialog() {
    angular.element(document.querySelector('#y4mFileInput'))[0].click();
  }

  uiAction(name) {
    let file;
    switch (name) {
      case "open-file":
        this.showFileInputDialog();
        return;
      case "open-y4m-file":
        this.showY4MFileInputDialog();
        return;
      case "open-crosswalk":
        file = "media/crosswalk_30.ivf";
        break;
      case "open-soccer":
        file = "media/soccer_30.ivf";
        break;
      case "open-tiger":
        file = "media/tiger_30.ivf";
        break;
      case "open-tiger-60":
        file = "media/tiger_60.ivf";
        break;
    }
    this.openFile(file, () => {
      this.playFrameAsync(1, () => {
        this.drawFrame();
      })
    });
  }

  createSharingLink() {
    let url = location.protocol + '//' + location.host + location.pathname;
    let args = {
      frameNumber: this.frameNumber
    };
    for (let name in this.options) {
      let option = this.options[name];
      if (option.doNotShare) {
        continue;
      }
      // Ignore default values.
      if (option.value == option.default) {
        continue;
      }
      args[name] = option.value;
    }
    let argList = [];
    for (let arg in args) {
      argList.push(arg + "=" + encodeURIComponent(args[arg]));
    }
    let argListString = argList.join("&");
    return url + "?" + argListString;
  }

  updateSharingLink() {
    this.sharingLink = this.createSharingLink();
  }

  fileIssue(label: string = "") {
    window.open("https://github.com/mbebenita/aomanalyzer/issues/new?labels=" + label + "&body=" + encodeURIComponent(this.createSharingLink()));
  }

  uiResetLayers() {
    for (let name in this.options) {
      this.options[name].value = this.options[name].default;
    }
    this.drawFrame();
  }

  uiResetColors() {
    for (let name in this.colorOptions) {
      this.colorOptions[name].value = this.colorOptions[name].default;
    }
    this.drawFrame();
    this.saveColorsOptions(true);
  }

  uiChangeColor() {
    this.saveColorsOptions();
    this.resetCanvases();
    this.drawFrame();
  }

  uiChange() {
    this.updateSharingLink();
    this.resetCanvases();
    this.drawFrame();
  }

  playPause() {
    if (this.playInterval) {
      this.$interval.cancel(this.playInterval);
      this.playInterval = 0;
      return;
    }
    this.playInterval = this.$interval(() => {
      if (!this.playFrame()) {
        this.$interval.cancel(this.playInterval);
        this.playInterval = 0;
        return;
      }
      this.drawFrame();
    }, 1);
  }

  playFrameAsync(count: number, step: () => void = null, stop: () => void = null) {
    this.$interval(() => {
      this.playFrame();
      step && step();
      if (--count == 0) {
        stop && stop();
      }
    }, 1, count);
  }

  playFrame(count: number = 1) {
    for (let i = 0; i < count; i++) {
      let s = performance.now();
      if (this.aom.readFrame()) {
        return false;
      }
      this.processFrame();
      this.lastDecodeFrameTime = performance.now() - s;
      this.frameNumber ++;

      let tileGridSize = this.aom.getTileGridSizeLog2();
      this.tileGridSize.cols = 1 << tileGridSize.cols;
      this.tileGridSize.rows = 1 << tileGridSize.rows;
    }
    return true;
  }

  frameStatistics = {
    bits: {
      values: [],
      show: true,
      description: "Plot frame bits over time.",
      detail: "Plot the number of bits per frame over time."
    },
    errors: {
      values: [],
      show: false,
      description: "Plot frame MSE over time.",
      detail: "Plot the MSE per frame over time. A .y4m file is required."
    }
  };

  lastAccounting: Accounting = null;

  processFrame() {
    let {cols, rows} = this.aom.getMIGridSize();
    let miTotalBits = 0;
    for (let c = 0; c < cols; c++) {
      for (let r = 0; r < rows; r++) {
        miTotalBits += this.aom.getMIProperty(MIProperty.GET_MI_BITS, c, r);
      }
    }
    this.frameStatistics.bits.values.push(miTotalBits);
    this.frameStatistics.errors.values.push(this.getFrameError());

    this.lastAccounting = this.accounting;
    this.accounting = this.getAccounting();
    this.updateBlockAccounting()
    this.updateFrameAccounting();
  }

  accounting: Accounting = null;

  getAccounting(): Accounting {
    var aom = this.aom;
    var accounting = new Accounting();
    let count = aom.getAccountingProperty(AccountingProperty.GET_ACCCOUNTING_SYMBOL_COUNT);
    let nameMap = [];
    for (let i = 0; i < count; i++) {
      let nameAddress = aom.getAccountingProperty(AccountingProperty.GET_ACCCOUNTING_SYMBOL_NAME, i);
      if (nameMap[nameAddress] === undefined) {
        nameMap[nameAddress] = aom.getString(nameAddress);
      }
      let name = nameMap[nameAddress];
      let bits = aom.getAccountingProperty(AccountingProperty.GET_ACCCOUNTING_SYMBOL_BITS, i);
      let x = aom.getAccountingProperty(AccountingProperty.GET_ACCCOUNTING_SYMBOL_CONTEXT_X, i);
      let y = aom.getAccountingProperty(AccountingProperty.GET_ACCCOUNTING_SYMBOL_CONTEXT_Y, i);
      accounting.symbols.push(new AccountingSymbol(name, bits, 1, x, y));
    }
    return accounting;
  }

  updateFrameAccounting() {
    let accounting = this.accounting;
    accounting.createFrameSymbols();
    let total = 0;
    forEachValue(accounting.frameSymbols, (symbol) => {
      total += symbol.bits;
    });
    this.uiAccountingFrameProperties = { };
    for (let name in accounting.frameSymbols) {
      let symbol = accounting.frameSymbols[name];
      let lastSymbol = this.lastAccounting ? this.lastAccounting.frameSymbols[name] : null;
      this.uiAccountingFrameProperties[symbol.name] = {
        description: symbol.name,
        value: [
          withCommas(symbol.bits),
          toPercent(symbol.bits / total),
          withCommas(symbol.samples),
          lastSymbol ? " " + (symbol.bits - lastSymbol.bits) : ""
        ]
      }
    }
  }

  updateBlockAccounting() {
    let accounting = this.accounting;
    if (!accounting) {
      return;
    }
    let mi = this.getMIUnderMouse();
    let blockSymbols = accounting.createBlockSymbols(mi);
    let total = 0;
    forEachValue(blockSymbols, (symbol) => {
      total += symbol.bits;
    });
    this.uiAccountingBlockProperties = { };
    for (let name in blockSymbols) {
      let symbol = blockSymbols[name];
      assert(symbol.bits >= 0);
      this.uiAccountingBlockProperties[symbol.name] = {
        description: symbol.name,
        value: [
          withCommas(symbol.bits),
          toPercent(symbol.bits / total),
          withCommas(symbol.samples)
        ]
      }
    }
  }

  uiZoom(value: number) {
    this.scale *= value;
    this.resetCanvases();
    this.drawFrame();
  }

  uiApply() {
    this.updateSharingLink();
    this.$scope.$apply();
  }

  uiReload() {
    // TODO: Don't reload entire file.
    this.fileBytes && this.openFileBytes(this.fileBytes);
    this.resetCanvases();
    this.playFrame();
    this.drawFrame();
  }

  uiPreviousFrame() {
    this.uiApply();
  }

  uiPlayPause() {
    this.playPause();
  }

  uiNextFrame() {
    this.playFrame();
    this.drawFrame();
  }

  uiToggleMenu() {
    this.$mdSidenav("left").toggle();
  }

  uiToggleInspector() {
    this.options.showInspector.value = !this.options.showInspector.value;
  }

  drawInfo() {
    let mousePosition = this.mousePosition.clone().divideScalar(this.scale).snap();
    let src = Rectangle.createRectangleCenteredAtPoint(mousePosition, 64, 64);
    let dst = new Rectangle(0, 0, this.zoomWidth * this.ratio, this.zoomWidth * this.ratio);

    this.zoomContext.clearRect(0, 0, dst.w, dst.h);
    if (this.options.showOriginalImage.value ||
        this.options.showDecodedImage.value ||
        this.options.showPredictedImage.value) {
      this.zoomContext.mozImageSmoothingEnabled = false;
      this.zoomContext.imageSmoothingEnabled = false;
      this.zoomContext.clearRect(dst.x, dst.y, dst.w, dst.h);
      this.zoomContext.drawImage(this.compositionCanvas,
        src.x, src.y, src.w, src.h,
        dst.x, dst.y, dst.w, dst.h);
    }
    this.drawLayers(this.zoomContext, src, dst);
    this.drawCrosshair(this.zoomContext, mousePosition, dst);
    this.drawFrameStatistics();
  }

  drawFrameStatistics() {
    let count = 0;
    for (let name in this.frameStatistics) {
      if (this.frameStatistics[name].show && this.frameStatistics[name].values.length) {
        count ++;
      }
    }
    let dst = new Rectangle(0, 0, this.zoomWidth, this.chartSeriesHeight * count);
    this.chartCanvas.style.width = dst.w + "px";
		this.chartCanvas.style.height = dst.h + "px";
    this.chartCanvas.width = dst.w * this.ratio;
		this.chartCanvas.height = dst.h * this.ratio;
    let ctx = this.chartContext;
    ctx.clearRect(0, 0, dst.w, dst.h);

    let barHMax = this.chartSeriesHeight * this.ratio;
    let barHPadding = 0;
    let barW = 6 * this.ratio;
    let barWPadding = 2 * this.ratio;
    let barWTotal = barW + barWPadding;
    let maxFrames = dst.w * this.ratio / barWTotal | 0;

    function drawSeries(data, count) {
      var min = Number.MAX_VALUE;
      var max = Number.MIN_VALUE;
      for (var i = 0; i < data.length; i++) {
        min = Math.min(min, data[i]);
        max = Math.max(max, data[i]);
      }
      data = data.slice(Math.max(data.length - count, 0));
      for (var i = 0; i < data.length; i++) {
        let h = (data[i] / max) * barHMax | 0;
        ctx.fillRect(i * barWTotal, barHMax - h, barW, h);
      }
      ctx.translate(0, barHMax + barHPadding);
    }

    let bits = this.frameStatistics.bits;
    if (bits.show) {
      ctx.fillStyle = this.colorOptions.bitsColor.value;
      drawSeries(bits.values, maxFrames);
    }

    let errors = this.frameStatistics.errors;
    if (errors.show) {
      ctx.fillStyle = "#E91E63";
      drawSeries(errors.values.map(x => x ? x.mse : 0), maxFrames);
    }
  }

  drawCrosshair(ctx: CanvasRenderingContext2D, center: Vector, dst: Rectangle) {
    ctx.save();
    ctx.globalCompositeOperation = "difference";
    ctx.lineWidth = this.crosshairLineWidth;
    ctx.strokeStyle = this.colorOptions.crosshairColor.value;

    // Draw Lines
    ctx.beginPath();
    let lineOffset = getLineOffset(this.crosshairLineWidth);
    ctx.translate(lineOffset, lineOffset)
    ctx.moveTo(dst.x, dst.y + dst.h / 2);
    ctx.lineTo(dst.x + dst.w, dst.y + dst.h / 2);
    ctx.moveTo(dst.x + dst.w / 2, 0);
    ctx.lineTo(dst.x + dst.w / 2, dst.y + dst.h);
    ctx.closePath();
    ctx.stroke();

    // Draw Dot
    ctx.beginPath();
    ctx.arc(dst.x + dst.w / 2, dst.y + dst.h / 2, this.crosshairLineWidth * 2, 0, Math.PI * 2, true);
    ctx.closePath();
    ctx.fill();

    // Draw Text
    ctx.fillStyle = "#FFFFFF";

    ctx.lineWidth = 4;
    let textHeight = 12 * this.ratio;
    let textPadding = 4 * this.ratio;
    ctx.font = "bold " + textHeight + "px sans-serif";

    let x, y, text;
    x = dst.x + dst.w / 2 + textPadding;
    y = textHeight + textPadding;
    text = String(center.y) + " (" + (center.y / BLOCK_SIZE | 0) + ")";
    ctx.fillText(text, x, y);

    x = textPadding;
    y = dst.y + dst.h / 2 + textHeight + textPadding;
    text = String(center.x) + " (" + (center.x / BLOCK_SIZE | 0) + ")";
    ctx.fillText(text, x, y);

    ctx.restore();

  }

  drawSuperBlockGrid(ctx: CanvasRenderingContext2D, src: Rectangle, dst: Rectangle) {
    let {cols, rows} = this.aom.getMIGridSize();
    let scale = dst.w / src.w | 0;
    let scaledFrameSize = this.frameSize.clone().multiplyScalar(scale);
    ctx.save();
    ctx.globalAlpha = 1;
    let lineOffset = getLineOffset(this.gridLineWidth);
    ctx.translate(lineOffset, lineOffset);
    ctx.translate(-src.x * scale, -src.y * scale);
    ctx.beginPath();
    for (let c = 0; c <= cols; c += 8) {
      let offset = c * BLOCK_SIZE * scale;
      ctx.moveTo(offset, 0);
      ctx.lineTo(offset, scaledFrameSize.h);
    }
    for (let r = 0; r <= rows; r += 8) {
      let offset = r * BLOCK_SIZE * scale;
      ctx.moveTo(0, offset);
      ctx.lineTo(scaledFrameSize.w, offset);
    }
    ctx.closePath();
    ctx.strokeStyle = this.colorOptions.gridColor.value;
    ctx.lineWidth = this.gridLineWidth;
    ctx.stroke();
    ctx.restore();
  }

  drawTileGrid(ctx: CanvasRenderingContext2D, src: Rectangle, dst: Rectangle) {
    let {cols, rows} = this.aom.getMIGridSize();
    let {cols: tileColsLog2, rows: tileRowsLog2} = this.aom.getTileGridSizeLog2();
    let scale = dst.w / src.w | 0;
    let scaledFrameSize = this.frameSize.clone().multiplyScalar(scale);
    ctx.save();
    ctx.globalAlpha = 1;
    let lineOffset = getLineOffset(this.tileGridLineWidth);
    ctx.translate(lineOffset, lineOffset);
    ctx.translate(-src.x * scale, -src.y * scale);
    ctx.beginPath();

    for (let c = 0; c <= 1 << tileColsLog2; c ++) {
      let offset = tileOffset(c, cols, tileColsLog2) * BLOCK_SIZE * scale;
      ctx.moveTo(offset, 0);
      ctx.lineTo(offset, scaledFrameSize.h);
    }
    for (let r = 0; r <= 1 << tileRowsLog2; r ++) {
      let offset = tileOffset(r, rows, tileRowsLog2) * BLOCK_SIZE * scale;
      ctx.moveTo(0, offset);
      ctx.lineTo(scaledFrameSize.w, offset);
    }
    ctx.closePath();
    ctx.strokeStyle = this.colorOptions.tileGridColor.value;
    ctx.lineWidth = this.tileGridLineWidth;
    ctx.stroke();
    ctx.restore();
  }

  drawFrame() {
    this.drawImages();
    this.options.showInspector.value && this.drawInfo();
    this.drawMain();
  }

  clearImage() {
    this.compositionContext.clearRect(0, 0, this.frameSize.w, this.frameSize.h);
    this.displayContext.clearRect(0, 0, this.frameSize.w * this.scale * this.ratio, this.frameSize.h * this.scale * this.ratio);
  }

  drawImages() {
    this.clearImage();
    this.options.showOriginalImage.value && this.drawOriginalImage();
    if (this.options.showDecodedImage.value && this.options.showPredictedImage.value) {
      this.drawDecodedImage();
      this.drawPredictedImage("difference");
      this.invertImage();
    } else {
      this.options.showDecodedImage.value && this.drawDecodedImage();
      this.options.showPredictedImage.value && this.drawPredictedImage();
    }
  }

  drawOriginalImage(compositeOperation: string = "source-over") {
    if (!this.y4mFile) {
      return;
    }
    let file = this.y4mFile;
    let frame = file.frames[this.frameNumber];
    let Yp = frame.y;
    let Ys = file.size.w;
    let Up = frame.cb;
    let Us = file.size.w >> 1;
    let Vp = frame.cr;
    let Vs = file.size.w >> 1;
    this.drawImage(file.buffer, Yp, Ys, Up, Us, Vp, Vs, compositeOperation);
  }

  drawDecodedImage(compositeOperation: string = "source-over") {
    let Yp = this.aom.getPlane(0);
    let Ys = this.aom.getPlaneStride(0);
    let Up = this.aom.getPlane(1);
    let Us = this.aom.getPlaneStride(1);
    let Vp = this.aom.getPlane(2);
    let Vs = this.aom.getPlaneStride(2);
    this.drawImage(this.aom.HEAPU8, Yp, Ys, Up, Us, Vp, Vs, compositeOperation);
  }

  drawPredictedImage(compositeOperation: string = "source-over") {
    let Yp = this.aom.getPredictedPlaneBuffer(0);
    let Ys = this.aom.getPredictedPlaneStride(0);

    let Up = this.aom.getPredictedPlaneBuffer(1);
    let Us = this.aom.getPredictedPlaneStride(1);

    let Vp = this.aom.getPredictedPlaneBuffer(2);
    let Vs = this.aom.getPredictedPlaneStride(2);
    this.drawImage(this.aom.HEAPU8, Yp, Ys, Up, Us, Vp, Vs, compositeOperation);
  }

  invertImage() {
    this.compositionContext.globalCompositeOperation = "difference";
    this.compositionContext.fillStyle = "#FFFFFF";
    this.compositionContext.fillRect(0, 0, this.frameSize.w, this.frameSize.h);
  }

  drawImage(H: Uint8Array, Yp, Ys, Up, Us, Vp, Vs, compositeOperation: string) {
    let I = this.imageData.data;

    let w = this.frameSize.w;
    let h = this.frameSize.h;


    let showY = Yp && this.options.showY.value;
    let showU = Up && this.options.showU.value;
    let showV = Vp && this.options.showV.value;

    for (let y = 0; y < h; y++) {
      for (let x = 0; x < w; x++) {
        let index = (Math.imul(y, w) + x) << 2;

        let Y = showY ? H[Yp + Math.imul(y, Ys) + x] : 128;
        let U = showU ? H[Up + Math.imul(y >> 1, Us) + (x >> 1)] : 128;
        let V = showV ? H[Vp + Math.imul(y >> 1, Vs) + (x >> 1)] : 128;

        let bgr = YUV2RGB(Y, U, V);

        let r = (bgr >>  0) & 0xFF;
        let g = (bgr >>  8) & 0xFF;
        let b = (bgr >> 16) & 0xFF;

        I[index + 0] = r;
        I[index + 1] = g;
        I[index + 2] = b;
        I[index + 3] = 255;
      }
    }

    if (this.imageData) {
      this.frameContext.putImageData(this.imageData, 0, 0);
      this.compositionContext.globalCompositeOperation = compositeOperation;
      this.compositionContext.drawImage(this.frameCanvas, 0, 0, this.frameSize.w, this.frameSize.h);
    }

    this.drawMain();
  }

  drawMain() {
    // Draw composited image.
    this.displayContext.mozImageSmoothingEnabled = false;
    this.displayContext.imageSmoothingEnabled = false;
    let dw = this.frameSize.w * this.scale * this.ratio;
    let dh = this.frameSize.h * this.scale * this.ratio;
    this.displayContext.drawImage(this.compositionCanvas, 0, 0, dw, dh);

    // Draw Layers
    let ctx = this.overlayContext;
    let ratio = window.devicePixelRatio || 1;
    ctx.clearRect(0, 0, this.frameSize.w * this.scale * ratio, this.frameSize.h * this.scale * ratio);

    let src = Rectangle.createRectangleFromSize(this.frameSize);
    let dst = src.clone().multiplyScalar(this.scale * this.ratio);

    this.drawLayers(ctx, src, dst);
  }

  drawLayers(ctx: CanvasRenderingContext2D, src: Rectangle, dst: Rectangle) {
    this.options.showMotionVectors.value && this.drawMotionVectors(ctx, src, dst);
    this.options.showDering.value && this.drawDering(ctx, src, dst);
    this.options.showBits.value &&  this.drawBits(ctx, src, dst);
    this.options.showMode.value && this.drawMode(ctx, src, dst);
    this.options.showSkip.value && this.drawSkip(ctx, src, dst);
    this.options.showTransformSplit.value && this.drawTransformSplit(ctx, src, dst);
    this.options.showBlockSplit.value && this.drawBlockSplit(ctx, src, dst);
    this.options.showSuperBlockGrid.value && this.drawSuperBlockGrid(ctx, src, dst);
    this.options.showTileGrid.value && this.drawTileGrid(ctx, src, dst);
  }

  drawMode(ctx: CanvasRenderingContext2D, src: Rectangle, dst: Rectangle) {
    let scale = dst.w / src.w;
    ctx.save();
    ctx.lineWidth = this.modeLineWidth;
    ctx.strokeStyle = this.colorOptions.modeColor.value;
    ctx.globalAlpha = 0.5;
    let lineOffset = getLineOffset(this.modeLineWidth);
    ctx.translate(lineOffset, lineOffset);
    ctx.translate(-src.x * scale, -src.y * scale);
    let lineWidth = 1;
    ctx.lineWidth = lineWidth;
    this.visitBlocks(BlockVisitorMode.BLOCK, (blockSize, coordinates, subCoordinates, bounds) => {
      bounds.multiplyScalar(scale);
      let i = this.aom.getMIProperty(MIProperty.GET_MI_MODE, coordinates.c, coordinates.r);
      mode(i, bounds);
    });

    function mode(m: AOMAnalyzerPredictionMode, bounds: Rectangle) {
      let x = bounds.x;
      let y = bounds.y;
      let w = bounds.w;
      let h = bounds.h;
      let hw = w / 2;
      let hh = h / 2;
      switch (m) {
        case AOMAnalyzerPredictionMode.V_PRED:
          drawLine(ctx, x + hw + lineOffset, y, 0, h);
          break;
        case AOMAnalyzerPredictionMode.H_PRED:
          drawLine(ctx, x, y + hh + lineOffset, w, 0);
          break;
        case AOMAnalyzerPredictionMode.D45_PRED:
          drawLine(ctx, x, y + h, w, -h);
          break;
        case AOMAnalyzerPredictionMode.D63_PRED:
          drawLine(ctx, x, y + h, hw, -h);
          break;
        case AOMAnalyzerPredictionMode.D135_PRED:
          drawLine(ctx, x, y, w, h);
          break;
        case AOMAnalyzerPredictionMode.D117_PRED:
          drawLine(ctx, x + hw, y, hw, h);
          break;
        case AOMAnalyzerPredictionMode.D153_PRED:
          drawLine(ctx, x, y + hh, w, hh);
          break;
        case AOMAnalyzerPredictionMode.D207_PRED:
          drawLine(ctx, x, y + hh, w, -hh);
          break;
        default:
          ctx.fillStyle = colors[m];
          ctx.fillRect(x, y, w, h);
          break;
      }
    }

    ctx.restore();
  }

  drawBlockSplit(ctx: CanvasRenderingContext2D, src: Rectangle, dst: Rectangle) {
    let scale = dst.w / src.w;
    ctx.save();
    ctx.lineWidth = this.splitLineWidth;
    ctx.strokeStyle = this.colorOptions.splitColor.value;
    ctx.globalAlpha = 1;
    let lineOffset = getLineOffset(this.splitLineWidth);
    ctx.translate(lineOffset, lineOffset);
    ctx.translate(-src.x * scale, -src.y * scale);
    let lineWidth = 1;
    ctx.lineWidth = lineWidth;
    this.visitBlocks(BlockVisitorMode.BLOCK, (blockSize, coordinates, subCoordinates, bounds) => {
      bounds.multiplyScalar(scale);
      drawSplit(ctx, bounds.x, bounds.y, bounds.w, bounds.h);
    });
    ctx.restore();
  }

  drawTransformSplit(ctx: CanvasRenderingContext2D, src: Rectangle, dst: Rectangle) {
    let scale = dst.w / src.w;
    ctx.save();
    ctx.lineWidth = this.transformLineWidth;
    ctx.strokeStyle = this.colorOptions.transformColor.value;
    ctx.globalAlpha = 1;
    let lineOffset = getLineOffset(this.transformLineWidth);
    ctx.translate(lineOffset, lineOffset);
    ctx.translate(-src.x * scale, -src.y * scale);
    let lineWidth = 1;
    ctx.lineWidth = lineWidth;
    this.visitBlocks(BlockVisitorMode.TRANSFORM, (blockSize, coordinates, subCoordinates, bounds) => {
      bounds.multiplyScalar(scale);
      drawSplit(ctx, bounds.x, bounds.y, bounds.w, bounds.h);
    });
    ctx.restore();
  }

  drawFillBlock(mode: BlockVisitorMode, ctx: CanvasRenderingContext2D, src: Rectangle, dst: Rectangle, setFillStyle: (coordinates: MICoordinates, subCoordinates: MICoordinates) => boolean) {
    let scale = dst.w / src.w;
    ctx.save();
    ctx.translate(-src.x * scale, -src.y * scale);
    this.visitBlocks(mode, (blockSize, coordinates, subCoordinates, bounds) => {
      bounds.multiplyScalar(scale);
      setFillStyle(coordinates, subCoordinates) && ctx.fillRect(bounds.x, bounds.y, bounds.w, bounds.h);
    });
    ctx.restore();
  }

  drawDering(ctx: CanvasRenderingContext2D, src: Rectangle, dst: Rectangle) {
    this.drawFillBlock(BlockVisitorMode.SUPER_BLOCK, ctx, src, dst, (coordinates, subCoordinates) => {
      let i = this.aom.getMIProperty(MIProperty.GET_MI_DERING_GAIN, coordinates.c, coordinates.r);
      if (i == 0) return false;
      ctx.fillStyle = "rgba(33,33,33," + (i / 4) + ")";
      return true;
    });
  }

  drawSkip(ctx: CanvasRenderingContext2D, src: Rectangle, dst: Rectangle) {
    this.drawFillBlock(BlockVisitorMode.BLOCK, ctx, src, dst, (coordinates, subCoordinates) => {
      let i = this.aom.getMIProperty(MIProperty.GET_MI_SKIP, coordinates.c, coordinates.r);
      if (i == 0) return false;
      ctx.fillStyle = this.colorOptions.skipColor.value;
      return true;
    });
  }

  getMIBlockBitsPerPixel(c: number, r: number): number {
    // Bits are stored at the 8x8 level, even if the block is split further.
    let blockSize = this.getMIBlockSize(c, r, AOMAnalyzerBlockSize.BLOCK_8X8);
    let blockArea = blockSize.w * blockSize.h;
    let miBits = this.getMIBits(c, r);
    return miBits / blockArea;
  }

  drawBits(ctx: CanvasRenderingContext2D, src: Rectangle, dst: Rectangle) {
    let {cols, rows} = this.aom.getMIGridSize();
    let miMaxBitsPerPixel = 0;
    for (let c = 0; c < cols; c++) {
      for (let r = 0; r < rows; r++) {
        let miBitsPerPixel = this.getMIBlockBitsPerPixel(c, r);
        miMaxBitsPerPixel = Math.max(miMaxBitsPerPixel, miBitsPerPixel);
      }
    }
    let gradient = tinygradient([
      {color: tinycolor(this.colorOptions.bitsColor.value).brighten(100), pos: 0},
      {color: tinycolor(this.colorOptions.bitsColor.value), pos: 1}
    ]);
    let colorRange = gradient.rgb(32).map(x => x.toString());
    ctx.globalAlpha = 0.75;
    this.drawFillBlock(BlockVisitorMode.BLOCK, ctx, src, dst, (coordinates, subCoordinates) => {
      let miBitsPerPixel = this.getMIBlockBitsPerPixel(coordinates.c, coordinates.r);
      let color = colorRange[((miBitsPerPixel / miMaxBitsPerPixel) * (colorRange.length - 1)) | 0];
      ctx.fillStyle = color;
      return true;
    });
  }

  getMotionVector(c: number, r: number, i: number): Vector {
    let mv = this.aom.getMIProperty(MIProperty.GET_MI_MV, c, r, i);
    let y = (mv >> 16);
    let x = (((mv & 0xFFFF) << 16) >> 16);
    return new Vector(x, y);
  }

  getFrameError(): ErrorMetrics {
    let file = this.y4mFile;
    if (!this.y4mFile) {
      return;
    }
    let frame = file.frames[this.frameNumber];

    let AYp = frame.y;
    let AYs = file.size.w;
    let AH = file.buffer;

    let BYp = this.aom.getPlane(0);
    let BYs = this.aom.getPlaneStride(0);
    let BH = this.aom.HEAPU8;

    let h = this.frameSize.h;
    let w = this.frameSize.w;

    let Ap = AYp;
    let Bp = BYp;
    let error = 0;
    for (let y = 0; y < h; y++) {
      for (let x = 0; x < w; x++) {
        let d = AH[Ap + x] - BH[Bp + x];
        error += d * d;
      }
      Ap += AYs;
      Bp += BYs;
    }
    return new ErrorMetrics(error, error / this.frameSize.area());
  }

  getMIError(mi: Vector): ErrorMetrics {
    let file = this.y4mFile;
    if (!this.y4mFile) {
      return;
    }
    let frame = file.frames[this.frameNumber];
    let AYp = frame.y;
    let AYs = file.size.w;
    let AH = file.buffer;

    let BYp = this.aom.getPlane(0);
    let BYs = this.aom.getPlaneStride(0);
    let BH = this.aom.HEAPU8;
    let size = this.getMIBlockSize(mi.x, mi.y);
    let Ap = AYp + mi.y * BLOCK_SIZE * AYs + mi.x * BLOCK_SIZE;
    let Bp = BYp + mi.y * BLOCK_SIZE * BYs + mi.x * BLOCK_SIZE;
    let error = 0;
    for (let y = 0; y < size.h; y++) {
      for (let x = 0; x < size.w; x++) {
        let d = AH[Ap + x] - BH[Bp + x];
        error += d * d;
      }
      Ap += AYs;
      Bp += BYs;
    }
    return new ErrorMetrics(error, error / size.area());
  }

  visitBlocks(mode: BlockVisitorMode, visitor: BlockVisitor) {
    let {cols, rows} = this.aom.getMIGridSize();
    let s = BLOCK_SIZE;
    let coordinates = new MICoordinates(0, 0);
    let subCoordinates = new MICoordinates(0, 0);
    var bounds = new Rectangle(0, 0, 0, 0);

    if (mode === BlockVisitorMode.BLOCK) {
      // Visit blocks >= 8x8.
      for (let i = 3; i < BLOCK_SIZES.length; i++) {
        let dc = 1 << (BLOCK_SIZES[i][0] - 3);
        let dr = 1 << (BLOCK_SIZES[i][1] - 3);
        for (let c = 0; c < cols; c += dc) {
          for (let r = 0; r < rows; r += dr) {
            let size = this.aom.getMIProperty(MIProperty.GET_MI_BLOCK_SIZE, c, r);
            let w = (1 << BLOCK_SIZES[size][0]);
            let h = (1 << BLOCK_SIZES[size][1]);
            if (size == i) {
              visitor(size, coordinates.set(c, r), subCoordinates.set(0, 0), bounds.set(c * s, r * s, w, h));
            }
          }
        }
      }
      // Visit blocks < 8x8.
      for (let c = 0; c < cols; c++) {
        for (let r = 0; r < rows; r++) {
          let size = this.aom.getMIProperty(MIProperty.GET_MI_BLOCK_SIZE, c, r);
          let w = (1 << BLOCK_SIZES[size][0]);
          let h = (1 << BLOCK_SIZES[size][1]);
          coordinates.set(c, r);
          switch (size) {
            case AOMAnalyzerBlockSize.BLOCK_4X4:
              visitor(size, coordinates, subCoordinates.set(0, 0), bounds.set(c * s,     r * s,     w, h));
              visitor(size, coordinates, subCoordinates.set(0, 1), bounds.set(c * s,     r * s + h, w, h));
              visitor(size, coordinates, subCoordinates.set(1, 0), bounds.set(c * s + w, r * s,     w, h));
              visitor(size, coordinates, subCoordinates.set(1, 1), bounds.set(c * s + w, r * s + h, w, h));
              break;
            case AOMAnalyzerBlockSize.BLOCK_8X4:
              visitor(size, coordinates, subCoordinates.set(0, 0), bounds.set(c * s,     r * s,     w, h));
              visitor(size, coordinates, subCoordinates.set(0, 1), bounds.set(c * s,     r * s + h, w, h));
              break;
            case AOMAnalyzerBlockSize.BLOCK_4X8:
              visitor(size, coordinates, subCoordinates.set(0, 0), bounds.set(c * s,     r * s,     w, h));
              visitor(size, coordinates, subCoordinates.set(1, 0), bounds.set(c * s + w, r * s,     w, h));
              break;
          }
        }
      }
    } else if (mode === BlockVisitorMode.TRANSFORM) {
      // Some code duplication here, to keep things simple.

      // Visit blocks >= 8x8.
      for (let i = 1; i < TRANSFORM_SIZES.length; i++) {
        let dc = 1 << (TRANSFORM_SIZES[i][0] - 3);
        let dr = 1 << (TRANSFORM_SIZES[i][1] - 3);
        for (let c = 0; c < cols; c += dc) {
          for (let r = 0; r < rows; r += dr) {
            let size = this.aom.getMIProperty(MIProperty.GET_MI_TRANSFORM_SIZE, c, r);
            let w = (1 << TRANSFORM_SIZES[size][0]);
            let h = (1 << TRANSFORM_SIZES[size][1]);
            if (size == i) {
              visitor(size, coordinates.set(c, r), subCoordinates.set(0, 0), bounds.set(c * s, r * s, w, h));
            }
          }
        }
      }
      // Visit blocks < 4x4.
      for (let c = 0; c < cols; c++) {
        for (let r = 0; r < rows; r++) {
          let size = this.aom.getMIProperty(MIProperty.GET_MI_TRANSFORM_SIZE, c, r);
          if (size != 0) {
            continue;
          }
          let w = (1 << TRANSFORM_SIZES[size][0]);
          let h = (1 << TRANSFORM_SIZES[size][1]);
          coordinates.set(c, r);
          switch (size) {
            case AOMAnalyzerBlockSize.BLOCK_4X4:
              visitor(size, coordinates, subCoordinates.set(0, 0), bounds.set(c * s,     r * s,     w, h));
              visitor(size, coordinates, subCoordinates.set(0, 1), bounds.set(c * s,     r * s + h, w, h));
              visitor(size, coordinates, subCoordinates.set(1, 0), bounds.set(c * s + w, r * s,     w, h));
              visitor(size, coordinates, subCoordinates.set(1, 1), bounds.set(c * s + w, r * s + h, w, h));
              break;
          }
        }
      }
    } else if (mode === BlockVisitorMode.SUPER_BLOCK) {
      for (let c = 0; c < cols; c += 8) {
        for (let r = 0; r < rows; r += 8) {
          let w = s * 8;
          let h = s * 8;
          visitor(0, coordinates.set(c, r), subCoordinates.set(0, 0), bounds.set(c * s, r * s, w, h));
        }
      }
    } else {
      console.error("Invalid Mode");
    }
  }

  drawMotionVectors(ctx: CanvasRenderingContext2D, src: Rectangle, dst: Rectangle) {
    let {cols, rows} = this.aom.getMIGridSize();
    let scale = dst.w / src.w;
    let scaledFrameSize = this.frameSize.clone().multiplyScalar(scale);
    ctx.save();
    ctx.globalAlpha = 1;
    let aColor = this.colorOptions.mv0Color.value;
    let bColor = this.colorOptions.mv1Color.value;
    ctx.fillStyle = aColor;
    ctx.lineWidth = scale / 2;

    ctx.translate(-src.x * scale, -src.y * scale);
    this.visitBlocks(BlockVisitorMode.BLOCK, (blockSize, coordinates, subCoordinates, bounds) => {
      bounds.multiplyScalar(scale);
      let o = bounds.getCenter();

      let a = this.getMotionVector(coordinates.c, coordinates.r, 0);
      let b = this.getMotionVector(coordinates.c, coordinates.r, 1);

      if (a.length() > 0) {
        ctx.globalAlpha = Math.min(0.3, a.length() / 128);
        ctx.fillStyle = aColor;
        ctx.fillRect(bounds.x, bounds.y, bounds.w, bounds.h);
      }

      if (b.length() > 0) {
        ctx.globalAlpha = Math.min(0.3, b.length() / 128);
        ctx.fillStyle = bColor;
        ctx.fillRect(bounds.x, bounds.y, bounds.w, bounds.h);
      }

      a.divideScalar(8 / scale);
      let va = o.clone().add(a);
      b.divideScalar(8 / scale);
      let vb = o.clone().add(b);

      // Draw small vectors with a ligher color.
      ctx.globalAlpha = Math.max(0.2, Math.min(a.length() + b.length(), 1));
      ctx.strokeStyle = aColor;
      drawVector(ctx, o, va);

      ctx.strokeStyle = bColor;
      drawVector(ctx, o, vb);

      // Draw Dot
      ctx.beginPath();
      ctx.arc(o.x, o.y, scale / 2, 0, Math.PI * 2, true);
      ctx.closePath();
      ctx.fill();


    });
    ctx.restore();
  }
}

function clamp(v, a, b) {
	if (v < a) {
		v = a;
	}
	if (v > b) {
		v = b;
	}
	return v;
}

function YUV2RGB(yValue, uValue, vValue) {
  let rTmp = yValue + (1.370705 * (vValue - 128));
  let gTmp = yValue - (0.698001 * (vValue - 128)) - (0.337633 * (uValue - 128));
  let bTmp = yValue + (1.732446 * (uValue - 128));
  let r = clamp(rTmp | 0, 0, 255) | 0;
  let g = clamp(gTmp | 0, 0, 255) | 0;
  let b = clamp(bTmp | 0, 0, 255) | 0;
  return (b << 16) | (g << 8) | (r << 0);
}

function getUrlParameters(): any {
  let url = window.location.search.substring(1);
  url = url.replace(/\/$/, ""); // Replace / at the end that gets inserted by browsers.
  let params = {};
  url.split('&').forEach(function (s) {
    let t = s.split('=');
    params[t[0]] = decodeURIComponent(t[1]);
  });
  return params;
};

angular
.module('AomInspectorApp', ['ngMaterial', 'color.picker'])
.config(['$mdIconProvider', function($mdIconProvider) {
    $mdIconProvider
      // .iconSet('social', 'img/icons/sets/social-icons.svg', 24)
      .defaultIconSet('img/icons/sets/core-icons.svg', 24);
}])
.filter('keyboardShortcut', function($window) {
  return function(str) {
    if (!str) return;
    let keys = str.split('-');
    let isOSX = /Mac OS X/.test($window.navigator.userAgent);
    let seperator = (!isOSX || keys.length > 2) ? '+' : '';
    let abbreviations = {
      M: isOSX ? '⌘' : 'Ctrl',
      A: isOSX ? 'Option' : 'Alt',
      S: 'Shift'
    };
    return keys.map(function(key, index) {
      let last = index == keys.length - 1;
      return last ? key : abbreviations[key];
    }).join(seperator);
  };
})
.controller('AppCtrl', ['$scope', '$interval', '$mdSidenav', AppCtrl])
.directive('selectOnClick', ['$window', function ($window) {
  return {
    restrict: 'A',
    link: function (scope, element, attrs) {
      element.on('click', function () {
        if (!$window.getSelection().toString()) {
          // Required for mobile Safari
          this.setSelectionRange(0, this.value.length)
        }
      });
    }
  };
}]);
;

window.Module = {
  noExitRuntime: true,
  preRun: [],
  postRun: [function() {
    // startVideo();
  }],
  memoryInitializerPrefixURL: "bin/",
  arguments: ['input.ivf', 'output.raw']
};