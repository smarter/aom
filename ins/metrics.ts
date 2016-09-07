class AccountingSymbol {
  constructor(public name: string, public bits: number, public samples: number, public x: number, public y: number) {
    // ...
  }
}

type AccountingSymbolMap = { [name: string]: AccountingSymbol };

class Accounting {
  symbols: AccountingSymbol [] = null;
  frameSymbols: AccountingSymbolMap = null;
  constructor(symbols: AccountingSymbol [] = []) {
    this.symbols = symbols;
  }
  createFrameSymbols() {
    if (this.frameSymbols) {
      return this.frameSymbols;
    }
    this.frameSymbols = Object.create(null);
    this.frameSymbols = Accounting.flatten(this.symbols);
    return this.frameSymbols;
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

  static getSortedSymbolNames(accountings: Accounting []): string [] {
    let set = {};
    accountings.forEach(accounting => {
      let frameSymbols = accounting.createFrameSymbols();
      for (let name in frameSymbols) {
        set[name] = undefined;
      }
    });
    let names = Object.keys(set);
    names.sort();
    return names;
  }

  static makeSymbolsDataTable(accountings: Accounting [], count) {
    var data = new google.visualization.DataTable();
    data.addColumn('string', "Frame");
    let symbolNames = Accounting.getSortedSymbolNames(accountings);
    symbolNames.forEach(symbolName => {
      data.addColumn('number', symbolName);
    });
    let offset = Math.max(0, accountings.length - count);
    accountings = accountings.slice(offset)
    for (let i = 0; i < accountings.length; i++) {
      let frameSymbols = accountings[i].createFrameSymbols();
      let row = symbolNames.map(symbolName => {
        if (frameSymbols[symbolName]) {
          return frameSymbols[symbolName].bits / 8;
        }
        return 0;
      }) as any [];
      row.unshift("Frame " + (i + offset));
      data.addRows([row]);
    }

    // var csv = google.visualization.dataTableToCsv(data);
    // console.log(csv);
    return data;
  }

  static makeTotalBitsDataTable(accountings: Accounting [], count) {
    var data = new google.visualization.DataTable();

    data.addColumn('string', "Frame");
    data.addColumn('number', "Bits");
    let offset = Math.max(0, accountings.length - count);
    accountings = accountings.slice(offset)
    for (let i = 0; i < accountings.length; i++) {
      let frameSymbols = accountings[i].createFrameSymbols();
      let total = 0;
      for (let symbolName in frameSymbols) {
        total += frameSymbols[symbolName].bits
      }
      data.addRows([["Frame " + (i + offset), total / 8]]);
    }
    return data;
  }

  static makeFrameErrorDataTable(frameErrors: ErrorMetrics [], count) {
    var data = new google.visualization.DataTable();
    data.addColumn('string', "Frame");
    data.addColumn('number', "Error");
    let offset = Math.max(0, frameErrors.length - count);
    frameErrors = frameErrors.slice(offset)
    for (let i = 0; i < frameErrors.length; i++) {
      data.addRows([["Frame " + (i + offset), frameErrors[i].mse]]);
    }
    return data;
  }
}

class FrameModeMetrics {
  counts: Uint32Array;
  constructor() {
    this.counts = new Uint32Array(AnalyzerPredictionMode.NEWMV + 1);
  }
  static makeFrameModeDataTable(metrics: FrameModeMetrics [], count: number) {
    var data = new google.visualization.DataTable();
    data.addColumn('string', "Frame");
    for (let i = 0; i <= AnalyzerPredictionMode.NEWMV; i++) {
      data.addColumn('number', AnalyzerPredictionMode[i]);
    }
    let offset = Math.max(0, metrics.length - count);
    metrics = metrics.slice(offset)
    for (let i = 0; i < metrics.length; i++) {
      let row: any [] = ["Frame " + (i + offset)];
      for (let j = 0; j <= AnalyzerPredictionMode.NEWMV; j++) {
        row.push(metrics[i].counts[j]);
      }
      data.addRows([row]);
    }
    return data;
  }
}

class FrameBlockSizeMetrics {
  counts: Uint32Array;
  constructor() {
    this.counts = new Uint32Array(AnalyzerBlockSize.BLOCK_64X64 + 1);
  }
  static makeFrameModeDataTable(metrics: FrameBlockSizeMetrics [], count: number) {
    var data = new google.visualization.DataTable();
    data.addColumn('string', "Frame");
    for (let i = 0; i <= AnalyzerBlockSize.BLOCK_64X64; i++) {
      data.addColumn('number', AnalyzerBlockSize[i]);
    }
    let offset = Math.max(0, metrics.length - count);
    metrics = metrics.slice(offset)
    for (let i = 0; i < metrics.length; i++) {
      let row: any [] = ["Frame " + (i + offset)];
      for (let j = 0; j <= AnalyzerBlockSize.BLOCK_64X64; j++) {
        row.push(metrics[i].counts[j]);
      }
      data.addRows([row]);
    }
    return data;
  }
}