---
name: mosaic-visualization
description: Mosaic is a JavaScript framework for creating scalable, interactive data visualizations that can handle millions to billions of data points. Use this skill when building interactive dashboards, linked visualizations, cross-filtering interfaces, or any data visualization that needs to scale beyond traditional browser limits. Triggers include requests for (1) scalable visualizations with DuckDB, (2) linked/coordinated views, (3) cross-filtering dashboards, (4) grammar of graphics with vgplot, (5) declarative visualization specs in JSON/YAML, (6) density plots, heatmaps, contours, or hex binning, (7) interactive charts with brushing, panning, zooming, or (8) Jupyter notebook visualizations with DuckDB backend.
---

# Mosaic Visualization Library

Mosaic is an extensible framework for linking databases and interactive views. It leverages DuckDB for scalable data processing and enables real-time exploration of massive datasets.

## Core Architecture

```
┌─────────────┐     ┌─────────────────┐     ┌──────────────┐
│   Clients   │────▶│   Coordinator   │────▶│   DuckDB     │
│ (vgplot,    │◀────│ (query manager, │◀────│ (WASM/Server)│
│  inputs)    │     │  cache, optim.) │     └──────────────┘
└─────────────┘     └─────────────────┘
        │                   │
        └───────┬───────────┘
                ▼
    ┌───────────────────────┐
    │  Params & Selections  │
    │ (reactive variables)  │
    └───────────────────────┘
```

## Installation & Setup

```bash
npm install @uwdata/vgplot
```

### Browser (DuckDB-WASM)
```javascript
import { coordinator, wasmConnector } from "@uwdata/vgplot";
coordinator().databaseConnector(wasmConnector());
```

### With DuckDB Server
```javascript
import { coordinator, socketConnector } from "@uwdata/vgplot";
coordinator().databaseConnector(socketConnector("ws://localhost:8001/"));
```

## Quick Start

```javascript
import * as vg from "@uwdata/vgplot";

// Load data
await vg.coordinator().exec([
  vg.loadParquet("stocks", "data/stocks.parquet")
]);

// Create visualization
const chart = vg.plot(
  vg.lineY(vg.from("stocks"), { x: "Date", y: "Close" }),
  vg.width(680),
  vg.height(200)
);
document.body.appendChild(chart);
```

## Package Overview

| Package | Purpose |
|---------|---------|
| `@uwdata/vgplot` | Main package - re-exports core, sql, inputs, plot |
| `@uwdata/mosaic-core` | Coordinator, Params, Selections, Connectors |
| `@uwdata/mosaic-sql` | SQL query builder API |
| `@uwdata/mosaic-inputs` | Menu, Search, Slider, Table widgets |
| `@uwdata/mosaic-spec` | JSON/YAML declarative specifications |

For most applications, import only `@uwdata/vgplot`.

## Params & Selections

### Params (scalar reactive values)
```javascript
import { Param } from "@uwdata/mosaic-core";
const myParam = Param.value(5);
// Use in SQL: vg.sql`column * ${myParam}`
```

### Selections (filter predicates)
```javascript
import { Selection } from "@uwdata/mosaic-core";
Selection.single()       // Single active clause
Selection.union()        // Combine via OR
Selection.intersect()    // Combine via AND
Selection.crossfilter()  // Intersect with cross-filtering
```

## SQL Builder

See `references/sql-api.md` for complete Query API, aggregate functions, window functions, and data loading.

```javascript
import { Query, count, sum, avg } from "@uwdata/mosaic-sql";

Query.from("table")
  .select("category", { total: sum("value"), avg: avg("value") })
  .where(gt("value", 0))
  .groupby("category")
  .orderby("total");
```

## vgplot Components

### Plot Structure
```javascript
vg.plot(
  // Marks (chart layers)
  vg.barY(vg.from("data"), { x: "category", y: "value" }),
  // Interactors
  vg.intervalX({ as: selection }),
  vg.highlight({ by: selection }),
  // Legends
  vg.colorLegend({ as: selection }),
  // Attributes
  vg.width(500),
  vg.height(300)
)
```

### Marks Reference

See `references/marks-reference.md` for complete mark options.

| Mark | Use Case |
|------|----------|
| `barX`, `barY` | Categorical bar charts |
| `rectX`, `rectY` | Continuous rectangles, histograms |
| `lineX`, `lineY` | Time series, connected points |
| `areaX`, `areaY` | Filled line charts |
| `dot`, `circle`, `hexagon` | Scatter plots |
| `text`, `textX`, `textY` | Labels and annotations |
| `ruleX`, `ruleY` | Reference lines |
| `tickX`, `tickY` | Tick marks |
| `cell`, `cellX`, `cellY` | Heatmap cells |
| `contour` | Density contour lines |
| `heatmap`, `raster` | Density images |
| `hexbin` | Hexagonal binning |
| `density`, `densityX`, `densityY` | KDE plots |
| `denseLine` | Line density |
| `regressionY` | Linear regression |
| `geo`, `sphere`, `graticule` | Maps |
| `voronoi`, `delaunayMesh`, `hull` | Delaunay/Voronoi |

### Interactors

| Interactor | Purpose |
|------------|---------|
| `intervalX`, `intervalY` | 1D brush selection |
| `intervalXY` | 2D rectangular brush |
| `toggle`, `toggleX`, `toggleY`, `toggleColor` | Point selection |
| `nearestX`, `nearestY` | Hover nearest point |
| `panZoom` | Pan and zoom navigation |
| `highlight` | Highlight selected points |

### Inputs

```javascript
vg.menu({ label: "Category", as: selection, from: "data", column: "category" })
vg.search({ label: "Search", as: selection, from: "data", column: "name" })
vg.slider({ label: "Value", as: param, min: 0, max: 100, step: 1 })
vg.table({ from: "data", filterBy: selection, height: 300 })
```

### Layout

```javascript
vg.vconcat(plot1, plot2)           // Vertical stack
vg.hconcat(plot1, plot2)           // Horizontal stack
vg.vspace(10), vg.hspace(10)       // Spacing
```

## Declarative Specifications (JSON/YAML)

```yaml
data:
  mydata: { file: data/file.parquet }
params:
  brush: { select: crossfilter }
vconcat:
  - plot:
      - mark: barY
        data: { from: mydata, filterBy: $brush }
        x: category
        y: { sql: "count()" }
      - select: intervalX
        as: $brush
    width: 500
```

Parse and render:
```javascript
import { parseSpec, astToDOM } from "@uwdata/mosaic-spec";
const ast = parseSpec(jsonSpec);
const element = await astToDOM(ast);
```

## Performance Optimizations

Mosaic automatically applies:
- **M4 optimization**: Line/area marks reduce points to ~4 per pixel
- **Pre-aggregation**: Materialized views for filter groups
- **Caching**: Query result caching
- **Prefetching**: Anticipatory data loading on selection activation

Control optimizations:
```javascript
vg.from("data", { optimize: false })  // Disable M4
vg.coordinator().configure({ preagg: { enabled: false } })  // Disable pre-aggregation
```

## Data Loading

```javascript
await vg.coordinator().exec([
  vg.loadParquet("table", "path/to/file.parquet"),
  vg.loadCSV("table", "path/to/file.csv", { delimiter: "," }),
  vg.loadJSON("table", "path/to/file.json"),
]);
```

With filtering:
```javascript
vg.loadParquet("aapl", "stocks.parquet", { where: "Symbol = 'AAPL'" })
```

## Common Patterns

### Cross-filtering Dashboard
```javascript
const brush = vg.Selection.crossfilter();
vg.vconcat(
  vg.plot(
    vg.rectY(vg.from("data", { filterBy: brush }), 
      { x: vg.bin("delay"), y: vg.count() }),
    vg.intervalX({ as: brush })
  ),
  vg.plot(
    vg.rectY(vg.from("data", { filterBy: brush }), 
      { x: vg.bin("time"), y: vg.count() }),
    vg.intervalX({ as: brush })
  )
)
```

### Overview + Detail
```javascript
const brush = vg.Selection.intersect();
vg.vconcat(
  vg.plot(
    vg.areaY(vg.from("data"), { x: "date", y: "value" }),
    vg.intervalX({ as: brush })
  ),
  vg.plot(
    vg.areaY(vg.from("data", { filterBy: brush }), 
      { x: "date", y: "value" }),
    vg.yDomain(vg.Fixed)
  )
)
```

### Linked Scatter Plot Matrix
```javascript
const brush = vg.Selection.single();
// Create grid of plots sharing the same selection
vg.hconcat(
  vg.plot(
    vg.dot(vg.from("data"), { x: "a", y: "b" }),
    vg.intervalXY({ as: brush }),
    vg.highlight({ by: brush })
  ),
  vg.plot(
    vg.dot(vg.from("data"), { x: "c", y: "d" }),
    vg.intervalXY({ as: brush }),
    vg.highlight({ by: brush })
  )
)
```

## Key Encoding Channels

| Channel | Description |
|---------|-------------|
| `x`, `y` | Position |
| `x1`, `x2`, `y1`, `y2` | Range positions |
| `fill`, `stroke` | Colors |
| `fillOpacity`, `strokeOpacity` | Opacity |
| `r` | Radius (dots) |
| `z` | Series grouping |
| `fx`, `fy` | Faceting |

## Fixed Domains

Prevent scale domain "jumps" during filtering:
```javascript
vg.xDomain(vg.Fixed)  // Calculate once, then fix
vg.yDomain([0, 100])  // Explicit fixed domain
```

## Error Handling

```javascript
// Custom logger
coordinator().logger = {
  log: console.log,
  info: console.info,
  warn: console.warn,
  error: console.error
};

// Suppress logging
coordinator().logger = null;
```
