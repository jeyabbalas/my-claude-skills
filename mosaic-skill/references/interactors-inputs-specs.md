# Mosaic Interactors, Inputs & Specs Reference

## Interactors

Interactors add interactive behavior to plots. They listen to input events and update Selections.

### Interval Selection

```javascript
// 1D horizontal brush
vg.intervalX({ 
  as: selection,
  field: "date",        // Field to generate predicates for
  pixelSize: 1          // Snap resolution
})

// 1D vertical brush
vg.intervalY({ as: selection })

// 2D rectangular brush
vg.intervalXY({ 
  as: selection,
  xfield: "x",          // Optional custom field names
  yfield: "y"
})
```

### Toggle Selection

```javascript
// Point selection (click/shift-click)
vg.toggle({ 
  as: selection,
  channels: ["x", "y", "color"]  // Fields to include in predicates
})

// Convenience variants
vg.toggleX({ as: selection })
vg.toggleY({ as: selection })
vg.toggleColor({ as: selection })
```

### Nearest Value

```javascript
// Select nearest point on hover
vg.nearestX({ 
  as: selection,
  channels: ["z"]       // Include series in predicate
})
vg.nearestY({ as: selection })
```

### Highlight

```javascript
// Visual highlighting based on selection
vg.highlight({ 
  by: selection,
  opacity: 0.1,         // Non-selected opacity
  fill: "#ccc",         // Non-selected fill
  stroke: "#ccc"        // Non-selected stroke
})
```

### Pan & Zoom

```javascript
// Pan (drag) and zoom (scroll)
vg.panZoom({
  x: xSelection,        // Selection for x-axis domain
  y: ySelection,        // Selection for y-axis domain
  xfield: "x",
  yfield: "y"
})

// Apply selections to domain attributes
vg.plot(
  vg.dot(...),
  vg.panZoom({ x: $xs, y: $ys }),
  vg.xDomain($xs),      // Bind to selection
  vg.yDomain($ys)
)
```

## Input Widgets

All inputs can be bound to Params or Selections via the `as` option.

### Menu

```javascript
vg.menu({
  label: "Category",
  as: selection,                // Param or Selection
  field: "category",            // Field for predicates
  
  // Data-driven options
  from: "tableName",
  column: "category",
  filterBy: otherSelection,     // Filter menu options
  
  // OR static options
  options: [
    { label: "All", value: null },
    { label: "Option A", value: "a" },
    { label: "Option B", value: "b" }
  ],
  value: "a"                    // Initial value
})
```

### Search

```javascript
vg.search({
  label: "Search",
  as: selection,
  from: "tableName",
  column: "name",
  type: "contains",             // contains, prefix, suffix, regexp
  filterBy: otherSelection
})
```

### Slider

```javascript
vg.slider({
  label: "Value",
  as: param,
  min: 0,
  max: 100,
  step: 1,
  value: 50,                    // Initial value
  
  // OR data-driven range
  from: "tableName",
  column: "value"               // Auto min/max from data
})
```

### Table

```javascript
vg.table({
  from: "tableName",
  filterBy: selection,
  columns: ["name", "category", "value"],  // Column subset
  width: { 
    name: 200, 
    category: 100, 
    value: 80 
  },
  align: { value: "right" },
  format: { value: d => d.toFixed(2) },
  height: 400,                  // Table height in pixels
  rowBatch: 100                 // Rows per scroll batch
})
```

## Legends

### Color Legend

```javascript
// Inline legend
vg.plot(
  vg.dot(...),
  vg.colorLegend({ as: selection })  // Interactive filtering
)

// Standalone legend (references named plot)
vg.plot(
  vg.dot(...),
  vg.name("myPlot")
)
vg.colorLegend({ for: "myPlot", as: selection })
```

### Opacity Legend

```javascript
vg.opacityLegend({ as: selection })
vg.opacityLegend({ for: "namedPlot" })
```

## Plot Attributes

```javascript
vg.plot(
  // Marks and interactors...
  
  // Dimensions
  vg.width(680),
  vg.height(400),
  vg.margin(60),
  vg.marginTop(20),
  vg.marginRight(40),
  vg.marginBottom(40),
  vg.marginLeft(60),
  
  // Axis configuration
  vg.xAxis("bottom"),           // top, bottom, null
  vg.yAxis("left"),             // left, right, null
  vg.xLabel("X Axis Label"),
  vg.yLabel("Y Axis Label"),
  vg.xLabelAnchor("center"),    // center, left, right
  vg.xTicks(5),                 // Number of ticks
  vg.xTickFormat(".2f"),        // D3 format string
  vg.xGrid(true),
  
  // Scale domains
  vg.xDomain([0, 100]),         // Fixed domain
  vg.xDomain(vg.Fixed),         // Calculate once, then fix
  vg.xDomain(selection),        // Bind to selection
  vg.colorDomain(["a", "b"]),
  vg.colorRange(["red", "blue"]),
  
  // Scale types
  vg.xScale("log"),             // linear, log, sqrt, pow, symlog
  vg.colorScale("ordinal"),     // ordinal, categorical, quantize
  
  // Faceting
  vg.fx("category"),
  vg.fy("group"),
  
  // Other
  vg.name("plotName"),          // For standalone legends
  vg.style("overflow: visible")
)
```

## Declarative Specifications

### JSON Format

```json
{
  "meta": {
    "title": "Dashboard Title",
    "description": "Description text"
  },
  "config": {
    "preagg": true
  },
  "data": {
    "tableName": {
      "file": "path/to/data.parquet",
      "select": ["col1", "col2"],
      "where": "col1 > 0"
    },
    "csvData": {
      "file": "data.csv",
      "delimiter": ","
    }
  },
  "params": {
    "value": 50,
    "brush": { "select": "crossfilter" },
    "highlight": { "select": "intersect" }
  },
  "vconcat": [
    {
      "input": "slider",
      "label": "Value",
      "as": "$value",
      "min": 0,
      "max": 100
    },
    {
      "plot": [
        {
          "mark": "dot",
          "data": { "from": "tableName", "filterBy": "$brush" },
          "x": "xCol",
          "y": "yCol",
          "fill": "category"
        },
        { "select": "intervalXY", "as": "$brush" },
        { "select": "highlight", "by": "$brush" }
      ],
      "width": 500,
      "height": 400
    }
  ]
}
```

### YAML Format

```yaml
data:
  stocks: { file: data/stocks.parquet }

params:
  brush: { select: crossfilter }

vconcat:
  - plot:
      - mark: lineY
        data: { from: stocks, filterBy: $brush }
        x: Date
        y: Close
        stroke: Symbol
      - select: intervalX
        as: $brush
    width: 680
    height: 200
```

### SQL Expressions in Specs

```yaml
# Param reference in expression
y: { sql: "value + $offset" }

# Column reference via $$
fill: { sql: "$$columnParam" }  # Param as column name

# Aggregates
y: { sql: "count()" }
y: { sql: "avg(value)" }
```

### Parsing & Rendering

```javascript
import { parseSpec, astToDOM, astToESM } from "@uwdata/mosaic-spec";

// Parse JSON/YAML to AST
const ast = parseSpec(spec);

// Generate live DOM elements
const element = await astToDOM(ast);
document.body.appendChild(element);

// Generate JavaScript code
const jsCode = astToESM(ast);
```

## API Context

Create isolated Mosaic instances:

```javascript
import { createAPIContext, Coordinator, socketConnector } from "@uwdata/vgplot";

const api = createAPIContext({
  coordinator: new Coordinator(socketConnector("ws://localhost:8001/")),
  namedPlots: new Map(),
  extensions: { customMethod: () => {} }
});

// Use api instead of global vg
api.plot(
  api.lineY(api.from("data"), { x: "date", y: "value" })
);
```

## Jupyter Widget

```python
import mosaic_widget as mw
from mosaic_widget import MosaicWidget

# Create widget from spec
spec = {"plot": [{"mark": "dot", ...}]}
widget = MosaicWidget(spec)

# With pandas DataFrame
import pandas as pd
df = pd.read_csv("data.csv")
widget = MosaicWidget(spec, data={"tableName": df})
```
