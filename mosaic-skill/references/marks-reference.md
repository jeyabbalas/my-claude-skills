# Mosaic Marks Reference

Complete reference for vgplot marks with options and examples.

## Mark Syntax

```javascript
vg.markName(vg.from("table", options), markOptions)
```

### Data Source Options

```javascript
vg.from("tableName", {
  filterBy: selection,    // Selection for filtering
  optimize: true          // Enable M4 optimization (line/area)
})
```

### Common Mark Options

All marks support these options:

| Option | Description |
|--------|-------------|
| `x`, `y` | Position channels |
| `fill` | Fill color |
| `stroke` | Stroke color |
| `fillOpacity` | Fill opacity (0-1) |
| `strokeOpacity` | Stroke opacity (0-1) |
| `strokeWidth` | Stroke width in pixels |
| `opacity` | Overall opacity |
| `clip` | Clip to plot frame |
| `title` | Tooltip text |
| `ariaLabel` | Accessibility label |
| `select` | Filter value (for annotation marks) |

## Basic Marks

### Bar (`barX`, `barY`)
Categorical bar charts with ordinal dimension.

```javascript
vg.barY(vg.from("data"), {
  x: "category",        // Ordinal x-axis
  y: "value",           // Quantitative y-axis
  fill: "steelblue",
  inset: 1,             // Gap between bars
  rx: 2                 // Corner radius
})
```

### Rect (`rectX`, `rectY`)
Rectangles with continuous dimensions.

```javascript
vg.rectY(vg.from("data"), {
  x: vg.bin("value"),   // Binned x-axis
  y: vg.count(),        // Count aggregation
  fill: "steelblue",
  inset: 0.5
})
```

### Dot (`dot`, `dotX`, `dotY`, `circle`, `hexagon`)
Scatter plots and point marks.

```javascript
vg.dot(vg.from("data"), {
  x: "xVar",
  y: "yVar",
  r: 3,                 // Radius
  fill: "category",     // Color by category
  stroke: "#fff",
  strokeWidth: 0.5,
  symbol: "circle"      // circle, cross, diamond, square, star, triangle, wye
})
```

### Text (`text`, `textX`, `textY`)
Labels and annotations.

```javascript
vg.text(vg.from("data"), {
  x: "xVar",
  y: "yVar",
  text: "label",
  fontSize: 12,
  fontWeight: "bold",
  textAnchor: "middle", // start, middle, end
  dy: -8,               // Vertical offset
  dx: 0                 // Horizontal offset
})
```

### Rule (`ruleX`, `ruleY`)
Reference lines spanning the plot.

```javascript
vg.ruleY([0])           // Horizontal line at y=0
vg.ruleX([threshold], { stroke: "red", strokeDasharray: "4,2" })
```

### Tick (`tickX`, `tickY`)
Short tick marks.

```javascript
vg.tickX(vg.from("data"), {
  x: "value",
  y: "category",        // Ordinal position
  stroke: "steelblue"
})
```

### Cell (`cell`, `cellX`, `cellY`)
Heatmap cells in ordinal dimensions.

```javascript
vg.cell(vg.from("data"), {
  x: "xCategory",
  y: "yCategory",
  fill: "value",        // Color by value
  inset: 0.5
})
```

## Connected Marks

### Line (`line`, `lineX`, `lineY`)
Connected line charts with M4 optimization.

```javascript
vg.lineY(vg.from("data"), {
  x: "date",
  y: "value",
  z: "series",          // Group by series
  stroke: "series",     // Color by series
  strokeWidth: 1.5,
  curve: "linear"       // linear, step, natural, monotone-x, etc.
})
```

### Area (`area`, `areaX`, `areaY`)
Filled area charts with M4 optimization.

```javascript
vg.areaY(vg.from("data"), {
  x: "date",
  y: "value",
  y1: 0,                // Baseline
  fill: "steelblue",
  fillOpacity: 0.5,
  curve: "monotone-x"
})
```

## Density Marks

### Density 1D (`densityX`, `densityY`)
Kernel density estimation.

```javascript
vg.densityY(vg.from("data"), {
  x: "value",
  fill: "steelblue",
  bandwidth: 20,        // Smoothing bandwidth
  bins: 1024,           // Number of bins
  type: "areaX"         // Mark type: area, line, dot
})
```

### Density 2D (`density`)
2D point density.

```javascript
vg.density(vg.from("data"), {
  x: "xVar",
  y: "yVar",
  fill: "density",      // Map density to color
  bandwidth: 20,
  pixelSize: 2
})
```

### Contour
Density contour lines.

```javascript
vg.contour(vg.from("data"), {
  x: "xVar",
  y: "yVar",
  fill: "density",
  stroke: "density",
  thresholds: 10,       // Number of contour levels
  bandwidth: 20
})
```

### Heatmap
Smoothed density raster.

```javascript
vg.heatmap(vg.from("data"), {
  x: "xVar",
  y: "yVar",
  fill: "density",
  bandwidth: 20,        // Default smoothing
  pixelSize: 2
})
```

### Raster
Custom raster visualization.

```javascript
vg.raster(vg.from("data"), {
  x: "xVar",
  y: "yVar",
  fill: vg.avg("value"),  // Aggregate per pixel
  bandwidth: 0,            // No smoothing
  interpolate: "nearest",  // none, linear, nearest, barycentric, random-walk
  pixelSize: 1
})
```

### Hexbin
Hexagonal binning.

```javascript
vg.hexbin(vg.from("data"), {
  x: "xVar",
  y: "yVar",
  fill: vg.count(),
  r: vg.count(),          // Size by count
  binWidth: 20,           // Hex size in pixels
  type: "hexagon"         // Can use other marks
})
```

### Dense Line
Line density estimation.

```javascript
vg.denseLine(vg.from("data"), {
  x: "xVar",
  y: "yVar",
  z: "series",
  bandwidth: 0,
  pixelSize: 1,
  normalize: true         // Arc-length normalization
})
```

## Regression

```javascript
vg.regressionY(vg.from("data"), {
  x: "predictor",
  y: "response",
  stroke: "steelblue",
  ci: 0.95,              // Confidence interval
  precision: 4           // Pixels between samples
})
```

## Geographic Marks

### Geo
Geographic features.

```javascript
vg.geo(vg.from("geodata"), {
  geometry: "geom",       // Column with GeoJSON
  fill: "value",
  stroke: "#fff",
  strokeWidth: 0.5
})
```

### Sphere and Graticule
```javascript
vg.sphere({ fill: "aliceblue" })
vg.graticule({ stroke: "#ccc", strokeOpacity: 0.5 })
```

## Delaunay/Voronoi

```javascript
vg.voronoi(vg.from("data"), { x: "x", y: "y", fill: "category" })
vg.voronoiMesh(vg.from("data"), { x: "x", y: "y", stroke: "#ccc" })
vg.delaunayLink(vg.from("data"), { x: "x", y: "y" })
vg.delaunayMesh(vg.from("data"), { x: "x", y: "y" })
vg.hull(vg.from("data"), { x: "x", y: "y", fill: "category" })
```

## Other Marks

### Frame
Plot frame rectangle.
```javascript
vg.frame({ stroke: "#ccc", fill: "none" })
```

### Arrow
Arrows between points.
```javascript
vg.arrow(vg.from("data"), {
  x1: "startX", y1: "startY",
  x2: "endX", y2: "endY",
  stroke: "steelblue"
})
```

### Link
Lines between points.
```javascript
vg.link(vg.from("data"), {
  x1: "x1", y1: "y1",
  x2: "x2", y2: "y2"
})
```

### Image
Images at positions.
```javascript
vg.image(vg.from("data"), {
  x: "x", y: "y",
  src: "imageUrl",
  width: 50, height: 50
})
```

### Vector/Spike
Vector field visualization.
```javascript
vg.vector(vg.from("data"), {
  x: "x", y: "y",
  length: "magnitude",
  rotate: "angle"
})
```

## Transforms

### Binning
```javascript
x: vg.bin("value")              // Auto-binned
x: vg.bin("value", { steps: 20 })
```

### Aggregates in Marks
```javascript
y: vg.count()
y: vg.sum("value")
y: vg.avg("value")
fill: vg.max("value")
```

### Stack
```javascript
vg.barY(vg.from("data"), {
  x: "category",
  y: vg.sum("value"),
  fill: "series",
  sort: { color: "width" }
})
```

## Direct Data (No Query)

For annotations without database queries:
```javascript
vg.ruleY([0, 50, 100])          // Array of values
vg.text([{ x: 10, y: 20, label: "Annotation" }], { 
  x: "x", y: "y", text: "label" 
})
```
