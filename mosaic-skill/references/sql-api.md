# Mosaic SQL API Reference

Complete reference for `@uwdata/mosaic-sql` query construction.

## Query Builder

```javascript
import { Query } from "@uwdata/mosaic-sql";

Query
  .with({ cte: subquery })           // Common table expressions
  .select("col1", { alias: expr })   // Columns and expressions
  .distinct(true)                     // Distinct values only
  .from("table")                      // Source tables/subqueries
  .sample(1000)                       // Row sampling
  .where(predicate)                   // Filter criteria
  .groupby("col1", "col2")           // Group by columns
  .having(predicate)                  // Post-aggregation filter
  .window({ name: windowDef })        // Named windows
  .qualify(predicate)                 // Post-window filter
  .orderby(desc("col1"))             // Sort order
  .limit(100)                         // Max rows
  .offset(50)                         // Skip rows
```

## Operators

```javascript
import { 
  eq, neq, lt, gt, lte, gte,        // Comparisons
  and, or, not,                      // Logical
  isNull, isNotNull,                 // Null checks
  isBetween, isNotBetween,           // Range
  add, sub, mul, div, mod,           // Arithmetic
  literal, column, sql               // Builders
} from "@uwdata/mosaic-sql";

// Examples
where(and(gt("value", 0), lt("value", 100)))
where(isBetween("date", "2020-01-01", "2020-12-31"))
select({ calc: add("a", mul("b", 2)) })
```

## Aggregate Functions

```javascript
import {
  count, sum, avg, min, max,
  median, mode,
  variance, stddev, skewness, kurtosis,
  quantile,                          // quantile("col", 0.5) = median
  first, last, argmin, argmax,
  corr, covariance,                  // Two-column aggregates
  entropy,
  arrayAgg, stringAgg
} from "@uwdata/mosaic-sql";

// With modifiers
count().distinct()                    // COUNT(DISTINCT *)
sum("value").filter(gt("x", 0))      // SUM(value) FILTER (WHERE x > 0)
```

## Window Functions

```javascript
import {
  row_number, rank, dense_rank, ntile,
  lag, lead, first_value, last_value,
  cume_dist, percent_rank
} from "@uwdata/mosaic-sql";

// Make any aggregate a window function
avg("value")
  .partitionby("category")
  .orderby("date")
  .rows([-3, 3])                      // 7-day moving average

// Frame specifications
.rows([null, 0])      // Unbounded preceding to current
.rows([0, null])      // Current to unbounded following
.range([-7, 0])       // Value-based range
```

## Date Functions

```javascript
import {
  year, month, day, hour, minute, second,
  dayofweek, dayofyear, week,
  date_trunc, date_part, epoch_ms,
  dateMonth, dateMonthDay, dateDay
} from "@uwdata/mosaic-sql";

// Examples
select({ yr: year("date"), mo: month("date") })
date_trunc("month", "date")           // Truncate to month
```

## SQL Template Literal

```javascript
import { sql, column } from "@uwdata/mosaic-sql";

// Raw SQL expressions
const expr = sql`log(${column("value")} + 1)`;

// With params (reactive)
import { Param } from "@uwdata/mosaic-core";
const threshold = Param.value(100);
const expr = sql`${column("value")} > ${threshold}`;
```

## Column References

```javascript
import { column, all } from "@uwdata/mosaic-sql";

column("name")              // Reference column "name"
column("table.name")        // Qualified reference
all()                       // SELECT *
```

## Data Loading Functions

```javascript
import { 
  loadCSV, loadJSON, loadParquet, loadObjects 
} from "@uwdata/mosaic-sql";

// Load from files
loadCSV("tableName", "path/to/file.csv", {
  delimiter: ",",
  header: true,
  columns: ["a", "b", "c"],
  select: ["a", "b"],
  where: "a > 0"
});

loadParquet("tableName", "path/to/file.parquet", {
  select: ["col1", "col2"],
  where: "date > '2020-01-01'"
});

loadJSON("tableName", "path/to/file.json");

// Load from JavaScript objects
loadObjects("tableName", [
  { foo: 1, bar: 2 },
  { foo: 3, bar: 4 }
]);
```

## Query Analysis

```javascript
import { collectColumns, collectParams } from "@uwdata/mosaic-sql";

const query = Query.from("table").select("a", { b: sum("c") });

collectColumns(query);  // Extract column references
collectParams(query);   // Extract param references
```

## Cast & Type Conversion

```javascript
import { cast } from "@uwdata/mosaic-sql";

cast("column", "INTEGER")
cast("column", "VARCHAR")
cast("column", "TIMESTAMP")
```

## Case Expressions

```javascript
import { caseWhen } from "@uwdata/mosaic-sql";

caseWhen([
  [gt("value", 100), literal("high")],
  [gt("value", 50), literal("medium")]
], literal("low"))  // ELSE clause
```

## Set Operations

```javascript
import { union, unionAll, intersect, except } from "@uwdata/mosaic-sql";

union(query1, query2)
unionAll(query1, query2)
intersect(query1, query2)
except(query1, query2)
```

## Subqueries

```javascript
// Subquery in FROM
Query.from({ sub: Query.from("table").where(condition) })
  .select("*")

// Subquery in WHERE
Query.from("t1")
  .where(isIn("id", Query.from("t2").select("id")))
```
