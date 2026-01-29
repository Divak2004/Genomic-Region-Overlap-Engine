# Genomic Region Overlap Engine
An engine that tries to identify which genomic tracks overlap a user-specified genomic region.

## Overlap Logic

The engine employs a multi-stage process to ensure both speed and biological accuracy:

1.  **Spatial Filtering:** Instead of a linear scan, the engine utilizes a SQLite R-Tree index. It treats genomic intervals as 1D spatial objects, allowing the database to instantly discard non-overlapping tracks. The query logic is defined as `Target_Start < Query_End` **AND** `Target_End > Query_Start`.
2.  **Coordinate Clipping:** Once candidate intervals are found, the engine clips them to the query window boundaries using `max(Query_Start, Target_Start)` and `min(Query_End, Target_End)`.
3.  **Union of Intervals:** To handle tracks with multiple disjoint or overlapping peaks within the same query window, the engine applies a merging algorithm, preventing double-counting when calculating how much overlap there is. 

---

## Design Assumptions

* **Sampling:** We assume that taking a random sample of 2,500 tracks, with 1,000 intervals per track will provide enough representation of the genome data for a mini-project.
* **Source availability:** We assume that the URLs for the metadata and the format are up so we can update our database whenever needed.
* **User assumption:** We assume that the user will not try to overrun our server by making invalid requests and/or flooding the server with an egregious number of requests.
* **Data quality:** We assume that most of the data is "nice" in that most fields that we extract from are filled and not malformed.

---

## Scaling to Millions of Regions

The current architecture is designed to handle high-volume datasets with the following performance characteristics:

* **Search Complexity (average $O(\log N)$):** By utilizing an R-Tree spatial index, query time scales logarithmically with the number of genomic intervals. This ensures that the engine remains responsive even as the database grows from thousands to millions of intervals.
* **Memory Efficiency:** The engine employs a "Streaming & Tally" approach. It does not load the entire database into memory; instead, it uses SQLite's row-iterator to process hits. This allows the system to run on low-memory environments (like a standard laptop) while querying massive datasets.
* **Stateless Concurrency:** The FastAPI backend is stateless, allowing it to scale horizontally behind a load balancer or vertically across multiple CPU cores using Uvicorn workers.

---

## Future Improvements (With More Time)
* **Deployment:** Deploy the project and server, allowing other users to be able to make queries.
* **API Authentication:** Implement JWT-based auth and rate-limiting to prepare the backend for public production use.
* **Database usage:** Research to find a database more suitable for storing and querying interval-based queries and supporting more writes (as SQLite does not support write-heavy purposes). There may be more writes as new tracks are made, old tracks need to be fixed etc. 
* **Track coverage:** Cover all the tracks with a potentially superior database, or implement another method in which all tracks can be accessed at the cost of speed.
* **Better filtering:** Add more tissue options (perhaps even an option that lets you filter the tissue once all tracks have been generated for a query), as well as other filtering options like data source.

---

## Installation & Usage

### 1. Setup
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 2. Unzip (or build) database
Unzip the database named `genomics_subset.db.zip`. Alternatively, run the following (takes around ~10-15 minutes):
```
python3 data/indexer.py
```

### 3. Start backend
```
python3 backend/main.py
```

### 4. Launch frontend
In a different terminal, run the following:
```
open frontend/index.html
```