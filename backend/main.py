import sqlite3
import uvicorn
import os
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from typing import Optional, List, Tuple

BACKEND_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(BACKEND_DIR)
DB_PATH = os.path.join(PROJECT_ROOT, "data", "genomics_subset.db")

app = FastAPI(title="Genomic Region Overlap Engine")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

"""Calculates the total unique base pairs covered by a list of intervals"""
def merge_and_sum(intersections: List[Tuple[int, int]]) -> int:
    if not intersections:
        return 0
    
    # Sort by start coordinate
    intersections.sort(key=lambda x: x[0])
    
    total_bp = 0
    curr_start, curr_end = intersections[0]
    
    for next_start, next_end in intersections[1:]:
        if next_start < curr_end:
            # Overlap found: extend the current boundary
            curr_end = max(curr_end, next_end)
        else:
            # No overlap: add the previous block and reset
            total_bp += (curr_end - curr_start)
            curr_start, curr_end = next_start, next_end
            
    # Add the final block
    total_bp += (curr_end - curr_start)
    return total_bp

def get_db():
    if not os.path.exists(DB_PATH):
        raise FileNotFoundError(f"Database not found at {DB_PATH}. Did you run indexer.py?")
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn

@app.get("/api/overlaps")
async def find_overlaps(
    chr: str, start: int, end: int, 
    tissue: Optional[str] = "All", 
    maxTracks: int = 50
):
    db = get_db()
    cursor = db.cursor()
    query_len = end - start
    
    query = """
    SELECT m.*, idx.start as t_start, idx.end as t_end
    FROM idx_intervals idx
    JOIN metadata m ON idx.id = m.id
    WHERE m.chrom = ? AND idx.start < ? AND idx.end > ?
    """
    params = [chr, end, start]

    if tissue and tissue != "All":
        query += " AND m.tissue = ?"
        params.append(tissue)

    cursor.execute(query, params)
    
    # Storage for track-wise intersections
    track_intersections = {} 
    track_metadata = {}

    for r in cursor:
        tid = r['track_id']
        
        # Calculate the actual intersection with the query window
        i_start = max(start, r['t_start'])
        i_end = min(end, r['t_end'])
        
        if tid not in track_intersections:
            track_intersections[tid] = []
            track_metadata[tid] = {
                "track_id": tid,
                "track_name": r['track_name'],
                "assay": r['assay'],
                "tissue": r['tissue'],
                "cell_type": r['cell_type'],
                "source": r['source']
            }
        
        track_intersections[tid].append((i_start, i_end))

    # Calculate final unique coverage for each track
    final_results = []
    for tid, intervals in track_intersections.items():
        # Merge overlaps to avoid double-counting
        unique_overlap_bp = merge_and_sum(intervals)
        
        res = track_metadata[tid]
        res["overlap_bp"] = unique_overlap_bp
        res["overlap_pct"] = round((unique_overlap_bp / query_len) * 100, 2)
        res["hit_count"] = len(intervals)
        res["target_interval_display"] = "; ".join([f"{int(s)}-{int(e)}" for s, e in intervals])
        
        final_results.append(res)

    sorted_results = sorted(final_results, key=lambda x: x['overlap_bp'], reverse=True)

    return {"results": sorted_results[:maxTracks]}

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=5000)