import gzip
import requests
import pandas as pd
import json
import os
from io import StringIO
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from concurrent.futures import ThreadPoolExecutor

app = FastAPI(title="Genomic Overlap Engine")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Global State ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "..", "data")
METADATA_URL = "https://tf.lisanwanglab.org/FILER2/metadata/tracks.metadata.tsv"
metadata_cache = None
BOUNDARY_INDEX = {}

def load_boundary_index():
    global BOUNDARY_INDEX
    index_path = os.path.join(DATA_DIR, "track_boundaries.json")
    if os.path.exists(index_path):
        with open(index_path, 'r') as f:
            BOUNDARY_INDEX = json.load(f)
        print(f"Loaded boundary index for {len(BOUNDARY_INDEX)} tracks.")

def load_metadata():
    global metadata_cache
    if metadata_cache is None:
        r = requests.get(METADATA_URL)
        metadata_cache = pd.read_csv(StringIO(r.text), sep='\t').dropna(subset=['processed_file_download_url'])
    return metadata_cache

def check_overlap(q_start, q_end, t_start, t_end):
    """Simple range intersection check."""
    return q_start < t_end and t_start < q_end

def fetch_and_parse_stream(track_info, query_chr, q_start, q_end):
    """Streams the .gz file and manually checks for overlaps."""
    url = track_info['processed_file_download_url']
    norm_query_chr = query_chr.lower().replace('chr', '')
    
    try:
        # Stream the file with a timeout
        response = requests.get(url, stream=True, timeout=5)
        if response.status_code != 200:
            return None

        overlapping_regions = []
        # Use gzip.open on the raw response stream
        with gzip.open(response.raw, 'rt') as f:
            for i, line in enumerate(f):
                # Performance safeguard: don't scan more than 20,000 lines per file
                if i > 20000: break 
                if line.startswith(('#', 'track', 'browser')): continue
                
                fields = line.strip().split('\t')
                if len(fields) < 3: continue
                
                # Normalize chromosome for comparison
                t_chrom = fields[0].lower().replace('chr', '')
                if t_chrom != norm_query_chr: continue
                
                try:
                    t_start, t_end = int(fields[1]), int(fields[2])
                    
                    # Optimization: if files are sorted and we passed the query, stop
                    if t_start > q_end: break

                    if check_overlap(q_start, q_end, t_start, t_end):
                        overlapping_regions.append({
                            "chrom": fields[0],
                            "start": t_start,
                            "end": t_end,
                            "extra": fields[3:6]
                        })
                        if len(overlapping_regions) >= 10: break
                except ValueError:
                    continue

        if not overlapping_regions:
            return None

        return {
            "track_name": track_info['track_name'],
            "assay": track_info.get('assay', 'N/A'),
            "overlap_count": len(overlapping_regions),
            "results": overlapping_regions
        }
    except Exception:
        return None

@app.get("/api/overlaps")
async def find_overlaps(chr: str, start: int, end: int, maxTracks: int = 20):
    meta = load_metadata()
    norm_chr = chr.lower().replace('chr', '')
    
    # 1. PRUNING: Use track_boundaries.json to identify candidates
    candidate_tracks = []
    for _, track in meta.iterrows():
        tid = track['identifier']
        if tid in BOUNDARY_INDEX:
            bounds = BOUNDARY_INDEX[tid]
            if norm_chr in bounds:
                t_min, t_max = bounds[norm_chr]
                # If query is completely outside the range, skip network call
                if end < t_min or start > t_max:
                    continue
            else:
                continue # Chromosome not in this track
        
        candidate_tracks.append(track.to_dict())
        if len(candidate_tracks) >= maxTracks:
            break

    # 2. SCANNING: Threaded streaming calls
    final_results = []
    # Lower worker count (10) for streaming to avoid overwhelming the server
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = [executor.submit(fetch_and_parse_stream, t, chr, start, end) for t in candidate_tracks]
        for f in futures:
            res = f.result()
            if res:
                final_results.append(res)

    return {
        "query": {"chr": chr, "start": start, "end": end},
        "tracks_indexed_as_candidates": len(candidate_tracks),
        "total_tracks_found": len(final_results),
        "results": final_results
    }

@app.on_event("startup")
async def startup_event():
    load_metadata()
    load_boundary_index()

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=5000)