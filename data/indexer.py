import pandas as pd
import sqlite3
import requests
import gzip
import os
import random
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

DATA_DIR = os.path.dirname(os.path.abspath(__file__))
DB_NAME = os.path.join(DATA_DIR, "genomics_subset.db")
SUBSET_COUNT = 2500
WORKERS = os.cpu_count() or 4

def get_indices(schema):
    """Maps FILER_BED_schema to the correct (chrom, start, end) column indices."""
    # Interaction tracks (Hi-C) store the interaction span in columns 5, 6, 7
    if schema == "bed4+19 interact":
        return 5, 6, 7
    # Standard BED format is 0, 1, 2
    return 0, 1, 2

def process_track(track_data):
    url = track_data['processed_file_download_url']
    schema = track_data.get('FILER_BED_schema', 'bed3')
    c_idx, s_idx, e_idx = get_indices(schema)
    
    # Target number of random samples per track
    RESERVOIR_SIZE = 1000
    reservoir = []
    
    meta = {
        "tid": track_data['identifier'],
        "name": str(track_data.get('track_name', 'Unnamed')),
        "tissue": str(track_data.get('tissue_category', 'Unknown')),
        "cell": str(track_data.get('cell_type', 'N/A')),
        "assay": str(track_data.get('assay', 'Unknown')),
        "source": str(track_data.get('data_source', 'Unknown'))
    }

    try:
        r = requests.get(url, stream=True, timeout=10)
        # We must decompress the stream to read line by line
        with gzip.open(r.raw, 'rt') as f:
            for i, line in enumerate(f):
                if line.startswith(('#', 'track')): continue
                
                # Reservoir Sampling Algorithm
                if len(reservoir) < RESERVOIR_SIZE:
                    reservoir.append(line)
                else:
                    # Randomly replace existing lines with decreasing probability
                    m = random.randint(0, i)
                    if m < RESERVOIR_SIZE:
                        reservoir[m] = line

        # Final Parsing of the randomly selected lines
        intervals = []
        for line in reservoir:
            fields = line.strip().split('\t')
            if len(fields) <= max(c_idx, s_idx, e_idx): continue
            
            try:
                intervals.append({
                    **meta,
                    "chrom": fields[c_idx],
                    "start": int(fields[s_idx]),
                    "end": int(fields[e_idx])
                })
            except ValueError:
                continue
        return intervals
    except Exception:
        return []

def build_index():
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)

    if os.path.exists(DB_NAME):
        os.remove(DB_NAME)

    conn = sqlite3.connect(DB_NAME)
    c = conn.cursor()
    
    # R-Tree virtual table for O(log N) spatial searching
    c.execute("CREATE VIRTUAL TABLE idx_intervals USING rtree(id, start, end)")
    # Relational table for biological metadata
    c.execute("""CREATE TABLE metadata (
        id INTEGER PRIMARY KEY, track_id TEXT, track_name TEXT, 
        chrom TEXT, tissue TEXT, cell_type TEXT, assay TEXT, source TEXT
    )""")
    conn.commit()

    print("Fetching and Merging FILER2 Metadata...")
    df_meta = pd.read_csv("https://tf.lisanwanglab.org/FILER2/metadata/tracks.metadata.tsv", sep='\t')
    df_formats = pd.read_csv("https://tf.lisanwanglab.org/FILER2/metadata/track.formats.4col.tsv", sep='\t')
    
    # Merge to ensure we know the schema/format of every track
    df = df_meta.merge(df_formats, left_on='file_format', right_on='FILER_BED_format', how='left')
    tracks = df.dropna(subset=['processed_file_download_url']).sample(n=SUBSET_COUNT).to_dict('records')

    with ProcessPoolExecutor(max_workers=WORKERS) as executor:
        futures = {executor.submit(process_track, t): t for t in tracks}
        with tqdm(total=len(tracks), desc="Sampling Tracks") as pbar:
            for future in as_completed(futures):
                data = future.result()
                if data:
                    for row in data:
                        c.execute("""INSERT INTO metadata 
                            (track_id, track_name, chrom, tissue, cell_type, assay, source) 
                            VALUES (?,?,?,?,?,?,?)""",
                            (row['tid'], row['name'], row['chrom'], row['tissue'], 
                             row['cell'], row['assay'], row['source']))
                        c.execute("INSERT INTO idx_intervals VALUES (?,?,?)", (c.lastrowid, row['start'], row['end']))
                    conn.commit()
                pbar.update(1)
    
    conn.close()
    print(f"✅ Database {DB_NAME} built with random genome-wide sampling.")

if __name__ == "__main__":
    build_index()