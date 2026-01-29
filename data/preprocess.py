import pandas as pd
import requests
import gzip
import json
import os
from io import BytesIO, StringIO
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

# Configuration
METADATA_URL = "https://tf.lisanwanglab.org/FILER2/metadata/tracks.metadata.tsv"
FORMATS_URL = "https://tf.lisanwanglab.org/FILER2/metadata/track.formats.4col.tsv"
OUTPUT_PATH = "data/track_boundaries.json"

def load_schema_map():
    print("Fetching format schemas...")
    r = requests.get(FORMATS_URL)
    df = pd.read_csv(StringIO(r.text), sep='\t')
    mapping = {}
    for _, row in df.iterrows():
        cols = row['FILER_BED_schema'].split(';')
        try:
            mapping[row['FILER_BED_type']] = {
                'chrom': next(i for i, c in enumerate(cols) if 'chrom' in c.lower() and 'start' not in c.lower() and 'end' not in c.lower()),
                'start': next(i for i, c in enumerate(cols) if 'start' in c.lower()),
                'end': next(i for i, c in enumerate(cols) if 'end' in c.lower())
            }
        except:
            mapping[row['FILER_BED_type']] = {'chrom': 0, 'start': 1, 'end': 2}
    return mapping

def get_track_bounds(args):
    row, schema_map = args
    url = row['processed_file_download_url']
    if pd.isna(url): return None
    
    col = schema_map.get(row.get('FILER_BED_type'), {'chrom': 0, 'start': 1, 'end': 2})
    
    try:
        # Range Request to grab just the start of the file
        headers = {"Range": "bytes=0-153600"} 
        response = requests.get(url, headers=headers, timeout=10)
        
        with gzip.open(BytesIO(response.content), 'rt') as f:
            for line in f:
                if line.startswith(('#', 'track', 'browser', '<!DOCTYPE')):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) <= max(col.values()): continue
                
                chrom = fields[col['chrom']].replace('chr', '')
                start = int(fields[col['start']])
                
                # Boundary Strategy: Save the first chromosome and starting position found.
                return (row['identifier'], {chrom: [start, start + 10_000_000]})
    except:
        return None

def run_preprocess():
    os.makedirs('data', exist_ok=True)
    
    # 1. Load existing progress if available (Checkpoint System)
    existing_data = {}
    if os.path.exists(OUTPUT_PATH):
        with open(OUTPUT_PATH, 'r') as f:
            existing_data = json.load(f)
        print(f"Found existing index with {len(existing_data)} tracks. Resuming...")

    schema_map = load_schema_map()
    print("Loading track metadata...")
    df = pd.read_csv(METADATA_URL, sep='\t')
    
    # Filter out tracks we already indexed
    tasks_df = df[~df['identifier'].isin(existing_data.keys())]
    
    if tasks_df.empty:
        print("All tracks already indexed!")
        return

    print(f"Indexing {len(tasks_df)} remaining tracks (Total: {len(df)}).")
    
    boundary_index = existing_data
    
    # 2. Process with Progress Bar
    with ThreadPoolExecutor(max_workers=20) as executor:
        tasks = [(row, schema_map) for _, row in tasks_df.iterrows()]
        
        # tqdm adds the progress bar
        for result in tqdm(executor.map(get_track_bounds, tasks), total=len(tasks), desc="Processing Tracks"):
            if result:
                identifier, bounds = result
                boundary_index[identifier] = bounds
                
                # Periodically save every 500 tracks so we don't lose work on crash
                if len(boundary_index) % 500 == 0:
                    with open(OUTPUT_PATH, 'w') as f:
                        json.dump(boundary_index, f)
    
    # Final Save
    with open(OUTPUT_PATH, 'w') as f:
        json.dump(boundary_index, f)
    
    print(f"Successfully finished! Total indexed tracks: {len(boundary_index)}")

if __name__ == "__main__":
    run_preprocess()