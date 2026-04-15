from __future__ import annotations
from pathlib import Path
from typing import List, Tuple
import time, zipfile, requests
import pandas as pd
from .io_utils import QueryGenome
from .phage_prediction import load_phastest_regions

PHASTEST_POST_URL="https://phastest.ca/phastest_api"
PHASTEST_GET_URL="https://phastest.ca/phastest_api"

def _submit_file(query: QueryGenome, email: str = "", contigs: bool = False) -> dict:
    params={}
    if contigs: params["contigs"]="1"
    if email: params["email"]=email
    with open(query.path,"rb") as fh:
        resp=requests.post(PHASTEST_POST_URL, params=params, data=fh.read(), timeout=120)
    resp.raise_for_status()
    return resp.json()

def _submit_accession(accession: str) -> dict:
    resp=requests.get(PHASTEST_GET_URL, params={"acc": accession}, timeout=120)
    resp.raise_for_status()
    return resp.json()

def _poll_job(job_id: str, poll_seconds: int = 30, max_polls: int = 120) -> dict:
    last={}
    for _ in range(max_polls):
        resp=requests.get(PHASTEST_GET_URL, params={"acc": job_id}, timeout=120)
        resp.raise_for_status()
        last=resp.json()
        if last.get("zip") or last.get("summary") or last.get("url"):
            return last
        time.sleep(poll_seconds)
    return last

def _download_results_zip(zip_url: str, dest: Path) -> Path:
    resp=requests.get(zip_url, timeout=300)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    return dest

def _find_region_tsv(unpack_dir: Path):
    candidates=list(unpack_dir.rglob("*.tsv"))+list(unpack_dir.rglob("*.txt"))
    for p in candidates:
        name=p.name.lower()
        if "region" in name or "summary" in name or "phage" in name:
            return p
    return None

def run_phastest_for_queries(queries: List[QueryGenome], outdir: Path, email: str = "", input_mode: str = "file", contigs: bool = False, poll_seconds: int = 30) -> Tuple[pd.DataFrame, list[dict]]:
    phastest_dir=outdir/"phastest_auto"; phastest_dir.mkdir(parents=True, exist_ok=True)
    region_tables=[]; job_rows=[]
    for query in queries:
        submit = _submit_accession(query.path.stem) if input_mode=="accession" else _submit_file(query, email=email, contigs=contigs)
        job_id=submit.get("job_id","")
        polled=_poll_job(job_id, poll_seconds=poll_seconds)
        zip_url=polled.get("zip","")
        zip_path=None; region_path=None
        if zip_url:
            zip_path=phastest_dir/f"{query.sample}__phastest_results.zip"
            _download_results_zip(zip_url, zip_path)
            unpack_dir=phastest_dir/f"{query.sample}__phastest_unpacked"; unpack_dir.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(zip_path,"r") as zf:
                zf.extractall(unpack_dir)
            region_path=_find_region_tsv(unpack_dir)
            if region_path and region_path.exists():
                try:
                    df=load_phastest_regions(region_path)
                    df["sample"]=query.sample
                    region_tables.append(df)
                except Exception:
                    pass
        job_rows.append({
            "sample":query.sample,"job_id":job_id,"status":polled.get("status", submit.get("status","")),
            "url":polled.get("url",""),"zip_url":zip_url,"zip_file":zip_path.name if zip_path else "",
            "parsed_region_file":region_path.name if region_path else ""
        })
    phage_regions=pd.concat(region_tables, ignore_index=True) if region_tables else pd.DataFrame(columns=["sample","contig","region_id","start","end","length","completeness","score","attl","attr","note"])
    return phage_regions, job_rows
