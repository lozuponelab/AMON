import os
import subprocess
from pathlib import Path

KEGG_CLOUD_URLS = {
    "ko":       "https://olucdenver-my.sharepoint.com/:t:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/IQAkn0t-N57yTrlSRoEQsqAHASHSkdJ4Xa-8o12JcprOUkw?download=1",
    "reaction": "https://olucdenver-my.sharepoint.com/:t:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/IQD_NFjvyvncTYU2HaJMme7qAa0FnJRZ2MbdIipgsqTCN8w?download=1",
    "compound": "https://olucdenver-my.sharepoint.com/:t:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/IQCnREW9k4mTR6MgRGPHOOUOAUZ7AX_g1C4XDr7Lv_fWLpc?download=1",
    "pathway":  "https://olucdenver-my.sharepoint.com/:t:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/IQCGf6MllCXKTbVUyJf0V5nvAYsltoH2tfHmnmlGEitUW-U?download=1",
    "enzyme":   "https://olucdenver-my.sharepoint.com/:t:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/IQA7SS36oChSSas37EbaErY9AaOp8YgyVfZSxqqzaNRYrUs?download=1",
}

def get_cache_dir():
    cache_dir = Path(
        os.environ.get(
            "AMON_KEGG_DIR",
            Path.home() / ".cache" / "amon" / "kegg",
        )
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def download_with_wget(url: str, dest: Path):
    print(f"Downloading KEGG file: {url}")
    try:
        subprocess.check_call(
            [
                "wget",
                "--quiet",
                "--content-disposition",
                "--trust-server-names",
                "-O",
                str(dest),
                url,
            ]
        )
    except FileNotFoundError:
        raise RuntimeError(
            "wget is required to download KEGG files from SharePoint.\n"
            "Please install wget or preload KEGG files and set AMON_KEGG_DIR."
        )


def get_kegg_files(force_download=False):
    cache_dir = get_cache_dir()
    paths = {}

    for name, url in KEGG_CLOUD_URLS.items():
        dest = cache_dir / f"{name}.txt"

        if force_download or not dest.exists():
            download_with_wget(url, dest)

        paths[name] = str(dest)

    return paths
