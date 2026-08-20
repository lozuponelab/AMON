import os
import subprocess
import tempfile
from pathlib import Path

KEGG_CLOUD_URLS = {
    "ko":       "https://olucdenver-my.sharepoint.com/:t:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/IQAkn0t-N57yTrlSRoEQsqAHASHSkdJ4Xa-8o12JcprOUkw?download=1",
    "reaction": "https://olucdenver-my.sharepoint.com/:t:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/IQD_NFjvyvncTYU2HaJMme7qAa0FnJRZ2MbdIipgsqTCN8w?download=1",
    "compound": "https://olucdenver-my.sharepoint.com/:t:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/IQCnREW9k4mTR6MgRGPHOOUOAUZ7AX_g1C4XDr7Lv_fWLpc?download=1",
    "pathway":  "https://olucdenver-my.sharepoint.com/:t:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/IQCGf6MllCXKTbVUyJf0V5nvAYsltoH2tfHmnmlGEitUW-U?download=1",
    "enzyme":   "https://olucdenver-my.sharepoint.com/:t:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/IQA7SS36oChSSas37EbaErY9AaOp8YgyVfZSxqqzaNRYrUs?download=1",
}

def get_cache_dir() -> Path:
    """
    Resolve the local directory used to cache downloaded KEGG reference files.

    Priority order:
        1. The AMON_KEGG_DIR environment variable, if set — this gives
           explicit control over where files are stored, which is strongly
           recommended on shared or quota-limited systems (e.g. HPC clusters),
           since the default location below can silently consume home
           directory quota.
        2. ~/.cache/amon/kegg as a fallback default. A warning is printed
           when this fallback is used so the side effect (writing to the
           user's home directory) is never silent.

    The directory is created (including any missing parent directories)
    if it does not already exist.

    Returns:
        Path: the resolved cache directory, guaranteed to exist on disk.
    """
    env_dir = os.environ.get("AMON_KEGG_DIR")

    if env_dir:
        # User explicitly configured a cache location - respect it as-is.
        cache_dir = Path(env_dir)
    else:
        # No explicit configuration - fall back to a per-user default,
        # but warn since this write location may not be obvious or
        # desirable (e.g. on systems with small home-directory quotas).
        cache_dir = Path.home() / ".cache" / "amon" / "kegg"
        print(
            f"AMON_KEGG_DIR not set — defaulting to {cache_dir}. "
            f"Set AMON_KEGG_DIR to control this location (recommended on "
            f"shared or quota-limited systems)."
        )

    # Ensure the directory (and any missing parents) exists before any
    # caller tries to read/write files inside it.
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def download_with_wget(url: str, dest: Path) -> None:
    """
    Download a single file from `url` to `dest` using wget, writing safely.

    The file is first downloaded to a temporary file in the same directory
    as `dest`, then atomically renamed into place with os.replace(). This
    guarantees that:
        - A crashed or interrupted download never leaves a corrupted or
          partially-written file at the real destination path.
        - Concurrent processes calling this function for the same `dest`
          (e.g. parallel pipeline runs sharing one cache directory) cannot
          race to write the same file and produce a corrupted result -
          each writes its own temp file, and only a complete download is
          ever renamed into the final location.

    Args:
        url: Source URL to download from.
        dest: Final destination path for the downloaded file.

    Raises:
        RuntimeError: if wget is not installed / not found on PATH.
        subprocess.CalledProcessError: if wget runs but exits non-zero
            (e.g. network failure, invalid URL, server error).
    """
    print(f"Downloading KEGG file: {url}")

    # Create the temp file in the same directory as `dest` (not /tmp) so
    # that os.replace() below is guaranteed to be an atomic rename rather
    # than a cross-filesystem copy, which is not atomic and could still
    # leave a partial file if interrupted.
    tmp_fd, tmp_path = tempfile.mkstemp(dir=dest.parent, suffix=".tmp")
    os.close(tmp_fd)  # we only needed mkstemp to reserve a unique filename

    try:
        subprocess.check_call(
            [
                "wget",
                "--quiet",
                "--content-disposition",
                "--trust-server-names",
                "-O",
                tmp_path,
                url,
            ]
        )
        # Only reached if wget succeeded - safe to publish the download
        # under its real filename.
        os.replace(tmp_path, dest)  # atomic on the same filesystem
    except FileNotFoundError:
        # subprocess couldn't find the `wget` executable at all.
        raise RuntimeError(
            "wget is required to download KEGG files from SharePoint.\n"
            "Please install wget or preload KEGG files and set AMON_KEGG_DIR."
        )
    finally:
        # Clean up the temp file if it's still present - this covers both
        # the wget-not-found case and any wget failure/exception, so we
        # never leave orphaned .tmp files behind in the cache directory.
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def get_kegg_files(force_download: bool = False) -> dict:
    """
    Ensure all required KEGG reference files are present locally, downloading
    any that are missing (or all of them, if force_download=True), and
    return their local paths.

    Files are cached under the directory resolved by get_cache_dir(), keyed
    by the names defined in KEGG_CLOUD_URLS (e.g. "ko", "reaction",
    "compound", "pathway", "enzyme"). Once downloaded, a file is reused on
    subsequent calls rather than re-downloaded, unless force_download=True.

    Args:
        force_download: if True, re-download every file even if a cached
            copy already exists locally. Useful for refreshing stale KEGG
            data, since there is currently no automatic staleness/version
            check on cached files.

    Returns:
        dict: mapping of file name (e.g. "ko") -> local file path (str)
        for every entry in KEGG_CLOUD_URLS.
    """
    cache_dir = get_cache_dir()
    paths = {}

    for name, url in KEGG_CLOUD_URLS.items():
        dest = cache_dir / f"{name}.txt"

        # Skip re-downloading if we already have a cached copy, unless
        # the caller explicitly asked to force a fresh download.
        if force_download or not dest.exists():
            download_with_wget(url, dest)

        paths[name] = str(dest)

    return paths