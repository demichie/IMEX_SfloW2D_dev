"""Create ensemble inputs while preserving binaries, symlinks and unchanged files."""
import os
from pathlib import Path
import re
import shutil
import tempfile
import pandas as pd

TOKENS = re.compile(r"ENSEMBLE_[A-Za-z][A-Za-z0-9_]*")
TEXT_SUFFIXES = {".inp", ".py", ".sh", ".txt", ".template"}

def replace_strings(working_dir, df, header, i_row):
    values = {"ENSEMBLE_" + str(k): str(df.at[i_row, k]) for k in header}
    for path in sorted(Path(working_dir).iterdir()):
        if path.is_symlink() or not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        raw = path.read_bytes()
        if b"ENSEMBLE_" not in raw:
            continue
        original = raw.decode("utf-8")
        unknown = set(TOKENS.findall(original)) - set(values) - {"ENSEMBLE_vx", "ENSEMBLE_vy"}
        if unknown:
            raise ValueError(f"Unresolved placeholders in {path}: {sorted(unknown)}")
        updated = TOKENS.sub(lambda m: values.get(m.group(), m.group()), original)
        if updated == original:
            continue
        fd, tmp = tempfile.mkstemp(prefix=".ensemble-", dir=str(path.parent))
        try:
            with os.fdopen(fd, "wb") as out:
                out.write(updated.encode("utf-8"))
            os.chmod(tmp, path.stat().st_mode & 0o777)
            os.replace(tmp, path)
        finally:
            if os.path.exists(tmp):
                os.unlink(tmp)

def main():
    df = pd.read_csv("samples.csv")
    for row in range(len(df.index)):
        directory = Path.cwd() / f"ensemble.{row:05d}"
        shutil.copytree("templatedir", directory, symlinks=True)
        try:
            replace_strings(directory, df, df.columns, row)
        except Exception:
            shutil.rmtree(directory)
            raise

if __name__ == "__main__":
    main()
