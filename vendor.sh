mkdir -p pypgcf/_vendor/pygblocks
curl -L https://files.pythonhosted.org/packages/cc/11/b0accfeb0c97914fcaf1b4e18dad8000c56c6c1b7a9176a2511bc777b341/itaxotools_pygblocks-0.1.0-py3-none-any.whl -o /tmp/itaxotools_pygblocks-0.1.0.whl
python - <<'PY'
import zipfile, os, shutil
whl = "/tmp/itaxotools_pygblocks-0.1.0.whl"
dst = "pypgcf/_vendor/pygblocks"
with zipfile.ZipFile(whl) as z:
    for name in z.namelist():
        if name.startswith("itaxotools/pygblocks/") and not name.endswith("/"):
            rel = name[len("itaxotools/pygblocks/"):]
            out = os.path.join(dst, rel)
            os.makedirs(os.path.dirname(out), exist_ok=True)
            with z.open(name) as src, open(out, "wb") as f:
                f.write(src.read())
PY
touch pypgcf/_vendor/__init__.py

# TODO: Need to create bioconda recipe for itaxotools-pygblocks and add it as a true dependency
# This will be added in pixi.toml, recipe.yaml and be removed from pyproject.toml
