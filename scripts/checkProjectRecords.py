"""Read-only checks for the short project recovery records; no scientific tests."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
LIMITS = {
    "AGENTS.md": 8192,
    "docs/PROJECT_GOALS.md": 12288,
    "docs/CURRENT_STATE.md": 8192,
    "docs/RECORDS_GUIDE.md": 12288,
}


def check(root):
    errors = []
    sizes = {}
    texts = {}
    documents = list(LIMITS) + [
        "docs/EXECUTION_INDEX.md",
        "docs/validation/record_layout_20261116/README.md",
    ]
    for name in documents:
        path = root / name
        if not path.is_file():
            errors.append(f"missing record: {name}")
            continue
        raw = path.read_bytes()
        sizes[name] = len(raw)
        if name in LIMITS and len(raw) > LIMITS[name]:
            errors.append(f"oversized record: {name}: {len(raw)} > {LIMITS[name]}")
        try:
            content = raw.decode("utf-8-sig")
        except UnicodeDecodeError:
            errors.append(f"not UTF-8: {name}")
            continue
        texts[name] = content
        for link in re.findall(r"\[[^\]\n]*\]\(([^)\n]+)\)", content):
            target = link.strip().strip("<>")
            if urlsplit(target).scheme or target.startswith("#"):
                continue
            local = unquote(target.split("#", 1)[0].split("?", 1)[0])
            if local and not (path.parent / local).exists():
                errors.append(f"missing link: {name}: {target}")

    state = texts.get("docs/CURRENT_STATE.md", "")
    metadata = {}
    patterns = {"last_completed": r"v\d+", "active": r"v\d+|none", "base_commit": r"[0-9a-f]{40}"}
    for key, pattern in patterns.items():
        matches = re.findall(r"<!-- records:" + key + r"=(" + pattern + r") -->", state)
        if len(matches) != 1:
            errors.append(f"expected one valid state marker: {key}")
        else:
            metadata[key] = matches[0]
    completed = metadata.get("last_completed")
    if completed:
        expected = f"validation/{completed}/NEXT_DECISION.md"
        if expected not in state or not (root / "docs" / expected).is_file():
            errors.append("completed version must link its actual NEXT_DECISION.md")
    active = metadata.get("active")
    if active and active != "none":
        expected = f"validation/{active}/PLAN.md"
        if expected not in state or not (root / "docs" / expected).is_file():
            errors.append("active iteration must link its actual PLAN.md")
    for entry in ("AGENTS.md", "docs/PROJECT_GOALS.md"):
        if "CURRENT_STATE.md" not in texts.get(entry, ""):
            errors.append(f"missing canonical recovery link: {entry}")
        if re.search(r"(?:\*\*)?\d{4}-\d{2}-\d{2}\s+v\d+(?:当前|最新)起点", texts.get(entry, "")):
            errors.append(f"version history returned to short entry: {entry}")
    base = metadata.get("base_commit")
    if base:
        result = subprocess.run(["git", "cat-file", "-e", base + "^{commit}"],
                                cwd=root, capture_output=True)
        if result.returncode:
            errors.append("recorded base commit is not available locally")

    manifest_path = root / "docs/validation/record_layout_20261116/migration_manifest.json"
    archived = 0
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for item in manifest["inputs"]:
            path = root / item["snapshot"]
            data = path.read_bytes() if path.is_file() else b""
            if len(data) != item["bytes"] or hashlib.sha256(data).hexdigest() != item["sha256"]:
                errors.append(f"historical snapshot changed or missing: {item['snapshot']}")
            else:
                archived += 1
    else:
        errors.append("migration manifest missing")
    return {"schema": "rah-project-record-check-v1", "passed": not errors,
            "bytes": sizes, "current_markers": metadata,
            "historical_snapshots_verified": archived, "errors": errors,
            "scientific_validation_performed": False,
            "process_liveness_rechecked": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, help="Optional new JSON receipt; will not overwrite")
    args = parser.parse_args()
    result = check(ROOT)
    output = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.report:
        with args.report.open("x", encoding="utf-8") as stream:
            stream.write(output)
    print(output)
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
