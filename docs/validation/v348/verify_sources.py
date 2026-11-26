"""核对 deg190 工况的一手来源，将发布说明、工作簿行和论文表号分别保存。"""
import argparse
import ast
from hashlib import sha1, sha256
import json
from pathlib import Path
import re
import tempfile
from unicodedata import normalize
from urllib.request import Request, urlopen

import numpy as np
from openpyxl import load_workbook
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[3]
CASE = ROOT / "cases/rfz_rec"
ARCHIVE = ROOT / "outputs/rfz_reference/author_data_20261124"
PREFIX = "Aerodynamics/Stage_1/EUCASS_2023"
SOURCE = ARCHIVE / "files" / PREFIX
RECORD = CASE / "deg190_source_record.json"
ORIGINAL = "612818947ed7679dcc2c4f85200d74c1b45ccc0a"


def read(file):
    return json.loads(file.read_bytes())


def hashes(file):
    data = file.read_bytes()
    git = sha1(b"blob " + str(len(data)).encode("ascii") + b"\0" + data)
    return sha256(data).hexdigest(), git.hexdigest()


def verify(original_tree, source=SOURCE, current_tree=None, paper=None):
    case = read(CASE / "deg190_reference.json")
    original = read(original_tree)
    current = read(current_tree) if current_tree else read(ARCHIVE / "github_tree.json")["tree"]
    assert original["sha"] == ORIGINAL and current["sha"] == case["sourceCommit"]
    assert not original["truncated"] and not current["truncated"]
    orig = {e["path"]: e["sha"] for e in original["tree"]}
    now = {e["path"]: e["sha"] for e in current["tree"]}
    field_hashes = read(ROOT / "docs/validation/v347/evidence/samples.json")["sourceSHA256"]
    payload = sha256()
    for name in ("X", "Y", "Z", "p", "q", "T"):
        stem = "deg190/"+name+"_data.npy"
        path = source / stem
        digest, blob = hashes(path)
        assert digest == field_hashes[name] and orig["Results/"+stem] == now[PREFIX+"/"+stem] == blob
        a = np.load(path, allow_pickle=False)
        assert a.shape == (72706,) and a.dtype == "<f8"
        payload.update(a.tobytes())

    readme = source / "readme"
    assert hashes(readme)[0] == case["readmeSHA256"]
    label, = re.findall(r"Mach\s+([0-9.]+)", readme.read_text(encoding="utf8"))
    workbook_name = "aerodynamic_coefficients.xlsx"
    workbook_path = source / workbook_name
    digest, blob = hashes(workbook_path)
    assert digest == "1d92cbb0c606b68c87560166484eb43be18042d27359591d89e1ae9cf37d9eec"
    assert orig["Results/"+workbook_name] == now[PREFIX+"/"+workbook_name] == blob
    workbook = load_workbook(workbook_path, read_only=True, data_only=True)
    sheet = workbook["Aerodynamic Coefficients"]
    headers = {"A":"Case", "B":"Time (s)", "D":"Velocity (m/s)", "E":"Temperature (K)", "F":"Pressure (Pa)",
               "H":"Dynamic Pressure (Pa)", "I":"Speed of Sound (m/s)", "J":"Mach Number (-)", "K":"Active Engines",
               "M":"Flight Phase", "N":"Geometry", "P":"AOA (deg)", "Q":"PFD (deg)"}
    for col, header in headers.items():
        assert sheet[col+"1"].value == header
    rows = [dict(zip(headers, [sheet[f"{col}{r}"].value for col in headers]), row=r)
            for r in range(2, sheet.max_row+1)]
    row, = [r for r in rows if (r["M"], r["N"], r["P"], r["Q"], r["K"]) ==
            ("Aerodynamic Glide", "RFZ-REC", 190, 0, "-")]
    assert (row["row"], row["A"]) == (4, 2) and round(row["D"]/row["I"], 2) == row["J"]
    workbook.close()

    paper = paper or ROOT / "outputs/rfz_reference/EUCASS2023-018.pdf"
    assert hashes(paper)[0] == case["paperSHA256"]
    page = normalize("NFKC", PdfReader(paper).pages[6].extract_text(extraction_mode="layout"))
    assert "Table 4: Test matrix" in page
    match, = re.findall(r"^\s*(\d+)\s+(1689\.7)\s+(368\.1)\s+(245\.4)\s+190\s+off\s+0\s+15\.26\s*$", page, re.M)
    case_number, speed, pressure, temperature = match
    assert (int(case_number), float(speed), float(pressure), float(temperature)) == (
        3, round(row["D"], 1), row["F"], round(row["E"], 1))

    script = source / "deg190/visualise_data.py"
    script_hash, script_blob = hashes(script)
    assert now[PREFIX+"/deg190/visualise_data.py"] == script_blob
    # 原始版本用 CRLF，后续版本用 LF。两份脚本按各自的 Git 对象独立核验。
    original_script = original_tree.with_name("original_visualise_data.py.txt")
    old_hash, old_blob = hashes(original_script)
    assert old_blob == orig["Results/deg190/visualise_data.py"]
    assert original_script.read_bytes().replace(b"\r\n", b"\n") == script.read_bytes()
    code = ast.parse(script.read_text(encoding="utf8"))
    assert not any(isinstance(n, ast.Name) and n.id == "pdyns" and isinstance(n.ctx, ast.Load) for n in ast.walk(code))
    return dict(schema="rfz-published-source-record-v1", dataset=case["dataset"], sourceCommit=case["sourceCommit"],
                archiveDOI="10.5281/zenodo.8146567", archiveCommit=ORIGINAL,
                caseFileSHA256=hashes(CASE / "deg190_reference.json")[0],
                readme=dict(path=PREFIX+"/readme", sha256=case["readmeSHA256"], nominalMach=float(label)),
                paper=dict(doi=case["paperDOI"], sha256=case["paperSHA256"], page=7, table=4, caseNumber=int(case_number),
                           speedMps=float(speed), staticPressurePa=float(pressure), staticTemperatureK=float(temperature)),
                workbook=dict(path=PREFIX+"/"+workbook_name, sha256=digest, sheet=sheet.title, row=row["row"],
                              caseNumber=row["A"], timeS=row["B"], timeKind="trajectory point",
                              geometry=row["N"], incidenceDeg=row["P"], finDeflectionDeg=row["Q"],
                              speedMps=row["D"], temperatureK=row["E"], pressurePa=row["F"],
                              dynamicPressurePa=row["H"], soundSpeedMps=row["I"], mach=row["J"]),
                visualizationScript=dict(path=PREFIX+"/deg190/visualise_data.py", sha256=script_hash, archiveSHA256=old_hash),
                samples=dict(count=72706, columns=["X","Y","Z","p","q","T"], units=["m","m","m","Pa","W/m2","K"],
                             sourceSHA256=field_hashes, orderedValuesLEFloat64SHA256=payload.hexdigest(),
                             sampleID="one-based source row", pressureDatum="unconfirmed",
                             heatFluxPositiveDirection="unconfirmed", kind="author CFD samples",
                             domain="published half-domain with capped nozzles"))


def restore_inputs(folder):
    """在独立缓存恢复一手资料。文件逐个校验后才发布，中断后可接续。"""
    case = read(CASE / "deg190_reference.json")
    commit = case["sourceCommit"]
    github = "https://api.github.com/repos/rfz-model/RFZ_data/git/trees/"
    raw = "https://raw.githubusercontent.com/rfz-model/RFZ_data/"
    def obtain(name, url, expected=None, limit=1024*1024):
        target = folder / name
        if not target.is_file():
            target.parent.mkdir(parents=True, exist_ok=True)
            with urlopen(Request(url, headers={"User-Agent":"RocketAeroHeat-research"}), timeout=20) as f:
                data = f.read(limit+1)
            assert len(data)<=limit and (not expected or sha256(data).hexdigest()==expected), url
            with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as tmp:
                tmp.write(data)
            Path(tmp.name).replace(target)
        if expected:
            assert hashes(target)[0] == expected, target
        return target
    orig = obtain("original_tree.json",github+ORIGINAL+"?recursive=1")
    tree = obtain("current_tree.json",github+commit+"?recursive=1")
    obtain("original_visualise_data.py.txt",raw+ORIGINAL+"/Results/deg190/visualise_data.py",
           "478f6be6706e2df99bae39301ebe332047c25473fdc261e5bde845b10a803db7")
    base = raw+commit+"/"+PREFIX+"/"
    expected = read(ROOT / "docs/validation/v347/evidence/samples.json")["sourceSHA256"]
    for k,h in expected.items():
        name = "deg190/"+k+"_data.npy"
        obtain(name,base+name,h)
    for name,h in (("readme",case["readmeSHA256"]),
                   ("deg190/visualise_data.py","533e4078accafec306a601c67a09c253a183cccd91ffa7a28dc3d5e7c9c4567a"),
                   ("aerodynamic_coefficients.xlsx","1d92cbb0c606b68c87560166484eb43be18042d27359591d89e1ae9cf37d9eec")):
        obtain(name,base+name,h)
    paper=obtain("EUCASS2023-018.pdf","https://www.eucass.eu/doi/EUCASS2023-018.pdf",case["paperSHA256"],40*1024**2)
    return orig,folder,tree,paper


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument("--original-tree", type=Path, help="使用已保存的作者资料及原始标签 Git 树")
    sources.add_argument("--cache", type=Path, help="在此目录读取或从固定公开来源恢复缺失资料")
    parser.add_argument("--write-new", type=Path, help="仅向新文件写入工况记录")
    args = parser.parse_args()
    result = verify(*restore_inputs(args.cache)) if args.cache else verify(args.original_tree)
    if args.write_new:
        with args.write_new.open("x", encoding="utf8", newline="\n") as f:
            json.dump(result,f,indent=2,ensure_ascii=False,allow_nan=False);f.write("\n")
    else:
        assert result == read(RECORD)
    print(json.dumps(dict(dataset=result["dataset"], paper=result["paper"], workbook=result["workbook"],
                          nominalMach=result["readme"]["nominalMach"], sourceRecordMatched=not bool(args.write_new))))
