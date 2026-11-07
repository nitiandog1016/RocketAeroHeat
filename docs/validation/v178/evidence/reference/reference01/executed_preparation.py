"""整理固定GDTK来源的CUBRC Run 14实验参考；不生成网格、不调用求解器。

模块一校验原始Git对象与SHA256，模块二按原表测站及Btu_IT转换，模块三
分开保存报告条件、代码实际物性和派生条件，模块四保留缺测/注释/越域点。
该程序准备具名层流简化对照，不授予真实20度箭体、RANS或任何物理资格。
"""
from __future__ import print_function
import argparse
import csv
from decimal import Decimal, getcontext
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import sys
import time
import traceback


ROOT = Path(__file__).resolve().parents[1]
COMMIT = "f53f4609a0331d48efee69a4e4f3c3598378cc03"
REPOSITORY = "gdtk-uq/gdtk"
BASE_URL = "https://raw.githubusercontent.com/{}/{}/".format(REPOSITORY, COMMIT)
PAPER_URL = "https://chimeracfd.com/professional/public/aiaa2004_0529.pdf"
PINNED_BLOBS = {
    "COPYING": "6018ea9c93f01ef16b26c686b4d32fec05e8ee10",
    "gpl.txt": "94a9ed024d3859793618152ea559a168bbcbb5e2",
    "contributors.txt": "8c666ac184a2b53b541a2fe2b5aec56b9e4df9e5",
    "cyl-flare.lua": "cfb788a0f798f69f1186a9c2b0d431812ad1764b",
    "idealN2.inp": "27bfc07e575be53890e2af07097c411635361105",
    "cylinder-extended-flare-heat-transfer.data": "07907c89f9b2713bbf645d61627c2e6076cd148e",
    "cylinder-extended-flare-pressure.data": "9b10d07833306749869f3e2fce798b91531f0434",
    "N2.lua": "098972692bb9d22a1f3258c54ce39b4baf5a6ab6",
    "prep_gas.lua": "2475a96fd02805574971b65933f1461940477cee",
    "ideal_gas.d": "f27e1abf8ddde203860873c4b0a946d485cf4d3c",
    "physical_constants.d": "dcb131945f7d455fe109ea8cdfb32231da739fb1",
    "cea_viscosity.d": "9797801418a196b1a1769440a171a8bd432b7c6f",
    "cea_therm_cond.d": "fd6a5807c7b1bc055541b086fd9a377a8b14e4f5",
}
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
getcontext().prec = 40
INCH_M = Decimal("0.0254")
FOOT_M = Decimal("0.3048")
BTU_IT_J = Decimal("1055.05585262")
POUND_KG = Decimal("0.45359237")
STANDARD_GRAVITY = Decimal("9.80665")
BTU_FLUX_SI = BTU_IT_J / FOOT_M ** 2
PSI_PA = POUND_KG * STANDARD_GRAVITY / INCH_M ** 2


def must(condition, text):
    if not condition:
        raise ValueError(text)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def text_number(value):
    return "" if value is None else format(value, "f")


def scalar(source, name):
    found = re.findall(r"(?m)^\s*" + re.escape(name) + r"\s*=\s*(" + NUMBER + r")", source)
    must(len(found) == 1, "原源码标量缺失或不唯一：" + name)
    return float(found[0])


def value_block(source, name):
    match = re.search(r"db\.N2\." + re.escape(name) + r"\s*=\s*\{(.*?)\n\}", source, re.S)
    must(match is not None, "N2数据块缺失：" + name)
    return match.group(1)


def assignments(block, keys):
    result = {}
    for name in keys:
        found = re.findall(r"\b" + re.escape(name) + r"\s*=\s*(" + NUMBER + r")", block)
        must(len(found) == 1, "物性系数缺失或不唯一：" + name)
        result[name] = float(found[0])
    return result


def verify_sources(dataset):
    """模块一：仅读固定小文件，离线核对Git blob与下载清单，不访问网络。"""
    manifest_path = dataset / "source_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    must(manifest["commit"] == COMMIT and manifest["repository"] == REPOSITORY, "上游版本身份不匹配")
    entries = {item["name"]: item for item in manifest["files"]}
    must(len(entries) == len(manifest["files"]) and set(entries) == set(PINNED_BLOBS), "原始文件集合不匹配")
    locks = {str(manifest_path.resolve()): sha(manifest_path)}
    texts = {}
    for name, expected_blob in PINNED_BLOBS.items():
        item = entries[name]
        must(item["localPath"] == "raw/" + name, "原始路径不能改指其他文件")
        must(item["url"] == BASE_URL + item["upstreamPath"], "下载URL必须绑定唯一commit")
        path = dataset / "raw" / name
        content = path.read_bytes()
        must(len(content) <= 512 * 1024 and len(content) == item["bytes"], "小文件字节数不符：" + name)
        git_blob = hashlib.sha1(("blob %d\0" % len(content)).encode("ascii") + content).hexdigest()
        must(git_blob == expected_blob == item["gitBlobSHA1"], "固定Git blob改变：" + name)
        digest = hashlib.sha256(content).hexdigest()
        must(digest == item["sha256"], "原始SHA256改变：" + name)
        locks[str(path.resolve())] = digest
        texts[name] = content.decode("utf-8")
    return manifest, entries, texts, locks


def measurement_rows(source, kind, entry, geometry):
    """模块二：保留注释传感器及缺测符号；越域不是删点或自动映射到末面。"""
    lines = source.splitlines()
    header = [line for line in lines if "Run 14" in line and "Sensor" in line]
    must(len(header) == 1, "实验运行列表不唯一")
    runs = [int(x) for x in re.findall(r"Run\s+(\d+)", header[0])]
    must(runs == [2, 3, 8, 9, 11, 12, 13, 14], "原表运行列发生变化")
    length_match = re.findall(r"#\s*L\s*=\s*(" + NUMBER + r")\s+inches", source)
    must(length_match == ["4.004"], "原测站归一化长度必须保持4.004英寸")
    length = Decimal(length_match[0]) * INCH_M
    factor = BTU_FLUX_SI if kind == "heatFlux" else PSI_PA
    original_unit = "Btu/(ft^2*s)" if kind == "heatFlux" else "psia"
    unit = "W/m^2" if kind == "heatFlux" else "Pa"
    sensor_pattern = r"^(HT\d+|P\d+)\s"
    rows = []
    for line_number, raw in enumerate(lines, 1):
        stripped = raw.strip()
        commented = stripped.startswith("#")
        line = stripped[1:].lstrip() if commented else stripped
        if not re.match(sensor_pattern, line):
            continue
        parts = line.split()
        must(len(parts) == 2 + len(runs), "测站列数异常：" + parts[0])
        sensor, ratio_text, value_text = parts[0], parts[1], parts[2 + runs.index(14)]
        ratio = Decimal(ratio_text)
        x = ratio * length
        value = None if set(value_text.lower()) == {"x"} else Decimal(value_text)
        must(value is None or value.is_finite(), "实验数值必须有限或原表缺测")
        inside = Decimal("0") <= x <= Decimal(str(geometry["domainEndX_m"]))
        if not inside:
            domain_status, region, radius = "outside-upstream-example-domain", "outside", ""
        elif x <= Decimal(str(geometry["cylinderEndX_m"])):
            domain_status, region = "inside-upstream-example-domain", "cylinder"
            radius = format(geometry["radius_m"], ".17g")
        else:
            domain_status, region = "inside-upstream-example-domain", "flare"
            radius = format(geometry["radius_m"] + (float(x) - geometry["cylinderEndX_m"]) *
                            math.tan(math.radians(geometry["flareAngle_deg"])), ".17g")
        record_status = "source-commented" if commented else "source-active"
        measurement_status = "missing-original-token" if value is None else "numeric-original-token"
        if commented:
            comparison = "excluded-by-source-comment-preserved"
        elif value is None:
            comparison = "unavailable-original-missing"
        elif not inside:
            comparison = "unmapped-outside-current-domain"
        else:
            comparison = "available-for-comparison-not-qualified"
        rows.append(dict(source_line=line_number, sensor=sensor, run=14, source_x_over_L=ratio_text,
                         source_value=value_text, source_unit=original_unit,
                         normalization_length_m=text_number(length), x_from_leading_edge_m=text_number(x),
                         project_X_m=text_number(-x), radius_m_on_upstream_geometry=radius, geometry_region=region,
                         value_SI=text_number(value * factor if value is not None else None), unit_SI=unit,
                         record_status=record_status, measurement_status=measurement_status,
                         domain_status=domain_status, comparison_status=comparison,
                         uncertainty_SI="", physical_qualification="not-assessed",
                         source_sha256=entry["sha256"], source_url=entry["url"]))
    return rows, length


def make_conditions(texts, entries):
    """模块三：解析固定源码真实分派；不把输入字样误当最终选中的物性。"""
    setup, species = texts["cyl-flare.lua"], texts["N2.lua"]
    vals = {key: scalar(setup, key) for key in ["p_inf", "u_inf", "T_inf", "T_vib", "T_wall", "L1", "L2", "R1", "alpha"]}
    geometry = dict(radius_m=vals["R1"] * .001, cylinderEndX_m=vals["L1"] * .001,
                    domainEndX_m=vals["L2"] * .001, flareAngle_deg=vals["alpha"],
                    source="cyl-flare.lua", physicalGeometry="axisymmetric hollow cylinder plus extended flare",
                    coordinate="原前缘x=0、流动沿+x；项目X=-x，原面测站不移动或裁到末面",
                    experimentalNormalizationLength_m=float(Decimal("4.004") * INCH_M),
                    geometricLabelsAreNotFlowSeparationLabels=True,
                    baseOrInteriorHeatingReferenceProvided=False)
    must(vals["alpha"] == 30 and vals["L1"] == 101.7 and vals["L2"] == 220, "具名几何发生变化")
    molar_mass = assignments(value_block(species, "M"), ["value"])["value"]
    gamma = assignments(value_block(species, "gamma"), ["value"])["value"]
    ru_match = re.search(r"R_universal\s*=\s*(" + NUMBER + r")", texts["physical_constants.d"])
    must(ru_match is not None, "固定引擎通用气体常数缺失")
    ru = float(ru_match.group(1))
    gas_r = ru / molar_mass
    cp = gamma * gas_r / (gamma - 1)
    mu_match = re.search(r"segment0\s*=\s*\{(.*?)\}", value_block(species, "ceaViscosity"), re.S)
    must(mu_match is not None, "N2首段CEA黏度缺失")
    cea = assignments(mu_match.group(1), ["T_lower", "T_upper", "A", "B", "C", "D"])
    sutherland = assignments(value_block(species, "sutherlandVisc"), ["mu_ref", "T_ref", "S"])
    options = texts["idealN2.inp"]
    prandtl = scalar(options, "Prandtl")
    must("viscosity = 'CEA'" in options and "thermal_conductivity = 'CEA'" in options, "示例物性选项改变")
    prep = texts["prep_gas.lua"]
    must(prep.index("if optsTable and optsTable.Prandtl then") < prep.index("elseif (optsTable and optsTable.thermal_conductivity"), "Prandtl优先分派改变")
    must("Q.k = _Cp*Q.mu/_Prandtl;" in texts["ideal_gas.d"], "实际固定Pr导热式改变")
    must("if ( T < _T_lowest ) {\n            return _mu_lowest;" in texts["cea_viscosity.d"], "原生CEA低温处理改变")
    p, t, u = vals["p_inf"], vals["T_inf"], vals["u_inf"]
    rho = p / gas_r / t
    effective_t = max(t, cea["T_lower"])
    must(effective_t <= cea["T_upper"], "当前只解析具名工况的首段CEA，不通用外推")
    mu = math.exp(cea["A"] * math.log(effective_t) + cea["B"] / effective_t + cea["C"] / effective_t ** 2 + cea["D"]) * 1e-7
    mu_suth = sutherland["mu_ref"] * (t / sutherland["T_ref"]) ** 1.5 * (sutherland["T_ref"] + sutherland["S"]) / (t + sutherland["S"])
    length = geometry["experimentalNormalizationLength_m"]
    return dict(schema="RAH.CUBRC.Run14.Reference.v1", caseID="CUBRC-Run14-hollow-cylinder-30deg-extended-flare",
                status="reference-preparation-only; laminar ideal-N2 control; not experimental qualification",
                repository=REPOSITORY, commit=COMMIT, release="v5.0.0", geometry=geometry,
                reported=dict(Mach=10.30, velocity_m_s=u, pressure_Pa=p, temperature_K=t,
                              density_kg_m3=8.810e-4, vibrationTemperature_K=vals["T_vib"], wallTemperature_K=vals["T_wall"],
                              Reynolds=None, ReynoldsNote="本次取得的原表1/配置未独立给出Re；不把派生Re写为实测",
                              source=dict(title="MacLean and Holden, AIAA 2004-0529, Table 1", url=PAPER_URL,
                                          sourceSHA256=None, paperCopiedIntoRepository=False,
                                          extraction="人工核对网页可读表1；未建立整论文原始文件SHA",
                                          codeAlsoRecords=["p", "T", "U", "Tw", "Tv"])),
                codeDerived=dict(independentInputs=["p_inf", "T_inf", "u_inf", "idealN2 model selected by fixed source"],
                                 gasConstant_J_kgK=gas_r, universalGasConstant_J_molK=ru, molarMass_kg_mol=molar_mass,
                                 gamma=gamma, cp_J_kgK=cp, density_kg_m3=rho,
                                 Mach=u / math.sqrt(gamma * gas_r * t),
                                 caloricallyPerfectTotalTemperature_K=t + u * u / (2 * cp),
                                 viscosity_Pas=mu, conductivity_W_mK=cp * mu / prandtl,
                                 Re_per_m=rho * u / mu, Re_L_originalTable=rho * u * length / mu,
                                 Re_D=rho * u * 2 * geometry["radius_m"] / mu,
                                 densityRelativeDifferenceFromReported=rho / 8.810e-4 - 1,
                                 derivationIsMeasurement=False, numericalRuntimeNotExecuted=True),
                transport=dict(requested=dict(viscosity="CEA", thermal_conductivity="CEA", Prandtl=prandtl),
                               actual=dict(viscosity="CEA with constant endpoint continuation",
                                           thermalConductivity="cp*mu/Pr; Prandtl option takes precedence over CEA conductivity",
                                           Prandtl=prandtl, ceaFirstSegment=cea, actualEvaluationTemperature_K=effective_t,
                                           freestreamBelowFittedRange=t < cea["T_lower"], lowTemperatureAccuracyQualified=False),
                               sourceBindings={key: entries[key]["sha256"] for key in ["N2.lua", "prep_gas.lua", "ideal_gas.d", "physical_constants.d", "cea_viscosity.d"]},
                               sutherlandDiagnostic=dict(selected=False, coefficients=sutherland, viscosity_Pas=mu_suth,
                                                         Re_L_originalTable=rho * u * length / mu_suth,
                                                         relativeDifferenceFromActualViscosity=mu_suth / mu - 1,
                                                         applicabilityQualified=False)),
                boundary=dict(inlet="uniform supersonic p/T/U at real leading edge and outer boundary",
                              upstreamBoundaryLayer="developed internally from physical leading edge; no supplied profile",
                              wall="no-slip fixed Tw=295.2 K", turbulence="none; laminar",
                              downstream="upstream example uses fixed outside pressure p_inf/5", downstreamPressure_Pa=p / 5,
                              baseInteriorWakeQualification=False, gasVibration="reported Tv=2467 K retained; upstream ideal-N2 example ignores it",
                              vibrationEnergySubstitutionPerformed=False, wallTemperatureSolved=False),
                units=dict(BtuConvention="International Table Btu (Btu_IT), explicit project conversion convention",
                           originalBtuVariantSpecified=False, Btu_IT_J=str(BTU_IT_J), foot_m=str(FOOT_M), inch_m=str(INCH_M),
                           heatFluxMultiplier_W_m2=str(BTU_FLUX_SI), pressureMultiplier_Pa=str(PSI_PA),
                           originalNormalizationInch="4.004", measurementXUsesOriginalLength=True,
                           originalCoefficientRoundoffOrUncertaintyImproved=False),
                limitations=["当前原始数据不提供逐点实验/读数不确定度；unknown保持未知", "报告p/T/rho并非完全相容，报告量不覆盖成派生量",
                             "CEA在来流温度返回200K端点常值；可复现源码不等于低温适用性通过",
                             "Prandtl优先，不能把该例说成同时采用CEA导热拟合", "忽略振动是具名简化模型差异，不获非平衡资格",
                             "两处活跃压力测站超出示例尾端，保留而不插值、裁点或外推几何", "内嵌成熟引擎的教程配置不等于本项目已运行或已验证",
                             "本轴对称层流候选不授真实20度三维RANS、底面、突起及全面积热环境资格"],
                experimentalQualification=False, newQualifiedArea_m2=0, nativeTasks=0, modifiesDefaultSolver=False)


def write_csv(path, rows):
    must(bool(rows), "不允许空实验测站表")
    with path.open("x", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    """模块四：输出独立新目录，保留准确执行源/输入锁/失败，不覆盖既有证据。"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, default=ROOT / "tests/data/eilmer_cylinder_flare")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    must(not args.output.exists(), "证据目录已存在，拒绝覆盖")
    args.output.mkdir(parents=True)
    shutil.copyfile(str(Path(__file__).resolve()), str(args.output / "executed_preparation.py"))
    start = time.perf_counter()
    locks = {str(Path(__file__).resolve()): sha(__file__)}
    report = dict(passed=False, commit=COMMIT, nativeTasks=0, newQualifiedArea_m2=0, sourceImmutable=True)
    try:
        manifest, entries, texts, source_locks = verify_sources(args.dataset)
        locks.update(source_locks)
        conditions = make_conditions(texts, entries)
        qname = "cylinder-extended-flare-heat-transfer.data"
        pname = "cylinder-extended-flare-pressure.data"
        heat, length = measurement_rows(texts[qname], "heatFlux", entries[qname], conditions["geometry"])
        pressure, plength = measurement_rows(texts[pname], "pressure", entries[pname], conditions["geometry"])
        must(length == plength == Decimal("0.1017016"), "两表原始长度必须独立一致")
        must(len(heat) == 55 and len(pressure) == 26, "包含原注释行的测站数不符")
        must(sum(r["record_status"] == "source-active" for r in heat) == 54, "54条活跃热流记录必须全部保留")
        must(sum(r["record_status"] == "source-active" for r in pressure) == 24, "24条活跃压力记录必须全部保留")
        outside = [r["sensor"] for r in pressure if r["record_status"] == "source-active" and r["geometry_region"] == "outside"]
        must(outside == ["P30", "P32"], "活跃越域测点必须保留原身份")
        must(next(r for r in heat if r["sensor"] == "HT56")["value_SI"] == "", "原缺失热流不能补值")
        must(next(r for r in pressure if r["sensor"] == "P33")["record_status"] == "source-commented", "有限注释点不能偷偷升级")
        must(next(r for r in heat if r["sensor"] == "HT36")["source_value"] == "21.01", "峰值原token改变")
        write_csv(args.output / "heat_flux_SI.csv", heat)
        write_csv(args.output / "pressure_SI.csv", pressure)
        conditions["measurements"] = dict(heatFile="heat_flux_SI.csv", pressureFile="pressure_SI.csv",
                                          heatRowsIncludingSourceComments=len(heat), pressureRowsIncludingSourceComments=len(pressure),
                                          activeHeatRows=54, activePressureRows=24, activeOutsidePressureSensors=outside,
                                          noInterpolatedOrRepairedMeasurements=True, allUncertaintiesUnknown=True)
        conditions["sourceManifestSHA256"] = sha(args.dataset / "source_manifest.json")
        write_json(args.output / "conditions.json", conditions)
        report.update(passed=True, rawFiles=len(entries), rawBytes=sum(e["bytes"] for e in entries.values()),
                      heatRows=len(heat), pressureRows=len(pressure), activeHeatRows=54, activePressureRows=24,
                      activeOutsidePressureSensors=outside, originalLength_m=float(length),
                      belowCEARange=conditions["transport"]["actual"]["freestreamBelowFittedRange"],
                      constantPrandtlActuallySelected=True, densityDifference=conditions["codeDerived"]["densityRelativeDifferenceFromReported"],
                      unknownUncertaintiesPreserved=True, commentsAndMissingTokensPreserved=True)
    except Exception as error:
        report["error"] = str(error)
        report["traceback"] = traceback.format_exc()
    finally:
        stable = {path: sha(path) == digest for path, digest in locks.items()}
        report["inputsStable"] = all(stable.values())
        report["passed"] = report["passed"] and report["inputsStable"]
        report["elapsedSeconds"] = time.perf_counter() - start
        write_json(args.output / "input_locks.json", dict(sha256=locks, stable=stable))
        report["outputsSHA256"] = {p.name: sha(p) for p in args.output.iterdir() if p.is_file()}
        write_json(args.output / "report.json", report)
        print(json.dumps({k: report.get(k) for k in ["passed", "heatRows", "pressureRows", "inputsStable", "elapsedSeconds", "error"]}, ensure_ascii=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
