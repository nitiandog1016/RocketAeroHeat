"""将项目算法清单排为含原生Word公式的报告，并绑定本地来源版本。"""

from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess
import sys
from datetime import datetime, timezone


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--math-package-dir', type=Path)
    parser.add_argument('--math-xsl', type=Path, default=Path('C:/Program Files/Microsoft Office/root/Office16/MML2OMML.XSL'))
    args = parser.parse_args()
    root = args.root.resolve()
    if args.math_package_dir:
        sys.path.insert(0, str(args.math_package_dir.resolve()))
    from latex2mathml.converter import convert
    from lxml import etree
    from docx import Document
    from docx.shared import Inches, Pt, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn

    source = root / 'docs/ENGINEERING_ALGORITHMS_REPORT_2026-11-12.md'
    target = root / 'docs/reports/20261112_algorithm_inventory'
    target.mkdir(parents=True, exist_ok=True)
    text = source.read_text(encoding='utf-8')
    transform = etree.XSLT(etree.parse(str(args.math_xsl)))
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Inches(8.2677), Inches(11.6929)
    sec.top_margin, sec.bottom_margin = Inches(.72), Inches(.68)
    sec.left_margin = sec.right_margin = Inches(.72)
    sec.header_distance = sec.footer_distance = Inches(.3)

    def font(style, size, bold=False):
        s = doc.styles[style]
        s.font.name = 'Calibri'
        s.font.size = Pt(size)
        s.font.bold = bold
        s.font.color.rgb = RGBColor(0, 0, 0)
        rpr = s.element.get_or_add_rPr()
        fonts = rpr.find(qn('w:rFonts'))
        if fonts is None:
            fonts = OxmlElement('w:rFonts'); rpr.append(fonts)
        fonts.set(qn('w:eastAsia'), '宋体' if style == 'Normal' else '微软雅黑')
        return s

    normal = font('Normal', 10.5)
    normal.paragraph_format.line_spacing = 1.17
    normal.paragraph_format.space_after = Pt(5)
    for name, size in [('Title', 22), ('Heading 1', 15), ('Heading 2', 12)]:
        s = font(name, size, name != 'Title')
        s.paragraph_format.space_before = Pt(14 if name != 'Title' else 0)
        s.paragraph_format.space_after = Pt(7)
        s.paragraph_format.keep_with_next = True
    font('Caption', 9)
    # 默认模板的段落或主题边线不进入标题。
    for s in ['Title', 'Heading 1', 'Heading 2']:
        pp = doc.styles[s].element.find(qn('w:pPr'))
        if pp is not None:
            for border in list(pp.findall(qn('w:pBdr'))):
                pp.remove(border)
    header = sec.header.paragraphs[0]
    header.text = 'RocketAeroHeat  工程算法与适用边界'
    header.style = 'Caption'
    footer = sec.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    footer.add_run('2026年11月12日   ')
    field = OxmlElement('w:fldSimple'); field.set(qn('w:instr'), 'PAGE')
    footer._p.append(field)

    table_count, eq_count = 0, 0
    equations = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1; continue
        if line == '$$':
            start = i + 1; i += 1; latex = []
            while i < len(lines) and lines[i].strip() != '$$':
                latex.append(lines[i]); i += 1
            if i == len(lines):
                raise ValueError(f'未闭合公式，行{start}')
            formula = '\n'.join(latex)
            # 该转换器不实现aligned环境；matrix保留同样行列和数学表达。
            math_input = formula.replace(r'\begin{aligned}', r'\begin{matrix}').replace(r'\end{aligned}', r'\end{matrix}')
            mml = convert(math_input)
            if 'merror' in mml:
                raise ValueError(f'MathML错误，行{start}')
            result = transform(etree.fromstring(mml.encode('utf-8')))
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_before = Pt(3)
            p.paragraph_format.space_after = Pt(7)
            p.paragraph_format.line_spacing = 1.0
            p._p.append(result.getroot())
            for run in p._p.findall('.//' + qn('m:r')):
                rpr = OxmlElement('w:rPr')
                ff = OxmlElement('w:rFonts'); ff.set(qn('w:ascii'), 'Cambria Math'); ff.set(qn('w:hAnsi'), 'Cambria Math')
                rpr.append(ff)
                size = OxmlElement('w:sz'); size.set(qn('w:val'), '21'); rpr.append(size)
                run.insert(0, rpr)
            eq_count += 1
            equations.append({'number': eq_count, 'markdownLine': start + 1, 'latex': formula})
            i += 1; continue
        if line.startswith('|'):
            records = []
            while i < len(lines) and lines[i].lstrip().startswith('|'):
                row = [x.strip() for x in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch(r':?-+:?', c.replace(' ', '')) for c in row):
                    records.append(row)
                i += 1
            assert records and all(len(r) == len(records[0]) for r in records)
            t = doc.add_table(rows=0, cols=len(records[0])); t.alignment = WD_TABLE_ALIGNMENT.CENTER
            t.autofit = False
            widths = [1.0, 2.25, 3.25] if len(records[0]) == 3 else [.55, 1.65, 1.7, 2.6]
            if records[0][0] == '索引':
                widths = [.5, 2.15, 3.85]
            elif records[0][0] == '记录':
                widths = [1.65, 2.1, 2.75]
            elif records[0][0] in ['资料', '分支']:
                widths = [1.45, 2.2, 2.85]
            if len(widths) != len(records[0]):
                widths = [6.5 / len(records[0])] * len(records[0])
            for col, width in zip(t.columns, widths):
                col.width = Inches(width)
            pr = t._tbl.tblPr
            borders = OxmlElement('w:tblBorders')
            for side in ['top', 'left', 'bottom', 'right', 'insideH', 'insideV']:
                edge = OxmlElement('w:' + side); edge.set(qn('w:val'), 'single')
                edge.set(qn('w:sz'), '4'); edge.set(qn('w:color'), 'D9D9D9'); borders.append(edge)
            pr.append(borders)
            margin = OxmlElement('w:tblCellMar')
            for side, value in [('top', 70), ('bottom', 70), ('left', 90), ('right', 90)]:
                e = OxmlElement('w:' + side); e.set(qn('w:w'), str(value)); e.set(qn('w:type'), 'dxa'); margin.append(e)
            pr.append(margin)
            for ri, record in enumerate(records):
                row = t.add_row()
                no_split = OxmlElement('w:cantSplit'); row._tr.get_or_add_trPr().append(no_split)
                if ri == 0:
                    row._tr.get_or_add_trPr().append(OxmlElement('w:tblHeader'))
                for ci, value in enumerate(record):
                    c = row.cells[ci]; c.width = Inches(widths[ci]); c.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
                    p = c.paragraphs[0]
                    p.paragraph_format.space_after = Pt(0); p.paragraph_format.line_spacing = 1.1
                    if ri == 0:
                        p.paragraph_format.keep_with_next = True
                    run = p.add_run(value); run.font.size = Pt(9)
                    if ri == 0:
                        run.bold = True; run.font.color.rgb = RGBColor(255, 255, 255)
                    fill = '253D50' if ri == 0 else 'F2F5F7' if ri % 2 == 0 else 'FFFFFF'
                    shading = OxmlElement('w:shd'); shading.set(qn('w:fill'), fill); c._tc.get_or_add_tcPr().append(shading)
            doc.add_paragraph().paragraph_format.space_after = Pt(1)
            table_count += 1
            continue
        if line.startswith('# '):
            doc.add_paragraph(line[2:], 'Title')
        elif line.startswith('## '):
            doc.add_paragraph(line[3:], 'Heading 1')
        elif line.startswith('### '):
            doc.add_paragraph(line[4:], 'Heading 2')
        else:
            p = doc.add_paragraph(line[2:] if line.startswith('- ') else line)
            if line.startswith('- '):
                p.paragraph_format.left_indent = Inches(.1)
                p.paragraph_format.space_after = Pt(3)
                for run in p.runs:
                    run.font.size = Pt(9)
        i += 1

    doc.core_properties.title = 'RocketAeroHeat工程算法与适用边界报告'
    doc.core_properties.subject = '当前源码公式、输入、适用边界与证据状态'
    doc.core_properties.author = 'RocketAeroHeat 项目'
    doc.core_properties.keywords = '工程气动热 公式 适用边界 来源 核验'
    output = target / 'RocketAeroHeat_工程算法与适用边界_20261112.docx'
    doc.save(output)
    # 只读绑定真正用于报告的源码和记录；不盘点或复制大型计算输出。
    paths = [
        'AGENTS.md', 'docs/PROJECT_GOALS.md', str(source.relative_to(root)),
        'airThermalProperties.m', 'phaseAwareAirState.m', 'phaseAwareAirCaloric.m', 'phaseAwareAirReferenceState.m',
        'data/thermo/phase_aware_dry_air_v1.json', 'airReferenceTemperature.m', 'airMeanFreePath.m', 'assessContinuumRegime.m',
        'calcCpMax.m', 'modifiedNewtonianPressure.m', 'solveExternalFlow.m', 'calcBoundaryLayerEdgeState.m',
        'solveNormalShock.m', 'solveThermallyPerfectNormalShock.m', 'solveObliqueShock.m',
        'solveThermallyPerfectObliqueShock.m', 'solveConicalShock.m', 'solveThermallyPerfectConicalShock.m',
        'solvePrandtlMeyerExpansion.m', 'solveThermallyPerfectPrandtlMeyerExpansion.m',
        'solvePhaseAwareWave.m', 'marchPhaseAwareWavePath.m', 'solveSphericalNoseEdge.m',
        'solvePhaseAwareNoseEdge.m', 'estimateBowShockSurface.m', 'calcStagnationHeating.m', 'calcBluntNoseHeating.m',
        'traceSurfaceStreamlines.m', 'calcBoundaryLayerState.m', 'calcTransition.m', 'calcTurbulentBoundaryLayerState.m',
        'calcFixedShapeTurbulentBoundaryLayerState.m', 'marchHeadBoundaryLayerSegment.m', 'calcSurfaceHeating.m',
        'calcPhaseAwareLaminarHeating.m', 'calcKinneyAcreageHeating.m', 'calcWeightedReferenceLength.m',
        'solveAxisymmetricLaminarMassProfile.m', 'solvePhaseAwareLaminarProfile.m', 'solvePhaseAwareNoseHeating.m',
        'solveSphericalNoseLaminarHeating.m', 'solveAxisymmetricBaseCorrelation.m', 'solveMeshBasePressure.m',
        'solveSharpConeBaseFlow.m', 'estimateSonicBoattailPressure.m', 'recoverShockEntropyEdgeState.m',
        'recoverEntropyLayerEdgeState.m', 'calcWallThermalResponse.m', 'solveLayeredWallConduction.m',
        'completeWallEnergyBalance.m',
    ]
    documents = [
        'WORKFLOW_METHOD_SCOPE', 'HYBRID_PHASE_AWARE_AIR_V206', 'HYBRID_PHASE_WAVES_V225',
        'SPHERICAL_NOSE_EDGE', 'HYBRID_PHASE_NOSE_EDGE_V226', 'HYBRID_PHASE_NOSE_HEATING_V227',
        'HYBRID_PHASE_CONTINUATION_V228', 'HYBRID_REFERENCE_ENTHALPY_V190', 'HYBRID_PHASE_LAMINAR_CHAIN_V224',
        'KINNEY_ACREAGE_V104', 'SPHERICAL_NOSE_HEATING', 'PRESSURE_DRIVEN_SURFACE_FLOW_V94',
        'THERMAL_PRESSURE_SURFACE_FLOW_V95', 'LABELLED_PRESSURE_SURFACE_FLOW_V96',
        'NONSIMILAR_BOUNDARY_LAYER_DRAFT', 'V61_ENERGY_AUDIT', 'AXISYMMETRIC_BASE_CORRELATION',
        'MESH_BASE_PRESSURE', 'BASE_RADIAL_HEATING_RESEARCH', 'HYBRID_SU2_PHASE_AIR_V220',
        'HYBRID_PHASE_EULER_NATIVE_V229', 'HYBRID_PHASE_EULER_ROCKET_V230', 'HYBRID_SU2_ROCKET_PHASE_V222',
        'HYBRID_PHASE_PRIMITIVE_RECOVERY_V231',
        'HYBRID_PHASE_EULER_DEVELOPMENT_V232',
        'HYBRID_ROUTE_SELECTION_V189', 'SONIC_BOATTAIL_PRESSURE', 'NASA_BOATTAIL_PRESSURE_COMPARISON',
        'NOSE_AFTERBODY_MATCHING_RESEARCH', 'ENTROPY_MASS_COUPLING', 'SHOCK_LAYER_TRANSPORT_V98_DRAFT',
        'GENERAL_SHOCK_SURFACE_MASS_V99', 'COUPLED_SHOCK_PATCH_V100', 'REGIONAL_HEATING_SOURCE_REVIEW',
        'REGIONAL_ALGORITHM_ROADMAP_2026-10-17', 'REGIONAL_ALGORITHM_RECHECK_V92',
        'ENGINEERING_SOURCE_TRIAGE_V93', 'CYLINDER_CROSSFLOW_SOURCE',
    ]
    paths += ['docs/' + p + '.md' for p in documents]
    paths += ['docs/validation/v230/NEXT_DECISION.md', 'docs/validation/v231/NEXT_DECISION.md', 'docs/validation/v232/NEXT_DECISION.md',
              'harness_tasks/head_closure_review_20261108/SOURCE_REVIEW.md',
              'harness_tasks/integration_recovery_20261110/REPORT.md',
              'harness_tasks/integration_recovery_20261110/candidates/calcBoundaryLayerState.m',
              'harness_tasks/integration_recovery_20261110/candidates/calcTransition.m',
              'harness_tasks/integration_recovery_20261110/candidates/marchHeadBoundaryLayerSegment.m',
              'docs/research/opencode_protuberance_d1372.md',
              'docs/research/d1372_verified_conditions.csv',
              str(Path(__file__).resolve().relative_to(root))]
    items = []
    for relative in sorted(set(paths)):
        file = root / relative
        if not file.is_file():
            raise FileNotFoundError(relative)
        payload = file.read_bytes()
        tracked = subprocess.run(['git', 'ls-files', '--error-unmatch', '--', relative], cwd=root, capture_output=True).returncode == 0
        item = {'path': relative, 'bytes': len(payload), 'sha256': hashlib.sha256(payload).hexdigest(),
                'trackedAtBuild': tracked, 'candidateOnly': relative.startswith('harness_tasks/')}
        # 仅冻结三个短篇来源说明及条件表，不复制并发候选源码或任何工程包。
        if relative in ['harness_tasks/head_closure_review_20261108/SOURCE_REVIEW.md',
                        'harness_tasks/integration_recovery_20261110/REPORT.md',
                        'docs/research/opencode_protuberance_d1372.md',
                        'docs/research/d1372_verified_conditions.csv']:
            frozen = target / 'source_notes' / (Path(relative).parent.name + '_' + Path(relative).name)
            frozen.parent.mkdir(exist_ok=True)
            frozen.write_bytes(payload)
            item['frozenEvidencePath'] = str(frozen.relative_to(root))
        items.append(item)
    report = {'createdUTC': datetime.now(timezone.utc).isoformat(),
              'headCommit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
              'scope': 'read-only source inventory and document production; no new physical computation or qualification',
              'equationBlocks': eq_count, 'tables': table_count, 'sources': items,
              'formulaIndex': equations, 'nativeWordMath': True,
              'docx': {'path': str(output.relative_to(root)), 'bytes': output.stat().st_size,
                       'sha256': hashlib.sha256(output.read_bytes()).hexdigest()}}
    (target / 'source_manifest.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'docx': str(output), 'equations': eq_count, 'tables': table_count, 'sources': len(items)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
