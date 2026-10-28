"""将已核对的有限步诊断接回既有五步骤工作台；不把小域状态写入整箭云图。"""
from pathlib import Path
import json,shutil,html,hashlib
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v84_validation'
def main():
    folder=OUT/'workbench01';assert not folder.exists()
    source=OUT/'workbench_base';shutil.copytree(source,folder)
    binding=json.loads((OUT/'view_source_binding.json').read_text('utf8'))
    # 原导出器生成真实接受状态；本步骤只扩展计算/结果侧栏，不改原面数据。
    shutil.copy2(ROOT/'web/rocket-wall-velocity.js',folder/'rocket-wall-velocity.js')
    audit=json.loads((OUT/'full02/independent_audit.json').read_text('utf8'))
    local=json.loads((OUT/'independent_audit_patch03.json').read_text('utf8'))['patchComparison']
    full=json.loads((OUT/'full02/full_comparison.json').read_text('utf8'))
    def summary(b):
        h=b['history'];return f"{b['method']}：{len(h)}次接受，末态RMS {h[-1]['RMS']:.3f}、最大 {h[-1]['maximum']:.3f}" if h else f"{b['method']}：0次接受"
    finite={'sourceLabel':'Q22→Q24有限步机制；从相同Q24开始整箭对照',
        'mechanismLabel':'内部对流曲率、波系采样变化与病态大更新共同作用',
        'patchComparisonLabel':'48核心、同0.0999约束六步：最小守恒变量幅度RMS 472952.097，最小原始变量幅度180922.009；后四步仍停滞',
        'fullComparisonLabel':'；'.join(summary(b) for b in audit['branches']),
        'scopeLabel':'局部诊断不是整箭收敛证明；当前云图为实际接受检查点，全部为研究值',
        'reportURL':'finite-step-diagnostics.html','newProductionSolvedAreaM2':0,'physicalSolution':False,
        'fullComparison':audit,'localComparison':local}
    d=json.loads((folder/'wall-velocity-study-data.json').read_text('utf8'));d['finiteStepResearch']=finite
    d['cases'][2]['label']='v84实际接受状态 · '+('原方向选择' if binding['selectedMethod']=='baseline' else '最小状态选择')
    raw=json.dumps(d,ensure_ascii=False,separators=(',',':'))
    (folder/'wall-velocity-study-data.json').write_text(raw,encoding='utf8');(folder/'wall-velocity-study-data.js').write_bytes(('window.ROCKET_WALL_VELOCITY_STUDY='+raw+';\n').encode('utf8'))
    rows=[]
    for branch in audit['branches']:
        for r in branch['history']:
            vals=[branch['method'],r['iteration'],r['alpha'],r['RMS'],r['maximum'],r['internalC'],r['originalNSDirection'],r['temperatureRange'],r['spatialBalance'],r['steadyBoundaryImbalance']]
            rows.append('<tr>'+''.join('<td>'+html.escape(f'{v:.7g}' if isinstance(v,float) else str(v))+'</td>' for v in vals)+'</tr>')
    stops=[]
    for branch in full['branches']:
        last=branch['steps'][-1];stops.append(f"{branch['method']}：{last['stopReason']}；最后内部C={last['internalCRelativeResidual']:.8g}，选择空间最低缺陷={last['selection']['minimumModelDefectInSpan']:.8g}。")
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>有限步诊断 · RocketAeroHeat</title>
<style>body{font:16px/1.65 system-ui,sans-serif;color:#183044;max-width:1200px;margin:32px auto;padding:0 24px;background:#f7f9fb}h1{font-size:26px}img{width:100%;background:white}table{border-collapse:collapse;background:white;width:100%;font-size:13px}td,th{border:1px solid #cbd5df;padding:8px;text-align:right}th{background:#e7eef4}a{color:#126caa}.table{overflow:auto}.note{background:#fff4d9;padding:16px;border-left:4px solid #c38a24}</style>
<a href="index.html">返回五步骤科研工作台</a><h1>真实20°箭体：有限步为何偏离正确导数</h1>
<p class="note">已定位机制并实现研究方向选择；局部六步改善明显，但后续仍停滞。整箭是否有效必须依据下表，不能由小域结果外推。全部结果未达到区域物理验证，新增正式求解面积为0。</p>
<p>原NS方程、边界、正性和80–2000 K温域保持。内部线性目标1e−3、原NS方向0.1、最大残差2e−6 /s和定常边界不平衡5e−10分别判断。</p>
<h2>机制与六步小复现</h2><p>总偏差扣除内部对流后仅余0.07922%。2980面采样区域改变；未改变区域时仍有平滑曲率。真实面238752的右激波在原接受增量15.6227%处跨过采样位置。</p>
<img src="finite_step_comparison.png" alt="六步残差与真实面采样切换对照"><p>48核心、312含邻域单元精确复现原NS内部行。固定同一Q24外侧状态；不是整箭替代边界。灰线为精确Newton局部参考，不能称为整箭LSQR的缩小实现。紫线与蓝线采用相同0.0999模型约束和六步预算；守恒变量度量20.330秒/2931次原残差，原始变量度量19.456秒/2954次。没有宣称计时加速。蓝线第3至6步RMS仅再下降0.01191%。</p>
<h2>整箭实际连续接受记录</h2><div class="table"><table><thead><tr><th>方法</th><th>本段步次</th><th>alpha</th><th>原RMS /s</th><th>最大残差 /s</th><th>内部C</th><th>原NS方向</th><th>全域温度/K</th><th>空间账</th><th>定常边界不平衡</th></tr></thead><tbody>'''+''.join(rows)+'''</tbody></table></div>
<p>'''+html.escape(' '.join(stops))+'''</p><p>两组每个新状态采用相同三段128步及15个完整模型列，最多17候选。中断前已完成检查点被保留；未保存的中断计算量未知，因此不宣称端到端加速。</p>
<p><a href="finite_step_comparison.pdf">下载科研图 PDF</a> · <a href="finite_step_record.json">下载完整诊断记录</a></p></html>'''
    (folder/'finite-step-diagnostics.html').write_text(page,encoding='utf8')
    (folder/'finite_step_record.json').write_text(json.dumps(finite,ensure_ascii=False,indent=2),encoding='utf8')
    for ext in ['png','pdf']:shutil.copy2(OUT/f'finite_step_comparison.{ext}',folder/f'finite_step_comparison.{ext}')
    same=[]
    for file in sorted((source/'cases').glob('*/results.json')):
        counterpart=folder/file.relative_to(source);assert file.read_bytes()==counterpart.read_bytes();same.append(str(file.relative_to(source)))
    result={'status':'exported-diagnosis-in-existing-five-step-workbench','originalFacetDatasetsByteIdentical':same,
        'cloudState':binding,'physicalSolution':False,'newProductionSolvedAreaM2':0}
    (folder/'finite_step_view_binding.json').write_text(json.dumps(result,indent=2),encoding='utf8');print(json.dumps(result))
if __name__=='__main__':main()
