# -*- coding: utf-8 -*-
"""把真实排序/Schur拒绝绑定到同Q研究工作台；未知指标保持null，不改场量。"""
from pathlib import Path
import os,json,hashlib,re
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';target=VAL/'block_linear_research_summary_v2.json';assert not target.exists()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
out=read(VAL/'block_linear_research_summary.json');parentSHA=out['sourceEvidenceSHA256'];assert sha(VAL/'krylov17/wall_velocity_krylov_evidence.mat')==parentSHA
for item in out['studies']:item.update(displayState='拒绝',note='内部GMRES完成或发生breakdown，原组合关系检查拒绝。')
action=VAL/'floworderfactor01/flow_ordering_actions_evidence.mat';binding=loadmat(action,simplify_cells=True,variable_names=['sourceBinding'])['sourceBinding'];assert binding['sourceFileSHA256']==parentSHA
flow=read(VAL/'floworderfactor01/flow_ordering_factor_summary.json');schur=read(VAL/'fieldschurreal01/field_schur_real_summary.json');process=read(VAL/'fieldschurreal01_process.json')
assert process['status']=='terminal' and process['exitCode']==0 and all(process['sourceStableAfter'].values())
before={p.replace('\\','/'):h for p,h in process['sourcesBefore'].items()}
assert before['outputs/v83_validation/krylov17/wall_velocity_krylov_evidence.mat']==parentSHA
assert before['outputs/v83_validation/blockilureal01/steady_graph_packed_evidence.mat']==sha(VAL/'blockilureal01/steady_graph_packed_evidence.mat')
assert before['tests/studyRocketFieldSchurActions.m']==sha(ROOT/'tests/studyRocketFieldSchurActions.m')
assert schur['candidateRejected'] and not schur['preconditionerAvailable'] and schur['errorIdentifier']=='FieldILU:FactorDefect'
# 旧宿主的中文异常消息发生编码损失；数字前缀仍由锁源格式明确给出量数/代数行。
match=re.match(r'^(\d)\D+(\d+)',schur['message']);assert match and int(match[1])==2;row=int(match[2])
patternFile=VAL/'blockfillreal01/expanded_pattern_evidence.mat';s=loadmat(patternFile,simplify_cells=True,variable_names=['pattern','sourceBinding']);assert s['sourceBinding']['sourceFileSHA256']==parentSHA
originalCV=int(s['pattern']['order'][row-1]);g=loadmat(VAL/'flowordergeometry01/flow_ordering_geometry.mat',simplify_cells=True,variable_names=['centroid','report']);nf=int(g['report']['sourceFaces']);layer=(originalCV-1)//nf+1;face=(originalCV-1)%nf+1
out['failedCellDetails']=dict(algebraicRow=row,originalControlVolume=originalCV,originalFace=face,layer=layer,centroidM=g['centroid'][originalCV-1].tolist(),
 numericalFieldGroup='density-and-total-energy',sourcePatternSHA256=sha(patternFile),localInverseDefectValue=None,
 warning='Only the failed 1e-8 gate and location are known. The rejected pivot matrix was not saved; no defect magnitude or condition number is invented.')
def item(label,case,file,residual,state,note,calls=0,papps=0):
 return dict(label=label,caseID=case,linearSystem='C=-Jgraph',arnoldiSteps=0,CApplications=calls,ILUApplications=papps,innerRelativeCResidual=residual,
 basisOrthogonalityError=None,modelRelationRelativeDifference=None,orthogonalityThreshold=1e-8,modelRelationThreshold=1e-8,candidateUsable=False,
 newActualNSEvaluations=0,newProductionSolvedAreaM2=0,evidenceFile=str(file),evidenceSHA256=sha(file),displayState=state,note=note)
out['studies'].append(item('给定来流重排 / 单次块ILU0','floworderfactor01',action,flow['summaries'][0]['relativeCResidual'],'单次解不稳定',
 '仅完成分解与三条固定右端，未开始GMRES。给定来流投影不是当地流线。',5,3))
out['studies'].append(item('三动量—密度/能量Schur','fieldschurreal01',VAL/'fieldschurreal01/field_schur_real_summary.json',None,'分解拒绝',
 f'密度/总能量2量块在代数行{row}的局部逆未过原1e−8；对应原CV {originalCV}、原面{face}、层{layer}。未形成可用预条件器或新流场。'))
out['status']='source-bound-internal-linear-studies-with-flow-order-and-field-Schur-rejections';out['diagnosticBuilderProcessID']=os.getpid();out['builderSHA256']=sha(__file__)
target.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8');print(json.dumps(dict(status=out['status'],studies=len(out['studies']),failedCell=out['failedCellDetails']),ensure_ascii=False))
