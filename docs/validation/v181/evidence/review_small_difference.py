"""限定只读复核：记录早期差异和官方算子检查的适用限制，不从残差猜唯一病因。"""
from pathlib import Path
import hashlib,json,subprocess,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];base=ROOT/'outputs/v181_validation'
source=Path('/opt/rocket-aeroheat/eilmer-v5.0.0/source');start=time.monotonic()
out=base/'difference01';out.mkdir(exist_ok=False)
assert subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()=='f53f4609a0331d48efee69a4e4f3c3598378cc03'
subprocess.run(['git','-C',str(source),'diff','--exit-code'],check=True,stdout=subprocess.DEVNULL)
readings={
 'src/lmr/fluidblock.d':[(1495,1600),(1620,1699)],
 'src/lmr/newtonkrylovsolver.d':[(1874,1897),(2212,2259),(2275,2360),(3194,3238),(3662,3734)],
 'src/lmr/commands/checkjacobian.d':[(145,235)],
 'src/lmr/commands/runsim.d':[(194,260)],
 'src/lmr/makefile':[(401,432)]}
evidence=[]
for name,ranges in readings.items():
    p=source/name;lines=p.read_text().splitlines();dest=out/(p.name+'.excerpt.txt')
    dest.write_text('\n\n'.join('\n'.join(f'{i+1}: {lines[i]}' for i in range(a-1,b)) for a,b in ranges)+'\n')
    evidence.append(dict(path=name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),ranges=ranges,excerpt=dest.name))
real=np.genfromtxt(base/'small01/lmrsim/diagnostics/nk-diagnostics',names=True,deletechars='')
old=np.genfromtxt(ROOT/'outputs/v180_validation/small01/lmrsim/diagnostics/nk-diagnostics',names=True,deletechars='')
variables=['mass-abs','x-mom-abs','y-mom-abs','z-mom-abs','total-energy-abs','tke-abs','omega-abs']
diff=[{key:float(abs(real[key][i]-old[key][i])/max(1.,abs(old[key][i]))) for key in variables} for i in range(8)]
r=dict(passed=True,scope='bounded-source-and-history-review-not-qualification',sourceEvidence=evidence,
       firstRecordedPhysicalResidualsIdentical=all(diff[0][k]==0 for k in variables),perStepResidualScaledDifference=diff,
       realLinearIterations=real['n-iters'].tolist(),complexLinearIterations=old['n-iters'].tolist(),
       sameStepComparisonPassed=False,uniqueCauseIdentified=False,
       selectedNextAction='One source-bound native diagnostic for same-state/order Jv and sparse Jacobian; not another 8-step rollout or a perturbation sweep',
       stockDiagnosticLimitations=['Uses config interpolation order for matrix but phase Jacobian order for Jv',
                                   'Calls Jv with raw configured perturbation, not production computePerturbationSize()',
                                   'Must bind R0, restored flow state and random vector before interpreting agreement'],
       installedNativeJacobianChecker=(source.parent/'install/bin/lmr-check-jacobian').exists(),
       defaultDiagnosticBuildCompilesFullTarget=True,
       sourceSearchErrors=['WSL rg unavailable','Guessed block_jacobian.d, run.d, simstate.d and lmrrun.d paths absent; corrected via observed files'],
       newNativeRuns=0,builds=0,physicalQualification=False,elapsedSeconds=time.monotonic()-start)
(out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:r[k] for k in ['passed','firstRecordedPhysicalResidualsIdentical','realLinearIterations','complexLinearIterations','uniqueCauseIdentified','elapsedSeconds']}))
