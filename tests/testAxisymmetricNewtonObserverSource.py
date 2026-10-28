"""移除明确列出的观察语句后，独立记录器必须与默认求解器的算法正文完全一致。"""
from pathlib import Path
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
core=ROOT/'regionalFlow/axisymmetric/advanceAxisymmetricSteady.m'
observer=core.with_name('advanceAxisymmetricSteadyObserved.m')
original=core.read_text(encoding='utf8')
observed=observer.read_text(encoding='utf8')


def replace_once(text,old,new):
    assert text.count(old)==1,old
    return text.replace(old,new,1)


observed=replace_once(observed,
    'function [result,diagnostics]=advanceAxisymmetricSteadyObserved(Q,mesh,gas,boundary,options,diagnosticOptions)',
    'function result=advanceAxisymmetricSteady(Q,mesh,gas,boundary,options)')
observed=replace_once(observed,'if nargin<6,diagnosticOptions=struct;end\ndiagnostics=prepareDiagnostics(diagnosticOptions);\n','')
observed=replace_once(observed,
    '            refreshedThisIteration=isempty(cache)||mod(iteration-1,options.jacobianRefreshEvery)==0;\n'
    '            if refreshedThisIteration',
    '            if isempty(cache)||mod(iteration-1,options.jacobianRefreshEvery)==0')
for line in ["                        recordCandidate('accepted',newnorm,'','');\n",
             "                    recordCandidate('valid-not-decreasing',newnorm,'','');\n",
             "                    recordCandidate('state-failure',NaN,exception.identifier,exception.message);\n"]:
    observed=replace_once(observed,line,'')
start=observed.index('if diagnostics.enabled\n')
end=observed.index('    function [r,i]=evaluate(state)',start)
observed=observed[:start]+observed[end:]
observed=observed[:observed.index('\nfunction d=prepareDiagnostics(options)')]


def canonical(text):
    return '\n'.join(line for line in text.splitlines() if line.strip() and not line.lstrip().startswith('%'))


assert canonical(observed)==canonical(original),'观察代码之外的默认算法正文发生改变'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
report=dict(status='passed-observer-source-body-equals-default-after-listed-observation-removal',
    defaultSolverSHA256=sha(core),observerSHA256=sha(observer),testSHA256=sha(Path(__file__)),
    onlyListedObservationHooksRemoved=True,completeResultEqualityRequiresRuntimeTest=True)
if len(sys.argv)>1:
    out=Path(sys.argv[1]);assert not out.exists();out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
print(json.dumps(report))
