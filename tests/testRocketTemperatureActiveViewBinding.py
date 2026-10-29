"""拒绝将旧基态、其他步的方向或被修改的来源贴到连续试算结果。"""
import json, sys, re
from pathlib import Path
import buildRocketTemperatureActiveWorkbench as workbench
if len(sys.argv) > 1:
    workbench.configure(sys.argv[1])
ROOT, OUT, STUDY, bind = workbench.ROOT, workbench.OUT, workbench.STUDY, workbench.bind
report_folder = Path(sys.argv[2]).resolve() if len(sys.argv)>2 else OUT
report_folder.relative_to(ROOT); report_folder.mkdir(parents=True, exist_ok=True)
from screenRocketFullFLeastSquares import sha

binding = json.loads((OUT / 'view_source_binding.json').read_text('utf8'))
chosen = binding['selected']
study = json.loads(STUDY.read_text('utf8'))
steps = study['steps'] if isinstance(study['steps'], list) else [study['steps']]
step = next(x for x in steps if x['step'] == chosen['step'])
state = ROOT / 'outputs/v83_validation' / chosen['displayCase'] / 'wall_velocity_krylov_evidence.mat'
candidate = step['directionFile']
valid = bind(STUDY, candidate, state, chosen['step'])
assert not valid['oldActiveSetInherited']
rejected = []
try:
    wrong_state = study['initialStateFile']
    if chosen.get('displayedStateUnchanged'):
        number = int(re.fullmatch(r'krylov(\d+)', chosen['displayCase']).group(1))
        wrong_state = ROOT/'outputs/v83_validation'/f'krylov{number-1}'/'wall_velocity_krylov_evidence.mat'
    bind(STUDY, candidate, wrong_state, chosen['step'])
except AssertionError:
    rejected.append('old-initial-base-as-current-result')
if chosen['step'] > 1:
    try:
        bind(STUDY, steps[0]['directionFile'], state, chosen['step'])
    except AssertionError:
        rejected.append('different-current-step-direction')
bad = report_folder / 'bad_active_source.json'; assert not bad.exists()
for item in steps:
    if item['step'] == chosen['step']:
        item['sourceBinding']['sourceFileSHA256'] = '0' * 64
bad.write_text(json.dumps(study), encoding='utf8')
try:
    bind(bad, candidate, state, chosen['step'])
except AssertionError:
    rejected.append('wrong-current-source-SHA')
assert len(rejected) == (3 if chosen['step'] > 1 else 2)
out = dict(status='passed-old-base-other-step-and-current-source-SHA-rejections', rejected=rejected,
           studySHA256=sha(STUDY), candidateSHA256=sha(candidate), scriptSHA256=sha(__file__),
           physicalSolution=False, newProductionSolvedAreaM2=0)
target = report_folder / 'view_binding_rejections.json'; assert not target.exists()
target.write_text(json.dumps(out, indent=2) + '\n', encoding='utf8')
print(json.dumps(out))
