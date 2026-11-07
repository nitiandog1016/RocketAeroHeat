"""原生完成证据独立于启动器exit0：终止标志、历史及最终载荷必须相容。"""
from pathlib import Path
import hashlib,json,re
import numpy as np


def completion_log_for_receipt(case,receipt_path,receipt):
    """直接回执读取自身日志；v178首次封存别名必须以原回执SHA找回准确来源。"""
    candidate=Path(receipt_path).parent/'process00/stdout.txt'
    if candidate.is_file():return candidate
    expected=receipt.get('originalRunReceiptSHA256')
    original=Path(case).parent/(Path(case).name+'_run')/'report.json'
    if expected and original.is_file() and hashlib.sha256(original.read_bytes()).hexdigest()==expected:
        source=json.loads(original.read_text())
        if source.get('case')==receipt.get('case') and source.get('commands')==receipt.get('commands'):
            return original.parent/'process00/stdout.txt'
    raise ValueError('missing-source-bound-native-completion-log')


def audit_completion(case,stdout):
    case=Path(case);text=Path(stdout).read_text(errors='replace')
    manifest=json.loads((case/'case_manifest.json').read_text())
    if 'explicit' in manifest.get('condition',{}):
        from auditEilmerExplicitCompletion import audit_explicit_completion
        return audit_explicit_completion(case,stdout,manifest)
    fatal=re.findall(r'^.*(?:OutOfMemoryError|Memory allocation failed|core\.exception\.).*$',text,re.M)
    if fatal:raise ValueError('native-fatal-error: '+fatal[0])
    reasons=re.findall(r'^STOP-REASON:\s*(.+)$',text,re.M)
    steps=re.findall(r'^FINAL-STEP:\s*(\d+)\s*$',text,re.M)
    if len(reasons)!=1 or len(steps)!=1:raise ValueError('missing-or-ambiguous-native-completion-markers')
    reason=reasons[0].strip();step=int(steps[0])
    if reason not in ('maximum-steps','maximum-wall-clock','absolute-global-residual-target','relative-global-residual-target'):
        raise ValueError('unaccepted-native-stop: '+reason)
    h=np.atleast_1d(np.genfromtxt(case/'lmrsim/diagnostics/nk-diagnostics',names=True,deletechars=''))
    if not len(h) or any(not np.isfinite(h[k]).all() for k in h.dtype.names):raise ValueError('missing-or-nonfinite-history')
    if step<=0 or h['step'][-1]!=step:raise ValueError('final-history-step-mismatch')
    meta=np.atleast_2d(np.loadtxt(case/'lmrsim/loads/loads-metadata',skiprows=1))
    if meta[-1,2]!=step:raise ValueError('final-loads-step-mismatch')
    m=json.loads((case/'case_manifest.json').read_text())
    if m.get('condition',{}).get('development',{}).get('profile')=='rah-nk-real-v1':
        # 实数新合同还需实际二进制自报；仅有配置名称不足以证明运行分派。
        if re.findall(r'^Number-type:\s*(\S+)\s*$',text,re.M)!=['real']:
            raise ValueError('wrong-native-number-type-for-real-profile')
    expected=m['mesh']['wallFaces'] if 'mesh' in m else m['wallFaces']
    files=list((case/'lmrsim/loads'/('%04d'%int(meta[-1,0]))).glob('*-wall.dat'))
    count=0
    for f in files:
        a=np.atleast_1d(np.genfromtxt(f,names=True,deletechars=''))
        if any(not np.isfinite(a[k]).all() for k in a.dtype.names):raise ValueError('nonfinite-final-loads')
        count+=len(a)
    if count!=expected:raise ValueError('incomplete-final-wall-coverage')
    return dict(passed=True,stopReason=reason,finalStep=step,wallFaces=count,physicsQualified=False,
                nativeLogSHA256=hashlib.sha256(Path(stdout).read_bytes()).hexdigest(),
                scope='Native completion evidence; does not establish convergence, gradient validity or physical qualification')
