# -*- coding: utf-8 -*-
"""原字节读取已冻结的真实C/F预检，转换为MATLAB输入，不重新选择代表或结构。"""
from pathlib import Path
import json, hashlib, os
import numpy as np
from scipy.io import savemat
ROOT = Path(__file__).resolve().parents[1]
folder = ROOT/'outputs/v83_validation/spatialreductionpreflight01'
source = folder/'spatial_reduction_topology.npz'
summary = json.loads((folder/'spatial_reduction_preflight.json').read_text(encoding='utf8'))
target = folder/'spatial_reduction_input.mat'
assert not target.exists()
digest = hashlib.sha256(source.read_bytes()).hexdigest()
assert digest == summary['topologyEvidenceSHA256']
s = np.load(source)
data = {key:s[key].astype(np.float64) for key in s.files}
data['topologyEvidenceSHA256'] = digest
data['sourceEvidenceSHA256'] = summary['sourceBinding']['sourceFileSHA256']
savemat(target, data, do_compression=True)
print(json.dumps(dict(status='converted-frozen-spatial-CF-topology-to-MATLAB-without-reselection', processID=os.getpid(),
    inputSHA256=digest, outputSHA256=hashlib.sha256(target.read_bytes()).hexdigest(), newActualNSEvaluations=0)))
