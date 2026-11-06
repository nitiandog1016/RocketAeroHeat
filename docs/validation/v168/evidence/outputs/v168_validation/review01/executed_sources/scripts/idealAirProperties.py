"""固定理想空气批处理入口，复用同一C++核并拒绝部分成功的状态数组。"""
import ctypes
from pathlib import Path

import numpy as np

from su2ThermalRuntime import read, sha

ROOT = Path(__file__).resolve().parents[1]
RUNTIME = ROOT/'tools/runtime/ideal_air/batch-v1'
FIELDS = ('T', 'rho', 'p', 'e', 'h', 'cp', 'sound')


class IdealAir:
    """模块一。每个准备任务绑定运行库与源码，调用期间无共享可变热状态。"""

    def __init__(self):
        manifest_path = RUNTIME/'runtime.json'
        manifest = read(manifest_path)
        assert manifest['model'] == 'Lemmon2000-ideal-Air-CoolProp7.2.0-OTH-v1'
        assert manifest['abi'] == 'PT-Ph-rhoe-7fields-v1'
        library = RUNTIME/manifest['library']
        assert library.parent == RUNTIME and sha(library) == manifest['sha256']
        for path, expected in manifest['sourceSHA256'].items():
            assert sha(Path(path)) == expected, path
        self.locks = dict(manifest['sourceSHA256'])
        self.locks.update({str(manifest_path): sha(manifest_path), str(library): sha(library),
                           str(Path(__file__).resolve()): sha(Path(__file__).resolve())})
        self.library = ctypes.CDLL(str(library))
        self.function = self.library.rah_ideal_air_batch
        array = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
        self.function.argtypes = [ctypes.c_int, ctypes.c_size_t, array, array, array,
                                  ctypes.c_char_p, ctypes.c_size_t]
        self.function.restype = ctypes.c_size_t

    def state(self, mode, first, second):
        """模块二。仅接受同尺寸或标量扩展，返回一维具名列，非法批次整体拒绝。"""
        modes = {'PT': 0, 'Ph': 1, 'rhoe': 2}
        if mode not in modes:
            raise ValueError('Unsupported ideal-air state pair')
        a = np.asarray(first, dtype=np.float64)
        b = np.asarray(second, dtype=np.float64)
        if a.ndim > 1 or b.ndim > 1 or (a.size != 1 and b.size != 1 and a.shape != b.shape):
            raise ValueError('State inputs must be scalar or equal one-dimensional arrays')
        a, b = np.broadcast_arrays(a, b)
        a = np.ascontiguousarray(a.reshape(-1))
        b = np.ascontiguousarray(b.reshape(-1))
        if not len(a) or not np.isfinite(a).all() or not np.isfinite(b).all():
            raise ValueError('Empty or nonfinite state input')
        output = np.empty((len(a), len(FIELDS)), dtype=np.float64)
        message = ctypes.create_string_buffer(512)
        failed = self.function(modes[mode], len(a), a, b, output, message, len(message))
        if failed:
            raise ValueError('Ideal-air row %d: %s' % (failed-1, message.value.decode('utf8')))
        if not np.isfinite(output).all():
            raise ValueError('Nonfinite ideal-air output')
        return {name: output[:, index] for index, name in enumerate(FIELDS)}
