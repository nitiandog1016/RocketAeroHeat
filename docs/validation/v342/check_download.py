"""用小文件覆盖安装包续传、忽略 Range、截断响应和错误区间的回归路径。"""
from hashlib import sha256
from io import BytesIO
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
import installGmshRuntime as runtime

DATA = b"published-wheel-bytes"
OFFSET = 7


def check(kind):
    calls = []
    def response(req, timeout):
        calls.append(req.get_header("Range"))
        assert timeout > 0 and calls[-1] == f"bytes={OFFSET}-{len(DATA)-1}"
        data = DATA[OFFSET:]
        if kind == "full":
            data = DATA
        elif kind == "partial" and len(calls) == 1:
            data = data[:-2]
        stream = BytesIO(data)
        stream.status = 200 if kind == "full" else 206
        start = 0 if kind == "wrong-range" else OFFSET
        stream.headers = {"Content-Range": f"bytes {start}-{len(DATA)-1}/{len(DATA)}"}
        return stream

    with TemporaryDirectory() as name, \
            patch.object(runtime.urllib.request, "urlopen", response):
        out = Path(name) / "official.whl"
        partial = out.with_name("download.part")
        partial.write_bytes(DATA[:OFFSET])
        if kind == "wrong-range":
            try:
                runtime.download_pinned(runtime.URL, out, partial, len(DATA), sha256(DATA).hexdigest())
            except ValueError:
                assert not out.exists() and partial.read_bytes() == DATA[:OFFSET]
            else:
                raise AssertionError("错误字节区间被接收")
        else:
            runtime.download_pinned(runtime.URL, out, partial, len(DATA), sha256(DATA).hexdigest())
            assert out.read_bytes() == DATA and not partial.exists()
        assert len(calls) == (2 if kind == "partial" else 1)


if __name__ == "__main__":
    for mode in ("range", "full", "partial", "wrong-range"):
        check(mode)
    print("Verified exact resume, full response, truncated retry and wrong-range rejection.")
