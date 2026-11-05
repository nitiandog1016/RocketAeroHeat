import subprocess,sys,json,os
from pathlib import Path
p=subprocess.Popen([sys.executable,"-c","import time; time.sleep(20)"],creationflags=subprocess.CREATE_NO_WINDOW)
Path("owned_pids.json").write_text(json.dumps([os.getpid(),p.pid]))
p.wait()
