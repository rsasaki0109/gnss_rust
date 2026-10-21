"""Measure one Linux solver subprocess without mixing prior children into peak RSS."""
import json
import os
from pathlib import Path
import sys
import time


def main():
    record=Path(sys.argv[1]);command=sys.argv[2:]
    assert command and sys.platform.startswith('linux')
    start=time.perf_counter();pid=os.fork()
    if pid==0:
        try:os.execv(command[0],command)
        except OSError as error:
            os.write(2,(str(error)+'\n').encode());os._exit(127)
    _,status,usage=os.wait4(pid,0);code=os.waitstatus_to_exitcode(status)
    record.write_text(json.dumps(dict(command=command,returncode=code,wall_seconds=time.perf_counter()-start,user_seconds=usage.ru_utime,system_seconds=usage.ru_stime,peak_RSS_KiB=usage.ru_maxrss,scope='one process, includes input parsing and CSV output',performance_accepted=False),indent=2)+'\n')
    raise SystemExit(code if code>=0 else 128-code)


if __name__=='__main__':main()
