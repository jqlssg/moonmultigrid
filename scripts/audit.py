"""Reproducible source-size and meaningful-history ledger (no Git writes)."""
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parents[1]
files = sorted(p for p in root.glob('*.mbt') if not p.stem.endswith(('_test','_wbtest')))
physical = 0
code = 0
for p in files:
    lines = p.read_text(encoding='utf-8-sig').splitlines()
    n = sum(bool(s.strip()) and not s.lstrip().startswith('//') for s in lines)
    physical += len(lines)
    code += n
    print('%4d physical / %4d nonblank noncomment  %s' % (len(lines),n,p.name))
print('TOTAL production: %d physical; %d nonblank noncomment; excludes tests, examples, generated interfaces' % (physical,code))
assert code >= 1500
history = subprocess.check_output(['git','log','--reverse','--format=%h %s'],cwd=str(root),universal_newlines=True).splitlines()
print('History (%d nonempty commits; substantive verification requires reviewing changes):' % len(history))
for entry in history:
    sha = entry.split()[0]
    names = subprocess.check_output(['git','diff-tree','--root','--no-commit-id','--name-only','-r',sha],cwd=str(root),universal_newlines=True).splitlines()
    assert names, 'empty commit'
    print(entry + ' | ' + ', '.join(names))
assert len(history) >= 20
