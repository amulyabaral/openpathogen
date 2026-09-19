# wasm32 has a 32-bit `long`; fastp keeps read/base totals in `long`, which
# overflows past 2^31 bases. Rewrite every `long` type (outside string
# literals and comments) as `long long`, which is 64-bit everywhere.
import re, sys, pathlib
WORD = re.compile(r'\blong(?:\s+long)?\b')
def fix_code(code):
    return WORD.sub('long long', code)
def fix_line(line, state):
    out, i, n = [], 0, len(line)
    buf = ''
    while i < n:
        if state['block']:
            j = line.find('*/', i)
            if j < 0: out.append(line[i:]); return ''.join(out)
            out.append(line[i:j+2]); i = j + 2; state['block'] = False; continue
        c = line[i]
        if line.startswith('//', i):
            out.append(fix_code(buf)); buf = ''; out.append(line[i:]); return ''.join(out)
        if line.startswith('/*', i):
            out.append(fix_code(buf)); buf = ''; state['block'] = True; continue
        if c in '"\'':
            out.append(fix_code(buf)); buf = ''
            j = i + 1
            while j < n and line[j] != c:
                j += 2 if line[j] == '\\' else 1
            out.append(line[i:j+1]); i = j + 1; continue
        buf += c; i += 1
    out.append(fix_code(buf))
    return ''.join(out)
for p in sys.argv[1:]:
    path = pathlib.Path(p)
    state = {'block': False}
    lines = path.read_text().splitlines(keepends=True)
    path.write_text(''.join(fix_line(l, state) for l in lines))
