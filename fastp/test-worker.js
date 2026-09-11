/* test-worker.js — smoke test for the wasm64 + pthreads fastp build.
 * Served from a cross-origin-isolated server (benchmark/harness/server.mjs):
 * threads need SharedArrayBuffer. Driven by fastp/test.html?flags=...; appends
 * the URL flags to the base command. Posts OUT/ERR lines, the exit code, and
 * the report.json size. */
self.onmessage = async (e) => {
  const post = t => self.postMessage(t);
  const JS_URL = new URL('fastp.js', self.location.href).href;
  const Module = {
    locateFile: p => new URL(p, JS_URL).href,
    // pthread workers must spawn from the glue itself, not this runner script
    mainScriptUrlOrBlob: JS_URL,
    print: t => post('OUT ' + t),
    printErr: t => post('ERR ' + t),
    onAbort: r => post('ABORT ' + r),
  };
  self.Module = Module;
  importScripts(JS_URL);
  const wait = setInterval(() => {
    if (!Module.FS) return;
    clearInterval(wait);
    const FS = Module.FS;
    post('FS OK (crossOriginIsolated=' + self.crossOriginIsolated + ')');
    // tiny built-in input
    let seq = '';
    for (let i = 0; i < 3000; i++) {
      const s = Array.from({length: 150}, () => 'ACGT'[i % 4]).join('');
      seq += `@r${i}/1\n${s}\n+\n${'I'.repeat(150)}\n`;
    }
    let seq2 = seq.replaceAll('/1', '/2');
    FS.mkdir('/in'); FS.mkdir('/out');
    FS.writeFile('/in/r1.fastq', seq);
    FS.writeFile('/in/r2.fastq', seq2);
    const args = ['--in1','/in/r1.fastq','--in2','/in/r2.fastq',
      '--out1','/out/o1.fq','--out2','/out/o2.fq','--json','/out/report.json','--thread','1',
      ...e.data.flags];
    post('$ fastp ' + args.join(' '));
    Module.ccall('fastp_run', 'number', ['string'], [args.join(' ')]);
    const poll = setInterval(() => {
      const state = (() => { try { return Module.ccall('fastp_state', 'number'); } catch { return 2; } })();
      if (state !== 2) return;
      clearInterval(poll);
      const code = (() => { try { return Module.ccall('fastp_code', 'number'); } catch { return '?'; } })();
      post('EXIT ' + code);
      try { post('report: ' + FS.readFile('/out/report.json', {encoding:'utf8'}).length + ' bytes'); }
      catch (err) { post('no report'); }
    }, 50);
  }, 300);
};
