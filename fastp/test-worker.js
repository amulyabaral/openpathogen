self.onmessage = async (e) => {
  const post = t => self.postMessage(t);
  const Module = {
    noExitRuntime: true,
    locateFile: p => new URL(p, self.location.href).href,
    print: t => post('OUT ' + t),
    printErr: t => post('ERR ' + t),
    onRuntimeInitialized: () => post('READY'),
    onAbort: r => post('ABORT ' + r),
  };
  self.Module = Module;
  importScripts(new URL('fastp.js', self.location.href).href);
  const wait = setInterval(() => {
    if (!Module.FS) return;
    clearInterval(wait);
    const FS = Module.FS;
    post('FS OK');
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
    try { Module.callMain(args); post('EXIT ok'); }
    catch (err) { post('CAUGHT ' + (err && err.message || String(err))); }
    try { post('report: ' + FS.readFile('/out/report.json', {encoding:'utf8'}).length + ' bytes'); }
    catch (err) { post('no report'); }
  }, 300);
};
