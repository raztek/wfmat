const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch();
  const p = await b.newPage();
  await p.goto('file://' + process.cwd() + '/spec.html', { waitUntil: 'networkidle' });
  await p.waitForFunction(() => document.querySelectorAll('.katex').length > 0);
  await p.evaluate(() => document.fonts.ready);
  const errs = await p.evaluate(() => document.querySelectorAll('.katex-error').length);
  console.log('katex nodes', await p.evaluate(() => document.querySelectorAll('.katex').length), 'errors', errs);
  await p.pdf({ path: 'spec.pdf', format: 'Letter', printBackground: true,
    margin: { top: '20mm', bottom: '18mm', left: '19mm', right: '19mm' },
    displayHeaderFooter: true,
    headerTemplate: '<div style="font-size:7pt;color:#888;width:100%;padding:0 19mm;">Wavefront MAT: Design and Specification · v' + (process.env.SPEC_VERSION || '') + '</div>',
    footerTemplate: '<div style="font-size:7pt;color:#888;width:100%;text-align:right;padding:0 19mm;">Page <span class="pageNumber"></span> of <span class="totalPages"></span></div>' });
  await b.close();
})();
