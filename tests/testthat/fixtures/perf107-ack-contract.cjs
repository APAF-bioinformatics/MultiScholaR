const vm = require('node:vm');
const fs = require('node:fs');
const assert = require('node:assert/strict');
const source = fs.readFileSync(process.argv[2], 'utf8');
let now = 0, next = 0, callbacks = new Map(), sent = [], handler;
let image = {complete: false, naturalWidth: 0, naturalHeight: 0, alt: 'old', src: 'old'};
let text = {textContent: 'old'};
let busy = false;
const context = {
  window: {}, Date: {now: () => now},
  setTimeout: (fn, delay) => { callbacks.set(++next, {fn, due: now + delay}); return next; },
  clearTimeout: id => callbacks.delete(id),
  document: {getElementById: id => id === 'text' ? text : {
    querySelector: () => image, classList: {contains: () => busy}
  }}
};
function tick(ms) {
  for (let end = now + ms; now < end; now += 25) {
    for (const [id, entry] of [...callbacks]) if (entry.due <= now) {
      callbacks.delete(id); entry.fn();
    }
  }
}
vm.runInNewContext(source, context);
assert.equal(context.window.MultiScholaRPerf107Ack, undefined);
context.window.Shiny = {
  addCustomMessageHandler: (name, fn) => { handler = fn; },
  setInputValue: (name, value) => sent.push(value)
};
tick(100);
const message = {request_id: 'one', content_id: 'new', text_id: 'text', plot_id: 'plot',
  expected_text: 'accepted result', plot_identity: 'new', input_id: 'ack'};
handler(message);
for (let i = 0; i < 31; i++) handler(message);
assert.equal(callbacks.size, 1, 'duplicate delivery must not multiply polling loops');
image = {complete: true, naturalWidth: 100, naturalHeight: 100, alt: 'old', src: 'old'};
text.textContent = message.expected_text;
tick(400);
assert.equal(sent.length, 0, 'old loaded plot is not the accepted result');
image.alt = 'new'; image.src = 'new'; busy = true;
tick(400);
assert.equal(sent.length, 0, 'recalculating plot must not acknowledge');
busy = false;
tick(150);
assert.equal(sent.length, 0, 'stable interval is required');
tick(200);
assert.equal(sent.length, 1);
assert.equal(sent[0].request_id, 'one');
handler(message); tick(400);
assert.equal(sent.length, 1, 'duplicate messages must not duplicate acknowledgement');
handler({...message, request_id: 'two', content_id: 'two', plot_identity: 'two'});
handler({...message, request_id: 'three', content_id: 'three', plot_identity: 'three'});
assert.equal(callbacks.size, 1, 'replacement request cancels previous polling');
image.alt = 'three'; image.src = 'three'; tick(400);
assert.equal(sent.length, 2);
assert.equal(sent[1].request_id, 'three');
console.log('ack lifecycle contract passed');
