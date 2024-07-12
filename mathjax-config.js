// Copied from https://docs.mathjax.org/en/latest/upgrading/v2.html#changes-in-the-mathjax-api
window.MathJax = {
  loader: {load: ['[tex]/physics']},
  tex: {
    packages: {'[+]': ['physics']},
    inlineMath: [['\\(', '\\)'], ['$', '$']],
    // displayMath: [['\\[', '\\]'], ['$$', '$$']],
    macros: {
      coloneqq: "\\mathrel{:=}",
      eqqcolon: "\\mathrel{=:}",
      ii: "\\mathrm{i}",
      ee: "\\mathrm{e}",
      Span: "\\operatorname{span}",
    },
  },
  processEscapes: true,
  options: {
    renderActions: {
      findScript: [10, function (doc) {
        for (const node of document.querySelectorAll('script[type^="math/tex"]')) {
          const display = !!node.type.match(/; *mode=display/);
          const math = new doc.options.MathItem(node.textContent, doc.inputJax[0], display);
          const text = document.createTextNode('');
          node.parentNode.replaceChild(text, node);
          math.start = {node: text, delim: '', n: 0};
          math.end = {node: text, delim: '', n: 0};
          doc.math.push(math);
        }
      }, '']
    }
  }
};
