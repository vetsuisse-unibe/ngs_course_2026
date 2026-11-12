// NGS Course 2026 site: callout boxes, contents box, copy buttons.
// Works on the HTML that GitHub Pages (kramdown) makes from the Markdown,
// so the exercise files stay plain Markdown.
(function () {
  var content = document.getElementById("content");
  if (!content) return;
  var isExercise = document.body.classList.contains("exercise");

  function startsWith(el, re) {
    return el && re.test((el.textContent || "").trim());
  }
  function isList(el) {
    return el && (el.tagName === "OL" || el.tagName === "UL");
  }
  function isCode(el) {
    return el && el.classList && el.classList.contains("highlighter-rouge");
  }
  function makeCallout(kind, title, first) {
    var box = document.createElement("div");
    box.className = "callout callout-" + kind;
    var t = document.createElement("div");
    t.className = "callout-title";
    t.textContent = title;
    box.appendChild(t);
    first.parentNode.insertBefore(box, first);
    return box;
  }
  // remove "Questions:" etc. from the start of a paragraph; drop it if nothing is left
  function stripLabel(p, re) {
    var node = p.firstChild;
    while (node && node.nodeType === 3 && !node.textContent.trim()) node = node.nextSibling;
    if (node && node.nodeType === 3) node.textContent = node.textContent.replace(re, "");
    else if (node && node.nodeType === 1 && re.test(node.textContent)) p.removeChild(node);
    if (!p.textContent.trim()) { p.remove(); return null; }
    return p;
  }

  if (isExercise) {
    var blocks = Array.prototype.slice.call(content.children);
    blocks.forEach(function (el) {
      if (!el.parentNode || el.parentNode !== content) return;

      // "Question:" / "Questions:" paragraph + following list(s)
      if (el.tagName === "P" && startsWith(el, /^Questions?:/)) {
        var label = /^Questions:/.test(el.textContent.trim()) ? "Questions" : "Question";
        var box = makeCallout("questions", label, el);
        var next = el.nextElementSibling;
        var p = stripLabel(el, /^\s*Questions?:\s*/);
        if (p) box.appendChild(p);
        while (next && (isList(next) || (isCode(next) && isList(next.nextElementSibling)))) {
          var after = next.nextElementSibling;
          box.appendChild(next);
          next = after;
        }
        return;
      }

      // "#### To remember" + list
      if (/^H[3-6]$/.test(el.tagName) && startsWith(el, /^To remember/)) {
        var rbox = makeCallout("remember", "Remember", el);
        var rnext = el.nextElementSibling;
        el.remove();
        while (isList(rnext)) { var ra = rnext.nextElementSibling; rbox.appendChild(rnext); rnext = ra; }
        return;
      }

      // "#### Please ask for help ..."
      if (/^H[3-6]$/.test(el.tagName) && startsWith(el, /^Please ask/)) {
        var ibox = makeCallout("info", "Help", el);
        var ip = document.createElement("p");
        ip.textContent = el.textContent;
        el.remove();
        ibox.appendChild(ip);
        return;
      }

      // "**For the brave.** ..." paragraph + following code
      if (el.tagName === "P" && el.firstElementChild && el.firstElementChild.tagName === "STRONG" &&
          /^For the brave/.test(el.firstElementChild.textContent) && el.textContent.trim().indexOf("For the brave") === 0) {
        var bbox = makeCallout("brave", "For the brave", el);
        var bnext = el.nextElementSibling;
        var bp = stripLabel(el, /^\s*For the brave[.:]?\s*/);
        if (bp) bbox.appendChild(bp);
        while (isCode(bnext)) { var ba = bnext.nextElementSibling; bbox.appendChild(bnext); bnext = ba; }
      }
    });

    // contents box from the h2 headings
    var heads = content.querySelectorAll("h2[id]");
    if (heads.length > 2) {
      var nav = document.createElement("nav");
      nav.className = "toc";
      nav.innerHTML = '<div class="toc-title">On this page</div><ul></ul>';
      var ul = nav.querySelector("ul");
      var links = {};
      heads.forEach(function (h) {
        var li = document.createElement("li");
        var a = document.createElement("a");
        a.href = "#" + h.id;
        a.textContent = h.textContent;
        li.appendChild(a);
        ul.appendChild(li);
        links[h.id] = a;
      });
      content.insertBefore(nav, content.firstChild);

      if ("IntersectionObserver" in window) {
        var current = null;
        var obs = new IntersectionObserver(function (entries) {
          entries.forEach(function (e) {
            if (e.isIntersecting && links[e.target.id]) {
              if (current) current.classList.remove("active");
              current = links[e.target.id];
              current.classList.add("active");
            }
          });
        }, { rootMargin: "0px 0px -70% 0px" });
        heads.forEach(function (h) { obs.observe(h); });
      }
    }
  }

  // wide tables scroll sideways on phones
  content.querySelectorAll("table").forEach(function (t) {
    var w = document.createElement("div");
    w.className = "table-wrap";
    t.parentNode.insertBefore(w, t);
    w.appendChild(t);
  });

  // copy buttons on every code block
  content.querySelectorAll("div.highlighter-rouge").forEach(function (block) {
    var pre = block.querySelector("pre");
    if (!pre) return;
    var btn = document.createElement("button");
    btn.type = "button";
    btn.className = "copy-btn";
    btn.textContent = "Copy";
    btn.addEventListener("click", function () {
      var text = pre.innerText.replace(/\n$/, "");
      var done = function () {
        btn.textContent = "Copied";
        btn.classList.add("copied");
        setTimeout(function () { btn.textContent = "Copy"; btn.classList.remove("copied"); }, 1500);
      };
      if (navigator.clipboard && window.isSecureContext) {
        navigator.clipboard.writeText(text).then(done, function () {});
      } else {
        var ta = document.createElement("textarea");
        ta.value = text;
        document.body.appendChild(ta);
        ta.select();
        try { document.execCommand("copy"); done(); } catch (e) {}
        document.body.removeChild(ta);
      }
    });
    block.appendChild(btn);
  });
})();
