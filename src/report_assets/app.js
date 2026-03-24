(function () {
  const reportElement = document.getElementById("seqcheck-report-data");
  const libElement = document.getElementById("seqspec-lib-data");
  const report = JSON.parse(reportElement.textContent);
  const libData = libElement ? JSON.parse(libElement.textContent) : null;
  const meta = report.meta || {};
  const results = report.results || [];
  const app = document.getElementById("app");
  const SEV = { error: 0, warning: 1, interpretation: 2, pass: 3 };
  const repositoryUrl = window.SEQCHECK_REPOSITORY || "https://github.com/sbooeshaghi/seqcheck";
  const seqcheckVersion = window.SEQCHECK_VERSION || "";
  const reportInputPath = window.SEQCHECK_REPORT_INPUT || "report.json";
  const reportOutputPath = window.SEQCHECK_REPORT_OUTPUT || "report.html";

  let filterText = "";
  let filterSev = null;

  function esc(value) {
    if (value == null) {
      return '<span class="val-null">\u2014</span>';
    }
    const div = document.createElement("div");
    div.textContent = String(value);
    return div.innerHTML;
  }

  function fmtVal(value) {
    if (value == null) {
      return '<span class="val-null">\u2014</span>';
    }
    if (typeof value === "boolean") {
      return value ? "true" : "false";
    }
    if (typeof value === "number") {
      return Number.isInteger(value) ? value.toLocaleString() : value.toPrecision(4);
    }
    return esc(value);
  }

  function worstSev(assessments) {
    let worst = "pass";
    for (const assessment of assessments || []) {
      const candidate = assessment.type || "pass";
      if ((SEV[candidate] || 3) < (SEV[worst] || 3)) {
        worst = candidate;
      }
    }
    return worst;
  }

  function sevLabel(severity) {
    return severity === "interpretation" ? "info" : severity;
  }

  function naturalDesc(result) {
    const assessments = (result.assessment || []).filter((assessment) => assessment.type !== "pass");
    if (assessments.length) {
      assessments.sort((left, right) => (SEV[left.type] || 3) - (SEV[right.type] || 3));
      return assessments[0].description;
    }
    if (result.assessment && result.assessment.length) {
      return result.assessment[0].description;
    }
    return result.check + " completed.";
  }

  function contextLine(result) {
    const parts = [result.check];
    if (result.files && result.files.length) {
      parts.push(result.files.join(", "));
    }
    if (result.regions && result.regions.length) {
      parts.push(result.regions.join(", "));
    }
    return parts.join(" \u00b7 ");
  }

  function findMetric(items, name) {
    return (items || []).find((item) => item.name === name) || null;
  }

  function inputCheckResult() {
    return results.find((result) => result.check === "input_check") || null;
  }

  function suppliedInputPaths() {
    const inputCheck = inputCheckResult();
    const supplied = findMetric(inputCheck ? inputCheck.observed : [], "supplied_input_paths");
    const value = supplied && supplied.data ? supplied.data.value : [];
    if (!Array.isArray(value)) {
      return [];
    }
    return value
      .map((row) => row && row.input_path)
      .filter((path) => typeof path === "string" && path.length > 0);
  }

  function buildInvocation() {
    const args = ["seqcheck", meta.command || "check", "--format", "json"];
    if (meta.spec) {
      args.push("-s", meta.spec);
    }
    if (meta.modality) {
      args.push("-m", meta.modality);
    }
    if (typeof meta.requested_reads === "number") {
      args.push("-n", String(meta.requested_reads));
    }
    suppliedInputPaths().forEach((path) => args.push(path));
    return args.join(" ");
  }

  function buildReportInvocation() {
    return ["seqcheck", "report", "-i", reportInputPath, "-o", reportOutputPath].join(" ");
  }

  function renderItem(item) {
    const data = item.data;
    if (!data) {
      return "";
    }
    if (data.kind === "scalar") {
      const unit = data.unit ? `<span class="sc-unit">${esc(data.unit)}</span>` : "";
      return `<div class="sc-row"><span class="sc-name">${esc(item.name)}</span><span class="sc-val">${fmtVal(data.value)}${unit}</span></div>`;
    }
    if (data.kind === "series") {
      if (!data.value || !data.value.length) {
        return `<div class="sc-row"><span class="sc-name">${esc(item.name)}</span><span class="sc-val val-null">empty</span></div>`;
      }
      const uid = "t" + Math.random().toString(36).slice(2, 8);
      let html = `<div class="tbl-toggle" data-tbl="${uid}">${esc(item.name)} \u2014 ${data.value.length} entries \u25b8</div>`;
      html += `<div class="tbl-inner" id="${uid}"><div class="tbl-wrap"><table class="dtbl"><thead><tr><th>key</th><th>value${data.unit ? " (" + esc(data.unit) + ")" : ""}</th></tr></thead><tbody>`;
      data.value.forEach((row) => {
        html += `<tr><td>${fmtVal(row.key)}</td><td>${fmtVal(row.value)}</td></tr>`;
      });
      html += "</tbody></table></div></div>";
      return html;
    }
    if (data.kind === "records") {
      if (!data.value || !data.value.length) {
        return `<div class="sc-row"><span class="sc-name">${esc(item.name)}</span><span class="sc-val val-null">empty</span></div>`;
      }
      const keys = [];
      const seen = new Set();
      data.value.forEach((record) => {
        Object.keys(record).forEach((key) => {
          if (!seen.has(key)) {
            seen.add(key);
            keys.push(key);
          }
        });
      });
      const uid = "t" + Math.random().toString(36).slice(2, 8);
      let html = `<div class="tbl-toggle" data-tbl="${uid}">${esc(item.name)} \u2014 ${data.value.length} records \u25b8</div>`;
      html += `<div class="tbl-inner" id="${uid}"><div class="tbl-wrap"><table class="dtbl"><thead><tr>${keys.map((key) => "<th>" + esc(key) + "</th>").join("")}</tr></thead><tbody>`;
      data.value.forEach((record) => {
        html += "<tr>" + keys.map((key) => "<td>" + fmtVal(record[key]) + "</td>").join("") + "</tr>";
      });
      html += "</tbody></table></div></div>";
      return html;
    }
    return "";
  }

  function buildMolSvg() {
    if (!libData || !libData.regions || !libData.regions.length) {
      return "";
    }

    const regions = libData.regions;
    const reads = libData.reads || [];
    const totalBp = libData.total_bp || regions.reduce((sum, region) => sum + region.len, 0);
    const pad = { left: 10, right: 10 };
    const molW = 860;
    const barY = 70;
    const barH = 22;
    const svgW = molW + pad.left + pad.right;
    const bpScale = molW / Math.max(totalBp, 1);
    const minPx = 18;
    const rects = [];

    let rawWidths = regions.map((region) => region.len * bpScale);
    let deficit = 0;
    rawWidths = rawWidths.map((px) => {
      if (px < minPx) {
        deficit += minPx - px;
        return minPx;
      }
      return px;
    });
    if (deficit > 0) {
      const shrinkTotal = rawWidths.filter((px) => px > minPx).reduce((sum, px) => sum + px, 0);
      if (shrinkTotal > 0) {
        rawWidths = rawWidths.map((px) => (px > minPx ? px - deficit * (px / shrinkTotal) : px));
      }
    }

    let currentX = pad.left;
    regions.forEach((region, index) => {
      const width = rawWidths[index];
      rects.push({ ...region, x: currentX, w: width });
      currentX += width;
    });

    function bpToX(bp) {
      for (const region of rects) {
        const regionStop = region.bp_start + region.len;
        if (bp <= regionStop) {
          const fraction = region.len === 0 ? 0 : (bp - region.bp_start) / region.len;
          return region.x + fraction * region.w;
        }
      }
      const last = rects[rects.length - 1];
      return last.x + last.w;
    }

    function seqTypeColor(sequenceType) {
      if (sequenceType === "onlist") {
        return "var(--reg-onlist)";
      }
      if (sequenceType === "random") {
        return "var(--reg-random)";
      }
      return "var(--reg-fixed)";
    }

    function seqTypeStroke(sequenceType) {
      if (sequenceType === "onlist") {
        return "var(--reg-onlist-stroke)";
      }
      if (sequenceType === "random") {
        return "var(--reg-random-stroke)";
      }
      return "var(--reg-fixed-stroke)";
    }

    let svg = "";
    svg += `<line x1="${pad.left}" y1="${barY - 1}" x2="${currentX}" y2="${barY - 1}" stroke="var(--border)" stroke-width="0.5"/>`;
    svg += `<line x1="${pad.left}" y1="${barY + barH + 1}" x2="${currentX}" y2="${barY + barH + 1}" stroke="var(--border)" stroke-width="0.5"/>`;

    rects.forEach((region) => {
      svg += `<rect class="region-rect" data-region="${esc(region.region_id)}" x="${region.x}" y="${barY}" width="${region.w}" height="${barH}" rx="2" fill="${seqTypeColor(region.sequence_type)}" stroke="${seqTypeStroke(region.sequence_type)}" stroke-width="0.5"/>`;
    });

    const labelY = barY + barH + 10;
    rects.forEach((region) => {
      const centerX = region.x + region.w / 2;
      if (region.w > 14) {
        svg += `<text class="region-label" x="${centerX}" y="${labelY}" text-anchor="end" transform="rotate(-40 ${centerX} ${labelY})">${esc(region.name)}</text>`;
      }
    });

    svg += `<text class="bp-label" x="${pad.left}" y="${barY - 3}" text-anchor="middle">0</text>`;
    rects.forEach((region) => {
      svg += `<text class="bp-label" x="${region.x + region.w}" y="${barY - 3}" text-anchor="middle">${region.bp_start + region.len}</text>`;
    });

    const readColors = ["#6366f1", "#059669", "#d97706", "#dc2626", "#7c3aed"];
    const posReads = reads.filter((read) => read.strand === "pos");
    const negReads = reads.filter((read) => read.strand === "neg");

    function drawRead(read, yBase, above, colorIndex) {
      const color = readColors[colorIndex % readColors.length];
      const x1 = bpToX(read.start);
      const x2 = bpToX(read.end);
      const arrowSize = 5;

      if (above) {
        svg += `<line x1="${x1}" y1="${yBase}" x2="${x2 - arrowSize}" y2="${yBase}" stroke="${color}" class="read-line"/>`;
        svg += `<polygon points="${x2},${yBase} ${x2 - arrowSize},${yBase - arrowSize} ${x2 - arrowSize},${yBase + arrowSize}" fill="${color}"/>`;
        svg += `<line x1="${x1}" y1="${yBase}" x2="${x1}" y2="${barY}" stroke="${color}" stroke-width="1" stroke-dasharray="2,2" opacity="0.4"/>`;
        svg += `<text class="read-label" x="${x1 + 3}" y="${yBase - 5}" fill="${color}">${esc(read.label)}</text>`;
      } else {
        svg += `<line x1="${x2}" y1="${yBase}" x2="${x1 + arrowSize}" y2="${yBase}" stroke="${color}" class="read-line"/>`;
        svg += `<polygon points="${x1},${yBase} ${x1 + arrowSize},${yBase - arrowSize} ${x1 + arrowSize},${yBase + arrowSize}" fill="${color}"/>`;
        svg += `<line x1="${x2}" y1="${yBase}" x2="${x2}" y2="${barY + barH}" stroke="${color}" stroke-width="1" stroke-dasharray="2,2" opacity="0.4"/>`;
        svg += `<text class="read-label" x="${x2 - 3}" y="${yBase + 13}" text-anchor="end" fill="${color}">${esc(read.label)}</text>`;
      }
    }

    let posY = barY - 14;
    posReads.forEach((read, index) => {
      drawRead(read, posY, true, index);
      posY -= 22;
    });

    let negY = barY + barH + 48;
    negReads.forEach((read, index) => {
      drawRead(read, negY, false, posReads.length + index);
      negY += 22;
    });

    const svgH = negY + 10;
    return `<svg class="mol-svg" viewBox="0 0 ${svgW} ${svgH}" width="100%" xmlns="http://www.w3.org/2000/svg">${svg}</svg>`;
  }

  function matches(result) {
    const severity = worstSev(result.assessment || []);
    if (filterSev && severity !== filterSev) {
      return false;
    }
    if (filterText) {
      const query = filterText.toLowerCase();
      const haystack = [
        result.check,
        ...(result.files || []),
        ...(result.reads || []),
        ...(result.regions || []),
        ...(result.assessment || []).map((assessment) => assessment.code + " " + assessment.description),
      ]
        .join(" ")
        .toLowerCase();
      if (!haystack.includes(query)) {
        return false;
      }
    }
    return true;
  }

  function render() {
    const generatedAt = window.SEQCHECK_GENERATED_AT || "";
    const sorted = results
      .map((result, index) => ({ result, index }))
      .sort((left, right) => {
        const severityDiff =
          (SEV[worstSev(left.result.assessment || [])] || 3) -
          (SEV[worstSev(right.result.assessment || [])] || 3);
        if (severityDiff !== 0) {
          return severityDiff;
        }
        const checkDiff = (left.result.check || "").localeCompare(right.result.check || "");
        if (checkDiff !== 0) {
          return checkDiff;
        }
        return left.index - right.index;
      });

    const counts = { error: 0, warning: 0, interpretation: 0, pass: 0 };
    sorted.forEach((entry) => {
      counts[worstSev(entry.result.assessment || [])] += 1;
    });

    let html = "";
    html += `<div class="hdr">
      <div class="hdr-title">seqcheck report</div>
      <div class="hdr-row">
        <span class="l">spec</span> ${esc(meta.spec)}
        <span class="sep">|</span>
        <span class="l">modality</span> ${esc(meta.modality)}
        <span class="sep">|</span>
        <span class="l">reads sampled</span> ${(meta.requested_reads || 0).toLocaleString()}
        <span class="sep">|</span>
        <span class="l">schema</span> ${esc(report.report_schema_version)}
        <span class="sep">|</span>
        <span class="l">generated</span> ${esc(generatedAt)}
      </div>
    </div>`;

    const invocation = buildInvocation();
    const reportInvocation = buildReportInvocation();

    html += `<div class="meta-section">
      <div class="meta-section-head">Run</div>
      <div class="meta-section-body">
        <div class="run-row"><span class="run-key">tool</span><span class="run-value"><a class="inline-link" href="${esc(repositoryUrl)}" target="_blank" rel="noreferrer">seqcheck</a>${seqcheckVersion ? " " + esc(seqcheckVersion) : ""}</span></div>
        <div class="run-row"><span class="run-key">report</span><span class="run-value">Generated from <span style="font-family:var(--mono)">${esc(meta.command || "check")}</span> JSON output.</span></div>
        <div class="run-row"><span class="run-key">source</span><span class="run-value"><a class="inline-link" href="${esc(repositoryUrl)}" target="_blank" rel="noreferrer">${esc(repositoryUrl)}</a></span></div>
        <div class="run-row"><span class="run-key">check</span><span class="run-value">Command used to create the JSON report.</span></div>
        <div class="run-command">${esc(invocation)}</div>
        <div class="run-row"><span class="run-key">report</span><span class="run-value">Command used to generate this HTML file.</span></div>
        <div class="run-command">${esc(reportInvocation)}</div>
      </div>
    </div>`;

    if (libData && libData.regions && libData.regions.length) {
      html += `<div class="mol-section">
        <div class="mol-section-head">Library Structure \u2014 ${esc(libData.assay_name)} (${esc(libData.modality)})</div>
        <div class="mol-body">${buildMolSvg()}</div>
        <div class="mol-legend">
          <span><span class="leg-swatch" style="background:var(--reg-fixed)"></span>fixed</span>
          <span><span class="leg-swatch" style="background:var(--reg-onlist)"></span>onlist</span>
          <span><span class="leg-swatch" style="background:var(--reg-random)"></span>random</span>
        </div>
      </div>`;
    }

    html += '<div class="summary">';
    [
      ["error", "err"],
      ["warning", "wrn"],
      ["interpretation", "int"],
      ["pass", "pas"],
    ].forEach(([severity, className]) => {
      const count = counts[severity] || 0;
      const active = filterSev === severity ? " active" : "";
      const zero = count === 0 ? " zero" : "";
      html += `<div class="s-pill ${className}${active}${zero}" data-sev="${severity}"><span class="n">${count}</span><span class="t">${sevLabel(severity)}</span></div>`;
    });
    html += "</div>";

    html += `<div class="search-bar"><input type="text" id="q" placeholder="Search checks, regions, descriptions\u2026" value="${esc(filterText)}"></div>`;

    const visible = sorted.filter((entry) => matches(entry.result));
    html += '<div class="result-list">';
    if (!visible.length) {
      html += '<div class="empty-msg">No results match the current filters.</div>';
    }
    visible.forEach((entry, index) => {
      const result = entry.result;
      const severity = worstSev(result.assessment || []);
      const uid = "ri" + index;
      html += `<div class="ri"><div class="ri-head" data-detail="${uid}">`;
      html += `<span class="ri-badge ${severity}">${sevLabel(severity)}</span>`;
      html += `<div class="ri-body"><div class="ri-desc">${esc(naturalDesc(result))}</div><div class="ri-context">${esc(contextLine(result))}</div></div>`;
      html += `<span class="ri-arrow" id="arr-${uid}">\u25b6</span></div>`;
      html += `<div class="ri-detail" id="${uid}">`;

      if (result.assessment && result.assessment.length) {
        html += '<div class="detail-section"><div class="detail-label">Assessments</div>';
        result.assessment
          .slice()
          .sort((left, right) => (SEV[left.type] || 3) - (SEV[right.type] || 3))
          .forEach((assessment) => {
            html += `<div style="margin-bottom:3px"><span class="ri-badge ${assessment.type}" style="font-size:9px;padding:1px 6px;width:auto;display:inline-block">${sevLabel(assessment.type)}</span> <span style="font-family:var(--mono);font-size:11px;color:var(--text-3)">${esc(assessment.code)}</span> <span style="font-size:12px;color:var(--text-2)">${esc(assessment.description)}</span></div>`;
          });
        html += "</div>";
      }

      if (result.expected && result.expected.length) {
        html += '<div class="detail-section"><div class="detail-label">Expected</div>';
        result.expected.forEach((item) => {
          html += renderItem(item);
        });
        html += "</div>";
      }

      if (result.observed && result.observed.length) {
        html += '<div class="detail-section"><div class="detail-label">Observed</div>';
        result.observed.forEach((item) => {
          html += renderItem(item);
        });
        html += "</div>";
      }

      html += "</div></div>";
    });
    html += "</div>";

    app.innerHTML = html;
    bind();
  }

  function bind() {
    document.querySelectorAll(".s-pill").forEach((element) => {
      element.addEventListener("click", () => {
        const severity = element.getAttribute("data-sev");
        filterSev = filterSev === severity ? null : severity;
        render();
      });
    });

    const query = document.getElementById("q");
    if (query) {
      query.addEventListener("input", () => {
        filterText = query.value;
        render();
        query.focus();
      });
      query.setSelectionRange(query.value.length, query.value.length);
    }

    document.querySelectorAll(".ri-head[data-detail]").forEach((element) => {
      element.addEventListener("click", () => {
        const uid = element.getAttribute("data-detail");
        const detail = document.getElementById(uid);
        const arrow = document.getElementById("arr-" + uid);
        if (detail) {
          detail.classList.toggle("open");
        }
        if (arrow) {
          arrow.classList.toggle("open");
        }
      });
    });

    document.querySelectorAll(".tbl-toggle[data-tbl]").forEach((element) => {
      element.addEventListener("click", (event) => {
        event.stopPropagation();
        const uid = element.getAttribute("data-tbl");
        const inner = document.getElementById(uid);
        if (inner) {
          inner.classList.toggle("open");
        }
      });
    });

    const tip = document.getElementById("mol-tip");
    if (!tip || !libData || !libData.regions) {
      return;
    }

    document.querySelectorAll(".region-rect").forEach((element) => {
      const regionId = element.getAttribute("data-region");
      const region = libData.regions.find((candidate) => candidate.region_id === regionId);
      if (!region) {
        return;
      }
      element.addEventListener("mouseenter", () => {
        let html = `<div class="tip-name">${esc(region.name)}</div>`;
        html += `<div>${esc(region.region_type)} \u00b7 ${esc(region.sequence_type)} \u00b7 ${esc(region.len)} bp</div>`;
        if (region.sequence && region.sequence.length <= 30) {
          html += `<div class="tip-dim">${esc(region.sequence)}</div>`;
        } else if (region.sequence) {
          html += `<div class="tip-dim">${esc(region.sequence.slice(0, 25))}\u2026</div>`;
        }
        tip.innerHTML = html;
        tip.classList.add("show");
      });
      element.addEventListener("mousemove", (event) => {
        tip.style.left = event.clientX + 12 + "px";
        tip.style.top = event.clientY - 10 + "px";
      });
      element.addEventListener("mouseleave", () => {
        tip.classList.remove("show");
      });
    });
  }

  render();
})();
