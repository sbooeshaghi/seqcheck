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
  const selected = {};

  let filterText = "";
  let filterSev = null;

  function esc(value) {
    if (value == null || value === "") {
      return '<span class="val-null">\u2014</span>';
    }
    const div = document.createElement("div");
    div.textContent = String(value);
    return div.innerHTML;
  }

  function escAttr(value) {
    const div = document.createElement("div");
    div.textContent = value == null ? "" : String(value);
    return div.innerHTML;
  }

  function fmtVal(value) {
    if (value == null || value === "") {
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
    const assessments = (result.assessment || []).filter(
      (assessment) => assessment.type !== "pass",
    );
    if (assessments.length) {
      assessments.sort(
        (left, right) => (SEV[left.type] || 3) - (SEV[right.type] || 3),
      );
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
    const supplied = findMetric(
      inputCheck ? inputCheck.observed : [],
      "supplied_input_paths",
    );
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

  function renderMetricItem(item) {
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
      html += `<div class="tbl-inner" id="${uid}"><div class="tbl-wrap"><table class="dtbl"><thead><tr>${keys
        .map((key) => "<th>" + esc(key) + "</th>")
        .join("")}</tr></thead><tbody>`;
      data.value.forEach((record) => {
        html +=
          "<tr>" +
          keys.map((key) => "<td>" + fmtVal(record[key]) + "</td>").join("") +
          "</tr>";
      });
      html += "</tbody></table></div></div>";
      return html;
    }
    return "";
  }

  function lengthLabel(minLen, maxLen) {
    if (minLen === maxLen) {
      return `${fmtVal(minLen)} bp`;
    }
    return `${fmtVal(minLen)}-${fmtVal(maxLen)} bp`;
  }

  function bpRangeLabel(start, end) {
    return `${fmtVal(start)}-${fmtVal(end)} bp`;
  }

  function kvList(rows) {
    return `<div class="kv-list">${rows
      .map(
        ([key, value, mono]) =>
          `<div class="kv-row"><div class="kv-key">${esc(
            key,
          )}</div><div class="kv-value${mono ? " mono" : ""}">${value}</div></div>`,
      )
      .join("")}</div>`;
  }

  function detailSection(title, body) {
    if (!body) {
      return "";
    }
    return `<section class="detail-section"><div class="detail-title">${esc(
      title,
    )}</div><div class="detail-content">${body}</div></section>`;
  }

  function detailTable(title, rows, columns) {
    if (!rows || !rows.length) {
      return "";
    }
    const keys =
      columns ||
      Array.from(
        rows.reduce((seen, row) => {
          Object.keys(row).forEach((key) => seen.add(key));
          return seen;
        }, new Set()),
      );
    return detailSection(
      title,
      `<div class="table-wrap"><table><thead><tr>${keys
        .map((key) => `<th>${esc(key)}</th>`)
        .join("")}</tr></thead><tbody>${rows
        .map(
          (row) =>
            `<tr>${keys.map((key) => `<td>${fmtVal(row[key])}</td>`).join("")}</tr>`,
        )
        .join("")}</tbody></table></div>`,
    );
  }

  function detailShell(header, sections) {
    return `<div class="detail-shell"><div class="detail-shell-head">${header}</div><div class="detail-shell-body">${sections.join(
      "",
    )}</div></div>`;
  }

  function pathLabel(pathNames) {
    return (pathNames || []).join(" / ");
  }

  function currentSelection() {
    if (!libData) {
      return null;
    }
    return selected[libData.modality] || null;
  }

  function visibleRegionsForRead(read) {
    return (libData.regions || []).filter(
      (region) => region.bp_start < read.end && region.bp_end > read.start,
    );
  }

  function visibleReadsForRegion(region) {
    return (libData.reads || []).filter(
      (read) => region.bp_start < read.end && region.bp_end > read.start,
    );
  }

  function overlappingReads(read) {
    return (libData.reads || []).filter(
      (other) =>
        other.read_id !== read.read_id &&
        Math.max(read.start, other.start) < Math.min(read.end, other.end),
    );
  }

  function overlapRows() {
    const rows = [];
    const reads = libData.reads || [];
    for (let i = 0; i < reads.length; i += 1) {
      for (let j = i + 1; j < reads.length; j += 1) {
        const left = reads[i];
        const right = reads[j];
        const start = Math.max(left.start, right.start);
        const end = Math.min(left.end, right.end);
        if (start < end) {
          rows.push({
            left_read: left.label || left.read_id,
            right_read: right.label || right.read_id,
            bp_range: bpRangeLabel(start, end),
            overlap_bp: end - start,
          });
        }
      }
    }
    return rows;
  }

  function regionTooltip(region) {
    let html = `<div class="tip-name">${esc(region.name)}</div>`;
    html += `<div>${esc(region.region_type)} \u00b7 ${esc(
      region.sequence_type,
    )} \u00b7 ${esc(lengthLabel(region.min_len, region.max_len))}</div>`;
    html += `<div class="tip-dim">${esc(bpRangeLabel(region.bp_start, region.bp_end))}</div>`;
    if (!region.is_leaf) {
      html += `<div class="tip-dim">${esc(region.child_region_ids.length)} child regions</div>`;
    }
    if (region.sequence) {
      const preview =
        region.sequence.length <= 40
          ? region.sequence
          : `${region.sequence.slice(0, 37)}\u2026`;
      html += `<div class="tip-dim">${esc(preview)}</div>`;
    }
    return html;
  }

  function readTooltip(read) {
    let html = `<div class="tip-name">${esc(read.label || read.read_id)}</div>`;
    html += `<div>${esc(read.read_id)} \u00b7 ${esc(read.strand)} \u00b7 ${esc(
      lengthLabel(read.min_len, read.max_len),
    )}</div>`;
    html += `<div class="tip-dim">${esc(bpRangeLabel(read.start, read.end))} anchored at ${esc(
      read.primer_id,
    )}</div>`;
    return html;
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

  function regionIsSelected(region, selection) {
    if (!selection || selection.kind !== "region") {
      return false;
    }
    return (region.path_region_ids || []).includes(selection.id);
  }

  function resultMatchesSelection(result, selection) {
    if (!selection) {
      return false;
    }
    if (selection.kind === "region") {
      return (result.regions || []).includes(selection.id);
    }
    if (selection.kind === "read") {
      return (result.reads || []).includes(selection.id);
    }
    return false;
  }

  function buildMolSvg() {
    if (!libData || !libData.regions || !libData.regions.length) {
      return "";
    }

    const leafRegions = libData.regions || [];
    const regionNodes = libData.region_nodes || [];
    const groupRegions = regionNodes.filter((region) => !region.is_leaf);
    const reads = libData.reads || [];
    const totalBp =
      libData.total_bp || leafRegions.reduce((sum, region) => sum + region.len, 0);
    const selection = currentSelection();

    const pad = { left: 10, right: 10 };
    const molW = 860;
    const barY = 70;
    const barH = 22;
    const svgW = molW + pad.left + pad.right;
    const bpScale = molW / Math.max(totalBp, 1);
    const minPx = 18;
    const groupGap = 12;
    const groupHeight = 8;
    const groupLevels = groupRegions.length
      ? Math.max(...groupRegions.map((region) => region.depth || 0)) + 1
      : 0;
    const groupTrackTop = barY - groupLevels * groupGap - 4;
    const rects = [];

    let rawWidths = leafRegions.map((region) => region.len * bpScale);
    let deficit = 0;
    rawWidths = rawWidths.map((px) => {
      if (px < minPx) {
        deficit += minPx - px;
        return minPx;
      }
      return px;
    });
    if (deficit > 0) {
      const shrinkTotal = rawWidths
        .filter((px) => px > minPx)
        .reduce((sum, px) => sum + px, 0);
      if (shrinkTotal > 0) {
        rawWidths = rawWidths.map((px) =>
          px > minPx ? px - deficit * (px / shrinkTotal) : px,
        );
      }
    }

    let currentX = pad.left;
    leafRegions.forEach((region, index) => {
      const width = rawWidths[index];
      rects.push({ ...region, x: currentX, w: width });
      currentX += width;
    });

    function bpToX(bp) {
      if (bp <= 0) {
        return pad.left;
      }
      for (const region of rects) {
        if (bp <= region.bp_end) {
          const fraction = region.len === 0 ? 0 : (bp - region.bp_start) / region.len;
          return region.x + fraction * region.w;
        }
      }
      const last = rects[rects.length - 1];
      return last ? last.x + last.w : pad.left;
    }

    let svg = "";
    svg += `<line x1="${pad.left}" y1="${barY - 1}" x2="${currentX}" y2="${barY - 1}" stroke="var(--border)" stroke-width="0.5"/>`;
    svg += `<line x1="${pad.left}" y1="${barY + barH + 1}" x2="${currentX}" y2="${barY + barH + 1}" stroke="var(--border)" stroke-width="0.5"/>`;

    groupRegions.forEach((region) => {
      const x = bpToX(region.bp_start);
      const width = Math.max(bpToX(region.bp_end) - x, 1);
      const y = groupTrackTop + (region.depth || 0) * groupGap;
      const selectedClass = regionIsSelected(region, selection) ? " selected" : "";
      svg += `<rect class="group-rect${selectedClass}" data-kind="region" data-modality="${esc(
        libData.modality,
      )}" data-id="${esc(region.region_id)}" x="${x}" y="${y}" width="${width}" height="${groupHeight}" rx="2" />`;
      if (width > 42) {
        svg += `<text class="group-label" x="${x + 3}" y="${y - 2}">${esc(region.name)}</text>`;
      }
    });

    rects.forEach((region) => {
      const selectedClass = regionIsSelected(region, selection) ? " selected" : "";
      svg += `<rect class="region-rect${selectedClass}" data-kind="region" data-modality="${esc(
        libData.modality,
      )}" data-id="${esc(region.region_id)}" x="${region.x}" y="${barY}" width="${region.w}" height="${barH}" rx="2" fill="${seqTypeColor(
        region.sequence_type,
      )}" stroke="${seqTypeStroke(region.sequence_type)}" stroke-width="0.5"/>`;
    });

    const labelY = barY + barH + 10;
    rects.forEach((region) => {
      const centerX = region.x + region.w / 2;
      if (region.w > 14) {
        svg += `<text class="region-label" x="${centerX}" y="${labelY}" text-anchor="end" transform="rotate(-40 ${centerX} ${labelY})">${esc(
          region.name,
        )}</text>`;
      }
    });

    svg += `<text class="bp-label" x="${pad.left}" y="${barY - 3}" text-anchor="middle">0</text>`;
    rects.forEach((region) => {
      svg += `<text class="bp-label" x="${region.x + region.w}" y="${barY - 3}" text-anchor="middle">${region.bp_end}</text>`;
    });

    const readColors = ["#1e40af", "#059669", "#d97706", "#dc2626", "#7c3aed"];
    const posReads = reads.filter((read) => read.strand === "pos");
    const negReads = reads.filter((read) => read.strand === "neg");

    function drawRead(read, yBase, above, colorIndex) {
      const color = readColors[colorIndex % readColors.length];
      const x1 = bpToX(read.start);
      const x2 = bpToX(read.end);
      const arrowSize = 5;
      const selectedClass =
        selection && selection.kind === "read" && selection.id === read.read_id
          ? " selected"
          : "";
      const groupAttrs = `class="read-group${selectedClass}" data-kind="read" data-modality="${esc(
        libData.modality,
      )}" data-id="${esc(read.read_id)}"`;

      if (above) {
        const y = yBase;
        svg += `<g ${groupAttrs}><line x1="${x1}" y1="${y}" x2="${Math.max(
          x1,
          x2 - arrowSize,
        )}" y2="${y}" stroke="${color}" class="read-line"/><polygon points="${x2},${y} ${
          x2 - arrowSize
        },${y - arrowSize} ${x2 - arrowSize},${y + arrowSize}" fill="${color}" stroke="${color}"/><line x1="${x1}" y1="${y}" x2="${x1}" y2="${barY}" stroke="${color}" stroke-width="1" stroke-dasharray="2,2" opacity="0.4"/><text class="read-label" x="${
          x1 + 3
        }" y="${y - 5}" fill="${color}">${esc(read.label || read.read_id)}</text></g>`;
      } else {
        const y = yBase;
        svg += `<g ${groupAttrs}><line x1="${x2}" y1="${y}" x2="${Math.min(
          x2,
          x1 + arrowSize,
        )}" y2="${y}" stroke="${color}" class="read-line"/><polygon points="${x1},${y} ${
          x1 + arrowSize
        },${y - arrowSize} ${x1 + arrowSize},${y + arrowSize}" fill="${color}" stroke="${color}"/><line x1="${x2}" y1="${y}" x2="${x2}" y2="${
          barY + barH
        }" stroke="${color}" stroke-width="1" stroke-dasharray="2,2" opacity="0.4"/><text class="read-label" x="${
          x2 - 3
        }" y="${y + 13}" text-anchor="end" fill="${color}">${esc(
          read.label || read.read_id,
        )}</text></g>`;
      }
    }

    let posY = groupTrackTop - 14;
    posReads.forEach((read, index) => {
      drawRead(read, posY, true, index);
      posY -= 22;
    });

    let negY = barY + barH + 48;
    negReads.forEach((read, index) => {
      drawRead(read, negY, false, posReads.length + index);
      negY += 22;
    });

    const svgH = Math.max(negY + 10, barY + barH + 50);
    return `<svg class="mol-svg" viewBox="0 0 ${svgW} ${svgH}" width="100%" xmlns="http://www.w3.org/2000/svg">${svg}</svg>`;
  }

  function selectorRow(config) {
    const {
      kind,
      id,
      label,
      sub,
      meta,
      active,
      depth = 0,
      nodeType = "",
      hasChildren = false,
    } = config;
    const padding = kind === "region" ? 10 + depth * 16 : 10;
    const marker = kind === "region" ? (hasChildren ? "\u25a1" : "\u2022") : "\u2192";
    const typeClass = nodeType ? ` ${nodeType}` : "";

    return `<button class="selector-row${active ? " active" : ""}${typeClass}" data-kind="${esc(
      kind,
    )}" data-modality="${esc(libData.modality)}" data-id="${esc(id)}"><span class="selector-main"><span class="selector-label-line" style="padding-left:${padding}px"><span class="selector-marker">${marker}</span><span class="selector-label">${esc(
      label,
    )}</span></span>${
      sub
        ? `<span class="selector-sub" style="padding-left:${padding + 16}px">${esc(sub)}</span>`
        : ""
    }</span><span class="selector-meta">${esc(meta)}</span></button>`;
  }

  function selectorHtml() {
    if (!libData) {
      return "";
    }
    const selection = currentSelection();
    const regionRows = (libData.region_nodes || []).map((region) =>
      selectorRow({
        kind: "region",
        id: region.region_id,
        label: region.name,
        sub: `${region.region_type} \u00b7 ${region.sequence_type}`,
        meta: region.is_leaf
          ? bpRangeLabel(region.bp_start, region.bp_end)
          : `${bpRangeLabel(region.bp_start, region.bp_end)} \u00b7 ${
              region.child_region_ids.length
            } children`,
        active:
          selection &&
          selection.kind === "region" &&
          selection.id === region.region_id,
        depth: region.depth || 0,
        nodeType: region.is_leaf ? "leaf" : "branch",
        hasChildren: !region.is_leaf,
      }),
    );

    const readRows = (libData.reads || []).map((read) =>
      selectorRow({
        kind: "read",
        id: read.read_id,
        label: read.label || read.read_id,
        sub: `${read.strand} \u00b7 ${read.primer_id}`,
        meta: bpRangeLabel(read.start, read.end),
        active: selection && selection.kind === "read" && selection.id === read.read_id,
      }),
    );

    return `<div class="selector-pane"><div class="selector-group"><div class="selector-head">regions</div><div class="selector-body">${
      regionRows.join("") || '<div class="empty-state">No regions.</div>'
    }</div></div><div class="selector-group"><div class="selector-head">reads</div><div class="selector-body">${
      readRows.join("") || '<div class="empty-state">No reads.</div>'
    }</div></div></div>`;
  }

  function modalitySummary() {
    const rows = [
      ["assay id", esc(libData.assay_id), true],
      ["modality", esc(libData.modality), true],
      ["library region", esc(libData.library_region_id), true],
      ["seqspec version", esc(libData.seqspec_version || ""), true],
      ["total length", esc(`${libData.total_bp} bp`), true],
      ["region count", esc((libData.region_nodes || []).length), true],
      ["read count", esc((libData.reads || []).length), true],
    ];
    const sections = [
      detailSection(
        "summary",
        `<div class="selection-note">Select a region or read to inspect its metadata. Matching seqcheck results remain below.</div>${kvList(
          rows,
        )}`,
      ),
      detailTable("sequence protocols", libData.sequence_protocols, ["protocol_id", "name"]),
      detailTable("sequence kits", libData.sequence_kits, ["kit_id", "name"]),
      detailTable("library protocols", libData.library_protocols, ["protocol_id", "name"]),
      detailTable("library kits", libData.library_kits, ["kit_id", "name"]),
      detailTable("overlapping reads", overlapRows(), [
        "left_read",
        "right_read",
        "bp_range",
        "overlap_bp",
      ]),
    ].filter(Boolean);
    return detailShell("modality", sections);
  }

  function regionDetails(region) {
    const children = (libData.region_nodes || []).filter(
      (node) => node.parent_region_id === region.region_id,
    );
    const visibleReads = visibleReadsForRegion(region);
    const rows = [
      ["region id", esc(region.region_id), true],
      ["name", esc(region.name), false],
      ["path", esc(pathLabel(region.path_names)), true],
      ["region type", esc(region.region_type), true],
      ["sequence type", esc(region.sequence_type), true],
      ["length", esc(lengthLabel(region.min_len, region.max_len)), true],
      ["bp range", esc(bpRangeLabel(region.bp_start, region.bp_end)), true],
      ["leaf", esc(region.is_leaf ? "true" : "false"), true],
      ["children", esc(region.child_region_ids.length), true],
      ["visible in reads", esc(visibleReads.length), true],
    ];
    const sections = [
      detailSection("metadata", kvList(rows)),
      region.sequence
        ? detailSection("sequence", `<pre class="region-seq">${esc(region.sequence)}</pre>`)
        : "",
      children.length
        ? detailTable(
            "child regions",
            children.map((child) => ({
              region_id: child.region_id,
              name: child.name,
              region_type: child.region_type,
              sequence_type: child.sequence_type,
              bp_range: bpRangeLabel(child.bp_start, child.bp_end),
              length: child.len,
            })),
            ["region_id", "name", "region_type", "sequence_type", "bp_range", "length"],
          )
        : "",
      visibleReads.length
        ? detailTable(
            "visible reads",
            visibleReads.map((read) => ({
              read_id: read.read_id,
              name: read.label || read.read_id,
              strand: read.strand,
              bp_range: bpRangeLabel(read.start, read.end),
            })),
            ["read_id", "name", "strand", "bp_range"],
          )
        : "",
      region.onlist
        ? detailTable("onlist", [region.onlist], [
            "file_id",
            "filename",
            "filetype",
            "urltype",
            "url",
            "md5",
          ])
        : "",
    ].filter(Boolean);
    return detailShell(`region \u00b7 ${esc(region.name)}`, sections);
  }

  function readDetails(read) {
    const regions = visibleRegionsForRead(read);
    const overlaps = overlappingReads(read);
    const rows = [
      ["read id", esc(read.read_id), true],
      ["name", esc(read.name), false],
      ["strand", esc(read.strand), true],
      ["primer id", esc(read.primer_id), true],
      ["length", esc(lengthLabel(read.min_len, read.max_len)), true],
      ["bp range", esc(bpRangeLabel(read.start, read.end)), true],
      ["visible regions", esc(regions.length), true],
      ["overlapping reads", esc(overlaps.length), true],
    ];
    const sections = [
      detailSection("metadata", kvList(rows)),
      detailTable(
        "visible regions",
        regions.map((region) => ({
          region_id: region.region_id,
          name: region.name,
          region_type: region.region_type,
          sequence_type: region.sequence_type,
          bp_range: bpRangeLabel(region.bp_start, region.bp_end),
        })),
        ["region_id", "name", "region_type", "sequence_type", "bp_range"],
      ),
      detailTable(
        "overlapping reads",
        overlaps.map((other) => ({
          read_id: other.read_id,
          name: other.label || other.read_id,
          strand: other.strand,
          bp_range: bpRangeLabel(other.start, other.end),
        })),
        ["read_id", "name", "strand", "bp_range"],
      ),
      detailTable("files", read.files, [
        "file_id",
        "filename",
        "filetype",
        "urltype",
        "url",
        "md5",
      ]),
    ].filter(Boolean);
    return detailShell(`read \u00b7 ${esc(read.label || read.read_id)}`, sections);
  }

  function selectionHtml() {
    if (!libData) {
      return "";
    }
    const selection = currentSelection();
    if (!selection) {
      return modalitySummary();
    }
    if (selection.kind === "region") {
      const region = (libData.region_nodes || []).find(
        (node) => node.region_id === selection.id,
      );
      if (region) {
        return regionDetails(region);
      }
    }
    if (selection.kind === "read") {
      const read = (libData.reads || []).find(
        (item) => item.read_id === selection.id,
      );
      if (read) {
        return readDetails(read);
      }
    }
    return modalitySummary();
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
        ...(result.assessment || []).map(
          (assessment) => assessment.code + " " + assessment.description,
        ),
      ]
        .join(" ")
        .toLowerCase();
      if (!haystack.includes(query)) {
        return false;
      }
    }
    return true;
  }

  function render(options = {}) {
    const restoreSearchFocus = Boolean(options.restoreSearchFocus);
    const restoreSearchPosition =
      typeof options.restoreSearchPosition === "number"
        ? options.restoreSearchPosition
        : null;
    const generatedAt = window.SEQCHECK_GENERATED_AT || "";
    const selection = currentSelection();
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
          <span><span class="leg-swatch outline"></span>nested region span</span>
        </div>
        <div class="detail-layout">${selectorHtml()}<div class="detail-pane">${selectionHtml()}</div></div>
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

    html += `<div class="search-bar"><input type="text" id="q" placeholder="Search checks, regions, descriptions\u2026" value="${escAttr(filterText)}"></div>`;

    const visible = sorted.filter((entry) => matches(entry.result));
    html += '<div class="result-list">';
    if (!visible.length) {
      html += '<div class="empty-msg">No results match the current filters.</div>';
    }
    visible.forEach((entry, index) => {
      const result = entry.result;
      const severity = worstSev(result.assessment || []);
      const uid = "ri" + index;
      const relatedClass = resultMatchesSelection(result, selection) ? " related" : "";
      html += `<div class="ri${relatedClass}"><div class="ri-head" data-detail="${uid}">`;
      html += `<span class="ri-badge ${severity}">${sevLabel(severity)}</span>`;
      html += `<div class="ri-body"><div class="ri-desc">${esc(
        naturalDesc(result),
      )}</div><div class="ri-context">${esc(contextLine(result))}</div></div>`;
      html += `<span class="ri-arrow" id="arr-${uid}">\u25b6</span></div>`;
      html += `<div class="ri-detail" id="${uid}">`;

      if (result.assessment && result.assessment.length) {
        html += '<div class="detail-section"><div class="detail-label">Assessments</div>';
        result.assessment
          .slice()
          .sort((left, right) => (SEV[left.type] || 3) - (SEV[right.type] || 3))
          .forEach((assessment) => {
            html += `<div style="margin-bottom:3px"><span class="ri-badge ${assessment.type}" style="font-size:9px;padding:1px 6px;width:auto;display:inline-block">${sevLabel(assessment.type)}</span> <span style="font-family:var(--mono);font-size:11px;color:var(--text-3)">${esc(
              assessment.code,
            )}</span> <span style="font-size:12px;color:var(--text-2)">${esc(
              assessment.description,
            )}</span></div>`;
          });
        html += "</div>";
      }

      if (result.expected && result.expected.length) {
        html += '<div class="detail-section"><div class="detail-label">Expected</div>';
        result.expected.forEach((item) => {
          html += renderMetricItem(item);
        });
        html += "</div>";
      }

      if (result.observed && result.observed.length) {
        html += '<div class="detail-section"><div class="detail-label">Observed</div>';
        result.observed.forEach((item) => {
          html += renderMetricItem(item);
        });
        html += "</div>";
      }

      html += "</div></div>";
    });
    html += "</div>";

    app.innerHTML = html;
    bind();

    if (restoreSearchFocus) {
      const query = document.getElementById("q");
      if (query) {
        query.focus();
        const position =
          restoreSearchPosition == null ? query.value.length : restoreSearchPosition;
        query.setSelectionRange(position, position);
      }
    }
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
        const cursorPosition = query.selectionStart;
        render({
          restoreSearchFocus: true,
          restoreSearchPosition: cursorPosition,
        });
      });
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

    document.querySelectorAll("[data-kind][data-modality][data-id]").forEach((element) => {
      element.addEventListener("click", (event) => {
        event.preventDefault();
        event.stopPropagation();
        const modality = element.getAttribute("data-modality");
        const kind = element.getAttribute("data-kind");
        const id = element.getAttribute("data-id");
        const existing = selected[modality];
        if (existing && existing.kind === kind && existing.id === id) {
          delete selected[modality];
        } else {
          selected[modality] = { kind, id };
        }
        render();
      });
    });

    const tip = document.getElementById("mol-tip");
    document.querySelectorAll(".region-rect[data-id], .group-rect[data-id]").forEach((element) => {
      const regionId = element.getAttribute("data-id");
      const region = (libData.region_nodes || []).find(
        (item) => item.region_id === regionId,
      );
      if (!region) {
        return;
      }
      element.addEventListener("mouseenter", () => {
        tip.innerHTML = regionTooltip(region);
        tip.classList.add("show");
      });
      element.addEventListener("mousemove", (event) => {
        tip.style.left = `${event.clientX + 12}px`;
        tip.style.top = `${event.clientY - 10}px`;
      });
      element.addEventListener("mouseleave", () => {
        tip.classList.remove("show");
      });
    });

    document.querySelectorAll(".read-group[data-id]").forEach((element) => {
      const readId = element.getAttribute("data-id");
      const read = (libData.reads || []).find((item) => item.read_id === readId);
      if (!read) {
        return;
      }
      element.addEventListener("mouseenter", () => {
        tip.innerHTML = readTooltip(read);
        tip.classList.add("show");
      });
      element.addEventListener("mousemove", (event) => {
        tip.style.left = `${event.clientX + 12}px`;
        tip.style.top = `${event.clientY - 10}px`;
      });
      element.addEventListener("mouseleave", () => {
        tip.classList.remove("show");
      });
    });
  }

  render();
})();
