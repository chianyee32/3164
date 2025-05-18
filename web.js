// File: web.js

// ─── 1. Configuration ──────────────────────────────────────────────────────────

// Always post back to the same origin
const UPLOAD_ENDPOINT = `${window.location.origin}/upload_prediction`;

// ─── 2. State ──────────────────────────────────────────────────────────────────

let allRows = [];
let currentSortColumn = null;
let sortAscending = true;

// ─── 3. DOM Ready ──────────────────────────────────────────────────────────────

document.addEventListener("DOMContentLoaded", () => {
  // File input
  const fileInput = document.getElementById("file_upload_prediction");
  fileInput.addEventListener("change", event =>
    uploadFile(event)
  );

  // Download filtered
  document
    .getElementById("download_filtered")
    .addEventListener("click", downloadFilteredData);
});

// ─── 4. Upload & Preview ───────────────────────────────────────────────────────

function uploadFile(event) {
  const file = event.target.files[0];
  if (!file || (!file.name.endsWith(".csv") && file.type !== "text/csv")) {
    alert("⚠️ Please upload a valid CSV file.");
    return;
  }

  // show spinner + hide sections
  const infoDiv       = document.getElementById("file_info_prediction");
  const loadingDiv    = document.getElementById("loading_indicator");
  const uploadSection = document.getElementById("upload-section");
  const tooltipBox    = document.getElementById("tooltip-box");
  const formatBtn     = document.getElementById("format_btn");
  const spinner       = document.getElementById("spinner");
  const preview       = document.getElementById("csv_preview");

  uploadSection.style.display = "none";
  tooltipBox.style.display    = "none";
  formatBtn.style.display     = "none";
  spinner.style.display       = "block";
  preview.innerHTML           = "";

  infoDiv.innerHTML = `📄 File Selected: <span style="color: teal;">${file.name}</span>`;
  loadingDiv.innerHTML = `⏳ <span style="color: #9b59b6;">Processing your file, please wait...</span>`;
  loadingDiv.style.display = "block";

  const formData = new FormData();
  formData.append("file", file);

  fetch(UPLOAD_ENDPOINT, { method: "POST", body: formData })
    .then(res => res.json())
    .then(data => {
      // hide spinner
      spinner.style.display    = "none";
      loadingDiv.style.display = "none";

      if (data.error) {
        // show error
        let msg = `❌ Error: ${data.error}`;
        if (data.details) msg += `\n🔍 Details:\n${data.details}`;
        infoDiv.innerText += `\n${msg}`;
        uploadSection.style.display = "block";
        return;
      }

      // success: add download link
      const link = document.createElement("a");
      link.href        = data.download_url;
      link.textContent = "⬇️ Download Prediction File";
      link.target      = "_blank";
      link.style.cssText =
        "display:inline-block;margin-top:15px;" +
        "background-color:#2ecc71;color:white;" +
        "padding:10px 20px;border-radius:8px;text-decoration:none;";
      infoDiv.appendChild(document.createElement("br"));
      infoDiv.appendChild(link);

      // fetch the CSV for preview
      return fetch(data.download_url).then(r => r.text());
    })
    .then(csvText => {
      if (!csvText) return;
      displayCSV(csvText);
      createCancerTypeCheckboxes(allRows);
      document.getElementById("csv_preview")
              .scrollIntoView({ behavior: "smooth" });
    })
    .catch(() => {
      // network / unexpected error
      spinner.style.display    = "none";
      loadingDiv.style.display = "none";
      infoDiv.innerText        = "❌ Upload failed.";
      uploadSection.style.display = "block";
    });
}

// ─── 5. CSV → Table ────────────────────────────────────────────────────────────

function displayCSV(csvText) {
  allRows = csvText.trim().split("\n").map(r => r.split(","));
  displayTable(allRows);
}

function displayTable(rows) {
  // hide tooltip if visible
  const tooltip = document.getElementById("tooltip-box");
  if (tooltip) tooltip.style.display = "none";

  const preview = document.getElementById("csv_preview");
  preview.innerHTML = "";

  const table = document.createElement("table");
  table.className = "styled-table";

  const thead = document.createElement("thead");
  const tbody = document.createElement("tbody");

  // human-readable headers
  const cmap = {
    DRUG_ID: "Drug ID",
    DRUG_NAME: "Drug Name",
    COSMIC_ID: "Cosmic ID",
    CCLE_Name: "Cell Line Name",
    CANCER_TYPE: "Cancer Type",
    Predicted_LN_IC50: "Predicted LN IC50"
  };

  rows.forEach((row, i) => {
    const tr = document.createElement("tr");
    row.forEach((cell, idx) => {
      const tag = i === 0 ? "th" : "td";
      const td  = document.createElement(tag);
      if (i === 0) {
        // header: clickable to sort
        const raw = row[idx];
        td.innerHTML = `${cmap[raw]||raw} ⬍`;
        td.style.cursor = "pointer";
        td.addEventListener("click", () => sortTableByColumn(idx));
      } else if (rows[0][idx] === "DRUG_NAME") {
        // make drug links to Wikipedia
        const a = document.createElement("a");
        a.href        = `https://en.wikipedia.org/wiki/${encodeURIComponent(cell)}`;
        a.textContent = cell;
        a.target      = "_blank";
        a.style.textDecoration = "underline";
        td.appendChild(a);
      } else {
        td.textContent = cell;
      }
      tr.appendChild(td);
    });
    (i === 0 ? thead : tbody).appendChild(tr);
  });

  table.appendChild(thead);
  table.appendChild(tbody);
  preview.appendChild(table);
}

// ─── 6. Sorting ────────────────────────────────────────────────────────────────

function sortTableByColumn(colIdx) {
  if (!allRows.length) return;
  const isNum = !isNaN(allRows[1][colIdx]);
  sortAscending = (currentSortColumn === colIdx) ? !sortAscending : true;
  currentSortColumn = colIdx;

  const sorted = allRows.slice(1).sort((a,b) => {
    const A = a[colIdx], B = b[colIdx];
    if (isNum) {
      return sortAscending ? A - B : B - A;
    }
    return sortAscending
      ? A.localeCompare(B)
      : B.localeCompare(A);
  });

  displayTable([allRows[0], ...sorted]);
}

// ─── 7. Filtering ─────────────────────────────────────────────────────────────

function createCancerTypeCheckboxes(rows) {
  const cancerIdx = rows[0].indexOf("CANCER_TYPE");
  if (cancerIdx < 0) return;

  const types = new Set(rows.slice(1).map(r => r[cancerIdx].trim()));
  const container = document.getElementById("cancer_filter_checkboxes");
  container.innerHTML = Array.from(types).sort().map(type => {
    const id = `chk_${type.replace(/\s+/g,"_")}`;
    return `
      <label>
        <input type="checkbox" id="${id}" value="${type}" checked
               onchange="applyCancerFilter()">
        ${type}
      </label><br>
    `;
  }).join("");

  document.getElementById("filter_controls").style.display = "block";
  applyCancerFilter();
}

function applyCancerFilter() {
  const checked = Array.from(
    document.querySelectorAll("#cancer_filter_checkboxes input:checked")
  ).map(i => i.value);

  const ci = allRows[0].indexOf("CANCER_TYPE");
  const filtered = [allRows[0], ...allRows.slice(1)
    .filter(r => checked.includes(r[ci].trim()))
  ];

  displayTable(filtered);
  document.getElementById("download_filtered")
          .style.display = (filtered.length>1 ? "inline-block" : "none");

  if (checked.length) displayTopDrugs(filtered);
  else document.getElementById("top_drugs_section").style.display = "none";
}

// ─── 8. Download Filtered ────────────────────────────────────────────────────

function downloadFilteredData() {
  const rows = Array.from(
    document.querySelectorAll(".styled-table tr")
  ).map(tr =>
    Array.from(tr.cells).map(td => `"${td.textContent}"`).join(",")
  );
  const blob = new Blob([rows.join("\n")], {type:"text/csv"});
  const a = document.createElement("a");
  a.href     = URL.createObjectURL(blob);
  a.download = "filtered_cancer_data.csv";
  a.click();
}

// ─── 9. Top-Drugs Charts ──────────────────────────────────────────────────────

function displayTopDrugs(rows) {
  const hdr = rows[0];
  const ci  = hdr.indexOf("CANCER_TYPE");
  const di  = hdr.indexOf("DRUG_NAME");
  const li  = hdr.indexOf("Predicted_LN_IC50");
  if (ci<0||di<0||li<0) return;

  // group & sort
  const map = {};
  rows.slice(1).forEach(r => {
    const c = r[ci].toUpperCase();
    const ln = parseFloat(r[li]);
    if (isNaN(ln)) return;
    if (!map[c]) map[c]=[];
    map[c].push({drug:r[di], ln});
  });

  const container = document.getElementById("top_drugs_container");
  container.innerHTML = "";
  const types = Object.keys(map);
  types.forEach((cancer, idx) => {
    const top10 = map[cancer]
      .sort((a,b)=>a.ln-b.ln)
      .slice(0,10);

    const wrapper = document.createElement("div");
    wrapper.style.cssText = `
      display:inline-block;
      vertical-align:top;
      width:${types.length>1?"48%":"96%"};
      margin:1%;
    `;

    const canvas = document.createElement("canvas");
    canvas.id = `chart_${idx}`;
    canvas.height = 250;
    wrapper.appendChild(canvas);
    container.appendChild(wrapper);

    new Chart(canvas.getContext("2d"), {
      type: "bar",
      data: {
        labels: top10.map(o=>o.drug),
        datasets: [{
          label: `Top 10 Sensitive (${cancer})`,
          data:  top10.map(o=>o.ln),
          backgroundColor: "#5d3b63"
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: { display: true, text: cancer, font:{size:18} }
        },
        scales: {
          y: { title: { display:true, text:"Pred LN IC50" } }
        }
      }
    });
  });

  document.getElementById("top_drugs_section").style.display = "block";
}

// ─── 10. Navigation ───────────────────────────────────────────────────────────

function showSection(section) {
  ["intro_section","about_section","upload-section"].forEach(id => {
    const el = document.getElementById(id);
    el.style.display = (id===section? "block":"none");
  });
  if (section==="upload") {
    document.getElementById("filter_controls").style.display="none";
  }
}
