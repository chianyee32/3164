let allRows = [];
let currentSortColumn = null;
let sortAscending = true;

document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("file_upload_prediction").addEventListener("change", function (event) {
    uploadFile(event, "file_info_prediction", "http://127.0.0.1:5000/upload_prediction");
  });

  document.getElementById("download_filtered").addEventListener("click", downloadFilteredData);
});

function uploadFile(event, fileInfoId, uploadUrl) {
  const file = event.target.files[0];
  if (!file || (!file.name.endsWith(".csv") && file.type !== "text/csv")) {
    alert("⚠️ Please upload a valid CSV file.");
    return;
  }

  const formData = new FormData();
  formData.append("file", file);

  const infoDiv = document.getElementById(fileInfoId);
  const loadingDiv = document.getElementById("loading_indicator");
  const uploadSection = document.getElementById("upload-section");
  const tooltipBox = document.getElementById("tooltip-box");
  const formatBtn = document.getElementById("format_btn");
  const spinner = document.getElementById("spinner");
  const preview = document.getElementById("csv_preview");

  uploadSection.style.display = "none";
  tooltipBox.style.display = "none";
  formatBtn.style.display = "none";
  spinner.style.display = "block";
  preview.innerHTML = "";
  infoDiv.innerHTML = `📄 File Selected: <span style="color: teal;">${file.name}</span>`;

  loadingDiv.innerHTML = `⏳ <span style="color: #9b59b6;">Processing your file, please wait...</span>`;
  loadingDiv.style.display = "block";

  fetch(uploadUrl, { method: "POST", body: formData })
    .then(res => res.json())
    .then(data => {
      spinner.style.display = "none";
      loadingDiv.style.display = "none";

      if (data.error) {
        let errorMessage = `❌ Error: ${data.error}`;
        if (data.details) {
          errorMessage += `\n\n🔍 Details:\n${data.details}`;
        }
        infoDiv.innerText += "\n" + errorMessage;
        uploadSection.style.display = "block";
        return;
      }

      const link = document.createElement("a");
      link.href = data.download_url;
      link.textContent = "⬇️ Download Prediction File";
      link.target = "_blank";
      link.style.cssText = "display:inline-block;margin-top:15px;background-color:#2ecc71;color:white;padding:10px 20px;border-radius:8px;text-decoration:none;";
      infoDiv.appendChild(document.createElement("br"));
      infoDiv.appendChild(link);

      fetch(data.download_url)
        .then(res => res.text())
        .then(csvText => {
          displayCSV(csvText);
          createCancerTypeCheckboxes(allRows);
          document.getElementById("csv_preview").scrollIntoView({ behavior: "smooth" });
        });
    })
    .catch(err => {
      spinner.style.display = "none";
      loadingDiv.style.display = "none";
      infoDiv.innerText = "❌ Upload failed.";
      uploadSection.style.display = "block";
    });
}

function displayCSV(csvText) {
  allRows = csvText.trim().split("\n").map(row => row.split(","));
  displayTable(allRows);
}

function displayTable(rows) {
  const tooltip = document.getElementById("tooltip-box");
  if (tooltip) tooltip.style.display = "none";

  const preview = document.getElementById("csv_preview");
  preview.innerHTML = "";

  const table = document.createElement("table");
  table.className = "styled-table";

  const thead = document.createElement("thead");
  const tbody = document.createElement("tbody");

  const columnNameMap = {
    "DRUG_ID": "Drug ID",
    "DRUG_NAME": "Drug Name",
    "COSMIC_ID": "Cosmic ID",
    "CCLE_Name": "Cell Line Name",
    "CCLE_Name": "Cell Line Name",
    "CANCER_TYPE": "Cancer Type",
    "Predicted_LN_IC50": "Predicted LN IC50"
  };

  const headerRow = rows[0];
  const includedIndices = headerRow.map((_, idx) => idx);

  rows.forEach((row, i) => {
    const tr = document.createElement("tr");
    includedIndices.forEach(colIndex => {
      const cell = row[colIndex];
      const tag = i === 0 ? "th" : "td";
      const td = document.createElement(tag);
      const header = headerRow[colIndex];

      if (i === 0) {
        td.innerHTML = (columnNameMap[header] || header) + " ⬍";
        td.style.cursor = "pointer";
        td.addEventListener("click", () => sortTableByColumn(colIndex));
      } else if (header === "DRUG_NAME") {
        const link = document.createElement("a");
        link.href = `https://en.wikipedia.org/wiki/${encodeURIComponent(cell)}`;
        link.textContent = cell;
        link.style.color = "#2c3e50";
        link.target = "_blank";
        link.style.textDecoration = "underline";
        td.appendChild(link);
      } else {
        td.textContent = cell;
      }
      tr.appendChild(td);
    });
    i === 0 ? thead.appendChild(tr) : tbody.appendChild(tr);
  });

  table.appendChild(thead);
  table.appendChild(tbody);
  preview.appendChild(table);
}

function sortTableByColumn(index) {
  if (!allRows.length) return;

  const isNumeric = !isNaN(allRows[1][index]);
  sortAscending = currentSortColumn === index ? !sortAscending : true;
  currentSortColumn = index;

  const sorted = [...allRows.slice(1)].sort((a, b) => {
    const aVal = a[index], bVal = b[index];
    return isNumeric
      ? (sortAscending ? aVal - bVal : bVal - aVal)
      : (sortAscending ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal));
  });

  displayTable([allRows[0], ...sorted]);
}

function createCancerTypeCheckboxes(dataRows) {
  const cancerSet = new Set();
  const cancerDiv = document.getElementById("cancer_filter_checkboxes");
  const filterControls = document.getElementById("filter_controls");

  const cancerIdx = dataRows[0].indexOf("CANCER_TYPE");
  if (cancerIdx === -1) return;

  dataRows.slice(1).forEach(row => {
    const type = row[cancerIdx];
    if (type) cancerSet.add(type.trim());
  });

  cancerDiv.innerHTML = "";
  const checkboxHTML = Array.from(cancerSet).sort().map(type => {
    const id = `chk_${type.replace(/\s+/g, "_")}`;
    return `<label><input type="checkbox" value="${type}" id="${id}" checked onchange="applyCancerFilter()"> ${type}</label><br>`;
  }).join('');
  
  cancerDiv.innerHTML = checkboxHTML;

  filterControls.style.display = "block";

  // ✅ Ensure checkboxes exist before calling this
  setTimeout(applyCancerFilter, 0);
}


function applyCancerFilter() {
  const selected = Array.from(document.querySelectorAll("#cancer_filter_checkboxes input:checked"))
                        .map(chk => chk.value.trim());

  const index = allRows[0].indexOf("CANCER_TYPE");
  if (index === -1) return;

  const filtered = [allRows[0], ...allRows.slice(1).filter(row => {
    const type = row[index]?.trim();
    return selected.includes(type);
  })];

  displayTable(filtered);

  const downloadBtn = document.getElementById("download_filtered");
  downloadBtn.style.display = filtered.length > 1 ? "inline-block" : "none";

  if (selected.length > 0) {
    displayTopDrugs(filtered);
  } else {
    document.getElementById("top_drugs_section").style.display = "none";
  }
}



function downloadFilteredData() {
  const table = document.querySelector(".styled-table");
  const rows = table.querySelectorAll("tr");

  const csv = Array.from(rows)
    .map(row => Array.from(row.cells).map(cell => `"${cell.textContent}"`).join(","))
    .join("\n");

  const blob = new Blob([csv], { type: "text/csv;charset=utf-8;" });
  const url = URL.createObjectURL(blob);

  const link = document.createElement("a");
  link.href = url;
  link.download = `filtered_cancer_data.csv`;
  link.click();
}

function showRequirements() {
  const tooltip = document.getElementById("tooltip-box");
  tooltip.style.display = tooltip.style.display === "block" ? "none" : "block";
}

function displayTopDrugs(dataRows) {
  
  const header = dataRows[0];
  const rows = dataRows.slice(1);
  const cancerIdx = header.indexOf("CANCER_TYPE");
  const drugNameIdx = header.indexOf("DRUG_NAME");
  const lnIC50Idx = header.indexOf("Predicted_LN_IC50");

  if (cancerIdx === -1 || drugNameIdx === -1 || lnIC50Idx === -1) return;

  const cancerDrugMap = {};

  rows.forEach(row => {
    const cancerType = row[cancerIdx]?.toUpperCase();
    const drugName = row[drugNameIdx];
    const lnIC50 = parseFloat(row[lnIC50Idx]);

    if (!cancerType || isNaN(lnIC50)) return;

    if (!cancerDrugMap[cancerType]) cancerDrugMap[cancerType] = [];
    cancerDrugMap[cancerType].push({ drugName, lnIC50 });
  });

  const container = document.getElementById("top_drugs_container");
  container.innerHTML = "";

  const chartCount = Object.keys(cancerDrugMap).length;

  Object.entries(cancerDrugMap).forEach(([cancerType, drugs], index) => {
    const sortedDrugs = drugs.sort((a, b) => a.lnIC50 - b.lnIC50).slice(0, 10);

    const wrapper = document.createElement("div");
    wrapper.style.display = "inline-block";
    wrapper.style.verticalAlign = "top";
    wrapper.style.width = chartCount === 1 ? "96%" : "48%";
    wrapper.style.margin = "1%";

    const canvas = document.createElement("canvas");
    canvas.id = `chart_${cancerType}_${index}`;
    canvas.height = 250;
    wrapper.appendChild(canvas);
    container.appendChild(wrapper);

    const ctx = canvas.getContext("2d");
    new Chart(ctx, {
      type: 'bar',
      data: {
        labels: sortedDrugs.map(d => d.drugName),
        datasets: [{
          label: `Top 10 Sensitive Drugs (${cancerType})`,
          data: sortedDrugs.map(d => d.lnIC50),
          backgroundColor: '#5d3b63'
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: `Cancer Type: ${cancerType}`,
            font: { size: 18 }
          }
        },
        scales: {
          y: {
            title: {
              display: true,
              text: 'Predicted LN IC50 (lower = more sensitive)'
            }
          }
        }
      }
    });
  });

  document.getElementById("top_drugs_section").style.display = "block";
}


function showSection(section) {
  const sections = ['intro_section', 'about_section', 'upload-section'];
  sections.forEach(id => {
    const el = document.getElementById(id);
    if (el) el.style.display = (id.includes(section)) ? 'block' : 'none';
  });

  if (section === 'upload') {
    document.getElementById('filter_controls').style.display = 'none';  // Optional reset
  }
}




