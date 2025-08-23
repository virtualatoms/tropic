import { MonomerSummary } from "@/lib/types";

export const monomerSummaryToCSV = (data: MonomerSummary[]): string => {
  const header = [
    "Monomer ID",
    "SMILES",
    "Ring Size",
    "Has Exp Data",
    "Has Comp Data",
    "Is Experimental Point",
    "ΔH (kJ/mol)",
    "ΔS (J/mol·K)",
    "Ceiling Temperature (°C)",
    "DOI",
    "Reference",
  ];

  const rows: string[] = [];

  for (const summary of data) {
    for (const reaction of summary.data) {
      const row = [
        summary.monomer.monomer_id,
        summary.monomer.smiles,
        summary.monomer.ring_size,
        summary.has_exp,
        summary.has_comp,
        reaction.is_experimental,
        reaction.delta_h,
        reaction.delta_s,
        reaction.ceiling_temperature,
        reaction.doi,
        reaction.formatted_reference,
      ]
        .map((value) => {
          const str =
            value === null || value === undefined ? "" : String(value);
          return `"${str.replace(/"/g, '""')}"`; // Handle commas and quotes
        })
        .join(",");
      rows.push(row);
    }
  }

  return [header.join(","), ...rows].join("\n");
};

export const triggerDownload = (
  content: string,
  fileName: string,
  mimeType: string,
) => {
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = fileName;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
};
