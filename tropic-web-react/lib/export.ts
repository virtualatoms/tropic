import { MonomerSummary, Reaction } from "@/lib/types";

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


export const reactionToCSV = (data: Reaction[]): string => {
  const header = [
    "Reaction ID",
    "Reaction Type",
    "Monomer ID",
    "Monomer SMILES",
    "Ring Size",
    "Repeating Units",
    "Is Experimental",
    "Temperature (K)",
    "Pressure (atm)",
    "State Summary",
    "Initiator SMILES",
    "Medium",
    "Solvent",
    "Method",
    "ΔH (kJ/mol)",
    "ΔS (J/mol·K)",
    "ΔH Uncertainty",
    "ΔS Uncertainty",
    "Ceiling Temperature (°C)",
    "DOI",
    "Year",
    "Reference",
    "Comment",
  ];

  const rows = data.map((r) => {
    const rowData = [
      r.reaction_id,
      r.type,
      r.monomer.monomer_id,
      r.monomer.smiles,
      r.monomer.ring_size,
      r.product.repeating_units,
      r.parameters.is_experimental,
      r.parameters.temperature,
      r.parameters.pressure,
      r.parameters.state_summary,
      r.parameters.initiator_smiles,
      r.parameters.medium,
      r.parameters.solvent,
      r.parameters.method,
      r.thermo.delta_h,
      r.thermo.delta_s,
      r.thermo.delta_h_std,
      r.thermo.delta_s_std,
      r.thermo.ceiling_temperature,
      r.metadata.doi,
      r.metadata.year,
      r.metadata.formatted_reference,
      r.metadata.comment,
    ];

    return rowData
      .map((value) => {
        const str =
          value === null || value === undefined ? "" : String(value);
        // Enclose in quotes and escape existing quotes
        return `"${str.replace(/"/g, '""')}"`;
      })
      .join(",");
  });

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
