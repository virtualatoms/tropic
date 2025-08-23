import { MonomerFilters, MonomerSummary, Reaction } from "./types";
import { API_ENDPOINT } from "./constants";

const buildMonomerQuery = (filters: MonomerFilters): string => {
  const params = new URLSearchParams();

  if (filters.smiles) {
    params.append("search", filters.smiles);
  }

  params.append("monomer__ring_size__gte", filters.ringSize[0].toString());
  if (filters.ringSize[1] < 15) {
    params.append("monomer__ring_size__lte", filters.ringSize[1].toString());
  }

  params.append(
    "monomer__molecular_weight__gte",
    filters.molWeight[0].toString(),
  );
  params.append(
    "monomer__molecular_weight__lte",
    filters.molWeight[1].toString(),
  );

  if (filters.hasComp !== "both") {
    params.append("has_comp", (filters.hasComp === "yes").toString());
  }
  if (filters.hasExp !== "both") {
    params.append("has_exp", (filters.hasExp === "yes").toString());
  }

  // This assumes your API can handle multiple functional_groups params
  filters.functionalGroups.forEach((group) => {
    params.append("functional_groups", group);
  });

  params.append("size", "1000");
  return params.toString();
};

export async function fetchMonomerSummaries(
  filters: MonomerFilters,
): Promise<MonomerSummary[]> {
  const query = buildMonomerQuery(filters);
  const response = await fetch(`${API_ENDPOINT}/monomer-summaries?${query}`);
  if (!response.ok) {
    throw new Error("Network response was not ok");
  }
  return (await response.json()) as MonomerSummary[];
}

export const fetchReactions = async (
  filters: MonomerFilters,
): Promise<Reaction[]> => {
  const query = buildMonomerQuery(filters);
  const response = await fetch(`${API_ENDPOINT}/reactions?${query}`);
  if (!response.ok) {
    throw new Error("Failed to fetch data for export");
  }
  return (await response.json()) as Reaction[];
};
