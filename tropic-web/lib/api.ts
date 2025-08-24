import {
  MonomerFilters,
  MonomerSummary,
  Reaction,
  ReactionFilters,
} from "./types";
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

  if (filters.functionalGroups.length > 0) {
    params.append(
      "monomer__functional_group__in",
      filters.functionalGroups.join(","),
    );
  }

  return params.toString();
};

export const buildReactionQuery = (filters: ReactionFilters): string => {
  const params = new URLSearchParams();

  if (filters.smiles) {
    params.append("search", filters.smiles);
  }

  if (filters.initiator) {
    params.append("parameters__initiator_smiles", filters.initiator);
  }

  if (filters.isExperimental !== "both") {
    params.append(
      "parameters__is_experimental",
      (filters.isExperimental === "yes").toString(),
    );
  }

  params.append("metadata__year__gte", filters.year[0].toString());
  params.append("metadata__year__lte", filters.year[1].toString());

  if (filters.medium.length > 0) {
    params.append("parameters__medium__in", filters.medium.join(","));
  }
  if (filters.reactionType.length > 0) {
    params.append("type__in", filters.reactionType.join(","));
  }
  if (filters.method.length > 0) {
    params.append("parameters__method__in", filters.method.join(","));
  }

  params.append("monomer__ring_size__gte", filters.ringSize[0].toString());
  if (filters.ringSize[1] < 15) {
    params.append("monomer__ring_size__lte", filters.ringSize[1].toString());
  }

  params.append(
    "monomer__molecular_weight__gte",
    filters.molecularWeight[0].toString(),
  );
  params.append(
    "monomer__molecular_weight__lte",
    filters.molecularWeight[1].toString(),
  );

  params.append("include_svg", "true");

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
  filters: ReactionFilters,
): Promise<Reaction[]> => {
  const query = buildReactionQuery(filters);
  const response = await fetch(`${API_ENDPOINT}/reactions?${query}`);
  if (!response.ok) {
    throw new Error("Failed to fetch data for export");
  }
  return (await response.json()) as Reaction[];
};
