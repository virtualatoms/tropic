import Link from "next/link";
import { Text, Table, Anchor, Image } from "@mantine/core";

export const ReactionIdCell = ({ value }: { value: string }) => (
  <Anchor href={`/reactions/${value}`} size="sm" component={Link}>
    {value}
  </Anchor>
);

export const MonomerIdCell = ({ value }: { value: string }) => (
  <Anchor href={`/monomers/${value}`} size="sm" component={Link}>
    {value}
  </Anchor>
);

export const HtmlHeader = ({ displayName }: { displayName: string }) => (
  <span dangerouslySetInnerHTML={{ __html: displayName }} />
);

export const ImageRenderer = ({ value }: { value: string }) => (
  <div
    style={{
      display: "flex",
      alignItems: "center",
      justifyContent: "center",
      width: "100%",
      height: "100%",
    }}
  >
    <Image src={`data:image/svg+xml;base64,${value}`} alt="Monomer structure" style={{ maxHeight: "100%", maxWidth: "100%", objectFit: "contain" }} />
  </div>
);

export const RefCell = ({ value }: { value: string | null }) => (
  <Anchor href={`https://doi.org/${value}`} target="_blank">
    <svg
      xmlns="http://www.w3.org/2000/svg"
      width="14"
      height="14"
      viewBox="0 0 24 24"
      style={{ display: "inline" }}
    >
      <path
        fill="none"
        stroke="currentColor"
        strokeLinecap="round"
        strokeLinejoin="round"
        strokeWidth="3"
        d="M13.5 10.5L21 3m-5 0h5v5m0 6v5a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h5"
      />
    </svg>
  </Anchor>
);

export function boolFormatter(params: { value: boolean }): string {
  return params.value ? "Yes" : "No";
}

export function decimalFormatter(value: number | null): string {
  return value !== null ? value.toFixed(2) : "";
}

export const idComparator = (valueA: string, valueB: string) => {
  const numA = parseInt(valueA.split("-")[1], 10);
  const numB = parseInt(valueB.split("-")[1], 10);
  return numA - numB;
};
