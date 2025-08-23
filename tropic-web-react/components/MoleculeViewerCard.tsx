import { useEffect, useRef } from "react";
import { Card } from "@mantine/core";

interface MoleculeViewerCardProps {
  xyz: string;
}

export function MoleculeViewerCard({ xyz }: MoleculeViewerCardProps) {
  const viewerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!viewerRef.current || !xyz) return;

    // Dynamically import 3dmol on the client
    import("3dmol").then(($3Dmol) => {
      viewerRef.current!.innerHTML = "";

      const config = { backgroundColor: "white" };
      const viewer = $3Dmol.createViewer(viewerRef.current!, config);

      viewer.addModel(xyz, "xyz");
      viewer.setStyle({}, { stick: { colorscheme: "Jmol" } });
      viewer.zoomTo();
      viewer.render();
      viewer.zoom(1.2, 0);
    });
  }, [xyz]);

  return (
    <Card withBorder shadow="sm" radius="md" p={0}>
      <div
        ref={viewerRef}
        style={{
          width: "100%",
          aspectRatio: "1 / 1",
          position: "relative",
        }}
      />
    </Card>
  );
}
