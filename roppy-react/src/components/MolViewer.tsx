declare global {
  interface Window {
    $3Dmol: any;
  }
}

export function create3DmolViewer(container: HTMLElement) {
  return new window.$3Dmol.GLViewer(container, {
    defaultcolors: window.$3Dmol.rasmolElementColors,
    backgroundColor: 'white',
  });
}

export default function MolViewer({ smiles }: { smiles: string }) {
  return (
    <div id="molviewer" style={{ width: '100%', height: '400px' }} />
  );
} 