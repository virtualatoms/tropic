import { Modal, Button, Group, Space } from "@mantine/core";
import { useDisclosure } from "@mantine/hooks";

type Props = {
  opened: boolean;
  onClose: () => void;
  onApply: () => void;
  // You can pass in drawing data / canvas refs here later
};

export default function DrawModal({ opened, onClose, onApply }: Props) {
  return (
    <Modal
      opened={opened}
      onClose={onClose}
      title="Draw Molecule"
      size={650}
      centered
    >
      <div
        style={{
          height: 300,
          border: "1px solid #ccc",
          borderRadius: 8,
          marginBottom: 20,
          display: "flex",
          justifyContent: "center",
          alignItems: "center",
          fontStyle: "italic",
          color: "#666",
        }}
      >
      </div>

      <Space h={20} />
      <Group justify="flex-end">
        <Button variant="default" color="red" onClick={onClose}>
          Cancel
        </Button>
        <Button onClick={onApply}>Apply</Button>
      </Group>
    </Modal>
  );
}